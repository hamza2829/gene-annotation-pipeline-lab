from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
import hashlib
import json
from pathlib import Path
import sqlite3
import time
from urllib.parse import parse_qsl, urlencode, urlsplit, urlunsplit

import requests


class ApiError(ValueError):
    pass


def request_key(method, url, params=None, body=None):
    parts = urlsplit(url)
    query = parse_qsl(parts.query) + [(str(key), str(value)) for key, value in (params or {}).items()]
    canonical = urlunsplit((parts.scheme, parts.netloc, parts.path, urlencode(sorted(query)), ""))
    raw = json.dumps([method.upper(), canonical, body], sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(raw.encode()).hexdigest()


class RecordedSession:
    offline = True

    def __init__(self, directory):
        self.directory = Path(directory)

    def request(self, method, url, params=None, json=None, **kwargs):
        key = request_key(method, url, params, json)
        path = self.directory / f"{key}.json"
        if not path.is_file():
            raise ApiError(f"No recorded response for {method} {url}; use live mode to fetch new symbols")
        import json as codec
        record = codec.loads(path.read_text(encoding="utf-8"))
        response = requests.Response()
        response.status_code = record["status"]
        response.headers.update(record["headers"])
        response.url = record["url"]
        response._content = codec.dumps(record["body"]).encode()
        return response


class CaptureSession:
    offline = False

    def __init__(self, directory):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.session = requests.Session()

    def request(self, method, url, params=None, json=None, **kwargs):
        key = request_key(method, url, params, json)
        if (self.directory / f"{key}.json").is_file():
            print(f"Reuse captured response: {method} {url}", flush=True)
            return RecordedSession(self.directory).request(method, url, params, json)
        print(f"Capture public response: {method} {url}", flush=True)
        response = self.session.request(method, url, params=params, json=json, **kwargs)
        if response.status_code < 400:
            key = request_key(method, url, params, json)
            import json as codec
            record = {"method": method, "url": url, "params": params, "request_body": json, "status": response.status_code,
                "headers": dict(response.headers), "body": response.json(), "retrieved_at": datetime.now(timezone.utc).isoformat(),
                "wire_sha256": hashlib.sha256(response.content).hexdigest(), "source": "public API response, not synthetic annotation"}
            (self.directory / f"{key}.json").write_text(codec.dumps(record, indent=2) + "\n", encoding="utf-8")
        return response


class HttpClient:
    def __init__(self, config, session=None, sleep=time.sleep, clock=time.monotonic):
        self.config, self.session, self.sleep, self.clock = config, session or requests.Session(), sleep, clock
        self.offline = getattr(self.session, "offline", False)
        self.ready = 0.0
        cache = Path(config["cache"])
        cache.parent.mkdir(parents=True, exist_ok=True)
        self.cache = sqlite3.connect(cache)
        self.cache.execute("CREATE TABLE IF NOT EXISTS http_cache (key TEXT PRIMARY KEY, saved REAL, body TEXT)")

    def close(self):
        self.cache.close()

    def delay(self, response, attempt):
        fallback = self.config["backoff"] * 2 ** attempt
        value = response.headers.get("Retry-After") if response is not None else None
        if value:
            try:
                fallback = float(value)
            except ValueError:
                try:
                    fallback = (parsedate_to_datetime(value) - datetime.now(timezone.utc)).total_seconds()
                except (ValueError, TypeError):
                    pass
        return min(self.config["max_backoff"], max(0, fallback))

    def json(self, method, url, params=None, body=None):
        key = ("replay:" if self.offline else "live:") + request_key(method, url, params, body)
        cached = self.cache.execute("SELECT saved, body FROM http_cache WHERE key=?", (key,)).fetchone()
        if cached and time.time() - cached[0] < self.config["cache_ttl"]:
            return json.loads(cached[1])
        for attempt in range(int(self.config["attempts"])):
            if not self.offline:
                self.sleep(max(0, self.ready - self.clock()))
            response = None
            try:
                response = self.session.request(method, url, params=params, json=body,
                    headers={"Accept": "application/json", "Content-Type": "application/json", "User-Agent": "gene-annotation-pipeline-lab/0.1"}, timeout=self.config["timeout"])
                self.ready = self.clock() + self.config["min_interval"]
                if response.status_code == 429 or response.status_code >= 500:
                    raise requests.HTTPError(f"HTTP {response.status_code}")
                if response.status_code >= 400:
                    raise ApiError(f"{method} {url}: HTTP {response.status_code}: {response.text[:180]}")
                try:
                    data = response.json()
                except ValueError as error:
                    raise ApiError(f"{method} {url}: invalid JSON response") from error
                self.cache.execute("INSERT OR REPLACE INTO http_cache VALUES (?,?,?)", (key, time.time(), json.dumps(data)))
                self.cache.commit()
                if response.headers.get("X-RateLimit-Remaining") == "0":
                    try:
                        self.ready += min(float(response.headers.get("X-RateLimit-Reset", 0)), self.config["max_backoff"])
                    except ValueError:
                        pass
                return data
            except requests.RequestException as error:
                if attempt + 1 == self.config["attempts"]:
                    status = f" HTTP {response.status_code}" if response is not None else ""
                    raise ApiError(f"{method} {url}: failed after {attempt + 1} attempts{status} ({type(error).__name__})") from error
                self.sleep(self.delay(response, attempt))
        raise ApiError("No request attempts configured")
