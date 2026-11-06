from datetime import datetime, timedelta, timezone
from email.utils import format_datetime
import json
from pathlib import Path
import sqlite3

import pytest
import requests

from gene_lab.http import ApiError, CaptureSession, HttpClient, RecordedSession, request_key

URL = "https://example.invalid/test"


def test_canonical_cache_key_query_order_and_method():
    assert request_key("GET", URL, {"a": 1, "b": 2}) == request_key("get", URL + "?b=2&a=1")
    assert request_key("GET", URL) != request_key("POST", URL, body={"a": 1})


def test_recorded_response_and_live_cache_without_network(config, record, no_network):
    saved = record("/lookup/symbol/homo_sapiens/BRCA1")
    client = HttpClient(config)
    first = client.json("GET", saved["url"], saved["params"])
    calls = len(no_network.calls)
    assert first["id"]
    assert client.json("GET", saved["url"], saved["params"]) == first
    assert len(no_network.calls) == calls
    client.close()


def test_expired_cache_refetches(config, no_network):
    no_network.get(URL, json={"value": 1})
    client = HttpClient(config)
    client.json("GET", URL)
    client.cache.execute("UPDATE http_cache SET saved=0")
    client.cache.commit()
    client.json("GET", URL)
    assert len(no_network.calls) == 2
    client.close()


def test_replay_never_uses_live_cache_namespace(config, record):
    saved = record("/lookup/symbol/homo_sapiens/BRCA1")
    live = HttpClient(config)
    live.json("GET", saved["url"], saved["params"])
    live.close()
    replay = HttpClient(config, RecordedSession(config["recordings"]))
    assert replay.json("GET", saved["url"], saved["params"])["id"]
    count = replay.cache.execute("SELECT count(*) FROM http_cache").fetchone()[0]
    assert count == 2
    replay.close()


@pytest.mark.parametrize("status", [429, 500, 503])
def test_retry_after_and_backoff(config, no_network, status):
    config.update(backoff=0.01, max_backoff=0.05)
    no_network.get(URL, status=status, headers={"Retry-After": "0.02"})
    no_network.get(URL, json={"ok": True})
    delays = []
    client = HttpClient(config, sleep=delays.append)
    assert client.json("GET", URL) == {"ok": True}
    assert 0.02 in delays and len(no_network.calls) == 2
    client.close()


def test_network_failure_retries_then_errors(config, no_network):
    no_network.get(URL, body=requests.Timeout("synthetic timeout"))
    client = HttpClient(config, sleep=lambda _: None)
    with pytest.raises(ApiError, match="failed after 3 attempts"):
        client.json("GET", URL)
    assert len(no_network.calls) == 3
    client.close()


@pytest.mark.parametrize("status", [400, 404, 403])
def test_permanent_http_error_does_not_retry(config, no_network, status):
    no_network.get(URL, status=status, json={"error": "synthetic failure"})
    client = HttpClient(config)
    with pytest.raises(ApiError, match=f"HTTP {status}"):
        client.json("GET", URL)
    assert len(no_network.calls) == 1
    client.close()


def test_invalid_json_has_url_and_no_cache(config, no_network):
    no_network.get(URL, body="not JSON", status=200)
    client = HttpClient(config)
    with pytest.raises(ApiError, match="invalid JSON"):
        client.json("GET", URL)
    assert client.cache.execute("SELECT count(*) FROM http_cache").fetchone()[0] == 0
    client.close()


def test_rate_limit_reset_paces_next_request(config, no_network):
    config.update(max_backoff=2, min_interval=0.25)
    no_network.get(URL, json={}, headers={"X-RateLimit-Remaining": "0", "X-RateLimit-Reset": "1"})
    no_network.get(URL + "/next", json={})
    delays = []
    client = HttpClient(config, sleep=delays.append, clock=lambda: 10)
    client.json("GET", URL)
    client.json("GET", URL + "/next")
    assert delays[-1] == 1.25
    client.close()


def test_retry_http_date_and_invalid_header(config):
    config.update(backoff=0.01, max_backoff=2)
    client = HttpClient(config)
    response = requests.Response()
    response.headers["Retry-After"] = format_datetime(datetime.now(timezone.utc) + timedelta(seconds=30))
    assert client.delay(response, 0) == 2
    response.headers["Retry-After"] = "invalid"
    assert client.delay(response, 2) == 0.04
    assert client.delay(None, 1) == 0.02
    client.close()


def test_capture_keeps_actual_response_and_hash_in_temp_directory(config, record, tmp_path):
    saved = record("/lookup/symbol/homo_sapiens/BRCA1")
    capture = CaptureSession(tmp_path / "captured")
    response = capture.request("GET", saved["url"], params=saved["params"])
    assert response.status_code == 200
    path = next((tmp_path / "captured").glob("*.json"))
    data = json.loads(path.read_text())
    assert data["body"] == saved["body"] and data["url"] == saved["url"]
    assert len(data["wire_sha256"]) == 64 and data["retrieved_at"]


def test_missing_recording_has_clear_error(tmp_path):
    with pytest.raises(ApiError, match="No recorded response"):
        RecordedSession(tmp_path).request("GET", URL)


def test_unit_test_socket_is_not_available(no_network):
    with pytest.raises(requests.ConnectionError):
        requests.get("https://unrecorded.invalid/no-network")


def test_cache_is_real_sqlite(config, no_network):
    no_network.get(URL, json={"ok": 1})
    client = HttpClient(config)
    client.json("GET", URL)
    client.close()
    with sqlite3.connect(config["cache"]) as connection:
        assert connection.execute("SELECT count(*) FROM http_cache").fetchone()[0] == 1
    assert Path(config["cache"]).is_file()
