import json
from pathlib import Path

import pytest
import responses
from responses import matchers

from gene_lab.config import load_config
from gene_lab.storage import Store


@pytest.fixture(autouse=True)
def no_network(request):
    if request.node.get_closest_marker("live"):
        yield None
        return
    with responses.RequestsMock(assert_all_requests_are_fired=False) as mock:
        for path in Path("data/recordings").glob("*.json"):
            if path.name == "index.json":
                continue
            record = json.loads(path.read_text())
            matching = []
            if record.get("params"):
                matching.append(matchers.query_param_matcher({key: str(value) for key, value in record["params"].items()}))
            if record.get("request_body") is not None:
                matching.append(matchers.json_params_matcher(record["request_body"]))
            mock.add(record["method"], record["url"], json=record["body"], status=record["status"], match=matching)
        yield mock


@pytest.fixture
def config(tmp_path):
    cfg = load_config("config.yaml")
    cfg.update(output=tmp_path / "run", cache=tmp_path / "cache.sqlite", min_interval=0, backoff=0, max_backoff=0.1)
    return cfg


@pytest.fixture
def store(config):
    database = Store(config["output"])
    yield database
    database.close()


@pytest.fixture
def record():
    def find(fragment):
        for path in Path("data/recordings").glob("*.json"):
            if path.name != "index.json":
                data = json.loads(path.read_text())
                if fragment in data["url"]:
                    return data
        raise AssertionError(f"Recorded public response missing: {fragment}")
    return find
