import os

import pytest

from gene_lab.http import HttpClient


@pytest.mark.live
@pytest.mark.skipif(os.getenv("RUN_LIVE") != "1", reason="Live public API test disabled by default; set RUN_LIVE=1 explicitly")
def test_live_ensembl_lookup(config):
    client = HttpClient(config)
    try:
        result = client.json("GET", config["ensembl_url"] + "/lookup/symbol/homo_sapiens/TP53", {"expand": 1})
        assert result["assembly_name"] == "GRCh38" and result["Transcript"]
    finally:
        client.close()
