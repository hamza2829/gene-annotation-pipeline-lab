import json
from pathlib import Path

import pytest
import yaml

from gene_lab import cli
from gene_lab.config import read_genes
from gene_lab.http import ApiError, HttpClient, RecordedSession
from gene_lab.pipeline import annotate_stage, associate_stage, fetch_stage, report_stage, verify_outputs
from gene_lab.providers import fetch_genes


def test_recorded_gene_protein_rna_and_transcript_fetch(config, store):
    client = HttpClient(config)
    genes = fetch_stage(read_genes("data/genes.txt"), client, config, store)
    assert len(genes) == 10
    assert len(store.read("transcripts")) > 10
    assert len(store.read("proteins")) > 0 and len(store.read("ncrna")) > 0
    assert set(store.read("provider_status")["provider"]) == {"Ensembl", "UniProt", "RNAcentral"}
    assert store.read("genes")["assembly"].eq("GRCh38").all()
    client.close()


def test_full_chain_on_genuine_recordings_and_synthetic_data(config, store):
    client = HttpClient(config, RecordedSession(config["recordings"]))
    symbols = read_genes("data/genes.txt")
    fetch_stage(symbols, client, config, store)
    annotated = annotate_stage(client, config, store)
    assert len(annotated) == 200 and any(row["vep_status"] == "ok" for row in annotated)
    assert all(row["overlapping_genes"] for row in annotated)
    associate_stage(config, store)
    path = report_stage(symbols, config, store)
    text = path.read_text()
    assert "data:image/png;base64" in text and "Toy single-cell" in text
    assert "SIM0043" in text and "not clinical classifications" in text
    assert "<script" not in text
    proof = verify_outputs(store.directory)
    assert proof["counts"]["variants"] == 200 and proof["association"]["planted_detected"]
    client.close()


def test_provider_failures_are_not_empty_successes(config, store):
    class Offline:
        def json(self, *args, **kwargs):
            raise ApiError("synthetic outage")
    assert fetch_genes(["UNKNOWN"], Offline(), config, store) == []
    statuses = store.read("provider_status")
    assert list(statuses["status"]) == ["error", "error", "unavailable"]


def test_provider_shape_and_assembly_errors(config, store):
    class Broken:
        def json(self, method, url, params=None):
            return {}
    fetch_genes(["SYN"], Broken(), config, store)
    assert "required gene coordinates" in store.read("provider_status").iloc[0]["message"]
    class WrongBuild:
        def json(self, method, url, params=None):
            if "lookup" in url:
                return {"id": "x", "seq_region_name": "1", "start": 1, "end": 2, "assembly_name": "GRCh37"}
            return {"results": []}
    fetch_genes(["SYN"], WrongBuild(), config, store)
    assert "Assembly mismatch" in store.read("provider_status").iloc[0]["message"]


def test_empty_protein_and_rna_results_and_bad_rna_shape(config, store, record):
    gene = record("/lookup/symbol/homo_sapiens/BRCA1")["body"]
    class Empty:
        def json(self, method, url, params=None):
            return gene if "lookup" in url else {"results": []}
    fetch_genes(["BRCA1"], Empty(), config, store)
    assert store.read("provider_status")["status"].tolist() == ["ok", "empty", "empty"]
    class BadRna:
        def json(self, method, url, params=None):
            return gene if "lookup" in url else {"results": []} if "uniprot" in url else {"results": "bad"}
    fetch_genes(["BRCA1"], BadRna(), config, store)
    assert store.read("provider_status").iloc[-1]["status"] == "error"


def test_failed_lookup_stops_all_with_finding_file(config, store):
    class Failed:
        offline = True
        def json(self, *args, **kwargs):
            raise ApiError("synthetic failure")
    with pytest.raises(ValueError, match="Some Ensembl"):
        fetch_stage(["UNKNOWN"], Failed(), config, store)
    assert (store.directory / "provider_status.csv").exists()


def config_file(tmp_path):
    raw = yaml.safe_load(Path("config.yaml").read_text())
    for key in ["recordings", "vcf", "expression"]:
        raw[key] = str(Path(raw[key]).resolve())
    raw.update(output=str(tmp_path / "cli-run"), cache=str(tmp_path / "cache.sqlite"), min_interval=0)
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(raw))
    return path


def test_cli_all_and_individual_commands(tmp_path):
    config = config_file(tmp_path)
    common = ["--config", str(config), "--genes", "data/genes.txt", "--offline"]
    assert cli.main(["all", *common]) == 0
    for command in ["fetch", "annotate", "associate", "report"]:
        assert cli.main([command, *common]) == 0
    summary = json.loads((tmp_path / "cli-run" / "association-summary.json").read_text())
    assert summary["samples"] == 2000 and summary["markers"] == 300


def test_cli_input_and_missing_stage_errors(tmp_path, capsys):
    assert cli.main(["all", "--config", str(tmp_path / "missing.yaml")]) == 1
    assert "ERROR:" in capsys.readouterr().err
    path = config_file(tmp_path)
    assert cli.main(["annotate", "--config", str(path), "--offline"]) == 1
    assert "Missing pipeline stage" in capsys.readouterr().err
    with pytest.raises(SystemExit) as error:
        cli.main([])
    assert error.value.code == 2


def test_recordings_have_provenance_and_are_not_invented():
    index = json.loads(Path("data/recordings/index.json").read_text())
    assert index["records"]
    for record in index["records"]:
        assert record["url"].startswith("https://") and len(record["wire_sha256"]) == 64 and record["retrieved_at"]
        assert (Path("data/recordings") / record["file"]).exists()
