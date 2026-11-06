from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from gene_lab.config import load_config, read_genes
from gene_lab.expression import summarize_expression
from gene_lab.storage import Store
from gene_lab.synthetic import make_expression, make_vcf


def test_hand_checked_toy_expression_and_missing_gene(tmp_path, store):
    path = tmp_path / "tiny.csv"
    path.write_text("gene,C1,C2,C3\nA,0,2,4\nB,0,0,3\n")
    rows = summarize_expression(path, ["A", "B", "MISSING"], store)
    assert rows[0]["mean"] == 2 and rows[0]["median"] == 2
    assert rows[0]["detected_fraction"] == pytest.approx(2 / 3)
    assert rows[1]["mean"] == 1
    assert rows[2]["status"] == "missing" and np.isnan(rows[2]["mean"])
    assert (store.directory / "expression.png").is_file()


@pytest.mark.parametrize("text,message", [("gene,C1,C1\nA,1,2\n", "unique cell"), ("gene,C1\nA,1\nA,2\n", "duplicate genes"),
    ("gene,C1\nA,-1\n", "nonnegative"), ("gene,C1\nA,word\n", "numeric"), ("gene,C1\nA,\n", "finite"), ("gene,C1\n", "empty")])
def test_bad_expression_matrices(tmp_path, store, text, message):
    path = tmp_path / "bad.csv"
    path.write_text(text)
    with pytest.raises(ValueError, match=message):
        summarize_expression(path, ["A"], store)


def test_generated_toy_counts_only_synthetic(tmp_path):
    path = tmp_path / "expression.csv"
    make_expression(["A", "B"], path)
    frame = pd.read_csv(path, index_col=0)
    assert frame.shape == (2, 24)
    assert all(name.startswith("CELL") for name in frame.columns)
    assert (frame.to_numpy() >= 0).all()


def test_gene_list_deduplicates_and_checks_symbols(tmp_path):
    path = tmp_path / "genes.txt"
    path.write_text("# synthetic test\ntp53\nTP53\nBRCA1\n\n")
    assert read_genes(path) == ["TP53", "BRCA1"]
    path.write_text("../../bad\n")
    with pytest.raises(ValueError, match="invalid gene"):
        read_genes(path)
    path.write_text("\n#empty\n")
    with pytest.raises(ValueError, match="empty"):
        read_genes(path)


@pytest.mark.parametrize("key,value", [("assembly", "GRCh37"), ("samples", 0), ("attempts", 11), ("alpha", 1), ("planted_index", 301), ("backoff", -1)])
def test_bad_yaml_values(tmp_path, key, value):
    raw = yaml.safe_load(Path("config.yaml").read_text())
    raw[key] = value
    path = tmp_path / "bad.yaml"
    path.write_text(yaml.safe_dump(raw))
    with pytest.raises(ValueError):
        load_config(path)


def test_yaml_mapping_keys_and_relative_paths(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text("- a\n")
    with pytest.raises(ValueError, match="mapping"):
        load_config(path)
    path.write_text("assembly: GRCh38\n")
    with pytest.raises(ValueError, match="config keys"):
        load_config(path)
    config = load_config("config.yaml")
    assert config["recordings"].is_absolute()


def test_storage_csv_sqlite_metadata_and_unknown_stage(store):
    store.write("vcf_issues", [{"line": 1, "message": "synthetic"}])
    assert store.read("vcf_issues").iloc[0]["message"] == "synthetic"
    assert (store.directory / "vcf_issues.csv").is_file()
    store.metadata("example", {"synthetic": True})
    with pytest.raises(ValueError, match="Unknown"):
        store.write("bad", [])
    with pytest.raises(ValueError, match="Unknown"):
        store.read("bad")
    with pytest.raises(ValueError, match="Missing pipeline stage"):
        store.read("association")


def test_reference_length_guard_and_local_generation(tmp_path):
    gene = {"symbol": "SYN", "chrom": "1", "start": 1, "end": 100}
    with pytest.raises(ValueError, match="length"):
        make_vcf([gene], {"SYN": "A"}, tmp_path / "out.vcf")
    rows = make_vcf([gene], {"SYN": "ACGT" * 25}, tmp_path / "out.vcf")
    assert len(rows) == 20 and rows[0].identifier == "SYN0001"


def test_store_reopens_sqlite(tmp_path):
    one = Store(tmp_path / "db")
    one.write("vcf_issues", [])
    one.close()
    two = Store(tmp_path / "db")
    assert two.read("vcf_issues").empty
    two.close()
