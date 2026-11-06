import gzip

import pytest

from gene_lab.http import ApiError
from gene_lab.variants import Variant, annotate, overlaps, read_vcf

HEADER = "##fileformat=VCFv4.2\n##reference=GRCh38\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def test_real_demo_has_200_synthetic_variants_in_ten_regions():
    variants, issues = read_vcf("data/synthetic.vcf")
    assert len(variants) == 200 and issues == []
    assert {row.change for row in variants} == {"SNV", "insertion", "deletion"}
    assert all(row.identifier.startswith("SYN") for row in variants)


def test_multiallelic_split_and_chr_prefix(tmp_path):
    path = tmp_path / "split.vcf"
    path.write_text(HEADER + "chr7\t100\tSYN1\tA\tG,AT\t.\tPASS\t.\n")
    variants, issues = read_vcf(path)
    assert len(variants) == 2 and issues == []
    assert variants[0].key == "7:100:A:G"
    assert variants[0].change == "SNV" and variants[1].change == "insertion"
    assert variants[0].vep_input() == "7 100 SYN1 A G . . ."


@pytest.mark.parametrize("row", ["1\t0\tx\tA\tG\t.\tPASS\t.", "1\tbad\tx\tA\tG\t.\tPASS\t.", "1\t10\tx\tQ\tG\t.\tPASS\t.", "1\t10"])
def test_bad_rows_include_filename_and_line(tmp_path, row):
    path = tmp_path / "bad.vcf"
    path.write_text(HEADER + row + "\n")
    with pytest.raises(ValueError, match=r"bad.vcf:4"):
        read_vcf(path)


@pytest.mark.parametrize("alt", [".", "*", "<DEL>", "A"])
def test_unsupported_and_reference_equal_alleles_are_reported(tmp_path, alt):
    path = tmp_path / "unsupported.vcf"
    path.write_text(HEADER + f"1\t100\tx\tA\t{alt}\t.\tPASS\t.\n")
    variants, issues = read_vcf(path)
    assert variants == [] and issues[0]["line"] == 4


def test_missing_bad_header_assembly_and_duplicate(tmp_path):
    path = tmp_path / "bad.vcf"
    for text, message in [("1\t1\tx\tA\tC", "missing #CHROM"), ("", "missing #CHROM"),
            (HEADER.replace("GRCh38", "GRCh37"), "assembly"), ("#CHROM\tbad\n", "column header"),
            (HEADER + "1\t10\tx\tA\tC\t.\tPASS\t.\n" * 2, "duplicate")]:
        path.write_text(text)
        with pytest.raises(ValueError, match=message):
            read_vcf(path)


def test_gzip_and_blank_lines(tmp_path):
    path = tmp_path / "tiny.vcf.gz"
    with gzip.open(path, "wt") as stream:
        stream.write(HEADER + "\n1\t1\tx\tAC\tGT\t.\tPASS\t.\n")
    assert read_vcf(path)[0][0].change == "MNV"
    assert Variant("1", 1, "x", "ATG", "A").change == "deletion"


def test_overlap_inclusive_reference_span_and_multiple_genes():
    genes = [{"symbol": "A", "chrom": "1", "start": 10, "end": 20}, {"symbol": "B", "chrom": "1", "start": 20, "end": 30}]
    assert overlaps(Variant("1", 20, "x", "A", "G"), genes) == ["A", "B"]
    assert overlaps(Variant("1", 8, "x", "AAA", "A"), genes) == ["A"]
    assert overlaps(Variant("2", 20, "x", "A", "AT"), genes) == []


def test_missing_and_failed_vep_answers_are_explicit(config, store):
    variants = [Variant("1", 10, "x", "A", "C")]
    class Empty:
        def json(self, *args, **kwargs):
            return []
    rows = annotate(variants, [], Empty(), config, store)
    assert rows[0]["consequence"] == "unavailable" and rows[0]["vep_status"] == "not_returned"
    class Failed:
        def json(self, *args, **kwargs):
            raise ApiError("synthetic failure")
    assert annotate(variants, [], Failed(), config, store)[0]["vep_status"] == "request_failed"
    class Malformed:
        def json(self, *args, **kwargs):
            return {}
    assert annotate(variants, [], Malformed(), config, store)[0]["vep_status"] == "request_failed"
