import numpy as np
import pytest

from gene_lab.association import adjust, allele_test, association, validate_inputs
from gene_lab.plots import make_plots, qq_values
from gene_lab.synthetic import genotypes

GENES = [{"symbol": "SYN_A", "chrom": "1", "start": 100, "end": 10000}, {"symbol": "SYN_B", "chrom": "2", "start": 300, "end": 20000}]


def test_planted_effect_detected_by_both_tests_and_corrections(config, tmp_path):
    matrix, phenotype, markers = genotypes(GENES, config, tmp_path)
    assert matrix.shape == (2000, 300) and phenotype.shape == (2000,)
    frame = association(matrix, phenotype, markers)
    planted = frame.iloc[42]
    assert planted["variant_id"] == "SIM0043" and planted["beta"] > 0
    assert planted["logistic_bonferroni"] < 0.05 and planted["logistic_bh"] < 0.05
    assert planted["chi2_bonferroni"] < 0.05 and planted["chi2_bh"] < 0.05
    assert frame.sort_values("logistic_p").iloc[0]["variant_id"] == "SIM0043"
    assert frame["logistic_status"].eq("ok").all()


def test_known_allele_count_chi_square():
    genotype = np.array([2, 2, 1, 1, 1, 1, 0, 0])
    phenotype = np.array([1, 1, 1, 1, 0, 0, 0, 0])
    statistic, pvalue, status = allele_test(genotype, phenotype)
    assert statistic == pytest.approx(4)
    assert pvalue == pytest.approx(0.0455002639)
    assert status == "low_expected_counts"


def test_bonferroni_and_bh_hand_example():
    bonf, bh = adjust([0.01, 0.04, 0.03])
    assert bonf == pytest.approx([0.03, 0.12, 0.09])
    assert bh == pytest.approx([0.03, 0.04, 0.04])


def test_missing_tests_still_count_in_multiplicity():
    bonf, bh = adjust([0.01, 0.04, 0.03, np.nan])
    assert bonf[:3] == pytest.approx([0.04, 0.16, 0.12])
    assert bh[:3] == pytest.approx([0.04, 0.04 * 4 / 3, 0.04 * 4 / 3])
    assert np.isnan(bh[-1]) and np.isnan(bonf[-1])
    assert len(adjust([])[0]) == 0
    with pytest.raises(ValueError, match="P-values"):
        adjust([1.1])


@pytest.mark.parametrize("matrix,phenotype", [(np.zeros((4, 2)), [1, 1]), (np.zeros(4), [0, 1, 0, 1]),
    (np.zeros((4, 0)), [0, 1, 0, 1]), (np.array([[3], [1], [0], [0]]), [0, 1, 0, 1]),
    (np.zeros((4, 1)), [0, 0, 0, 0]), (np.array([[np.nan], [1], [0], [0]]), [0, 1, 0, 1])])
def test_invalid_genotypes_and_phenotypes(matrix, phenotype):
    with pytest.raises(ValueError):
        validate_inputs(np.asarray(matrix), np.asarray(phenotype))


def test_monomorphic_and_separation_are_not_successful_fits():
    phenotype = np.array([0] * 20 + [1] * 20)
    matrix = np.column_stack([np.zeros(40), phenotype * 2])
    markers = [{"variant_id": "mono", "chrom": "1", "pos": 1, "gene": "A"}, {"variant_id": "sep", "chrom": "1", "pos": 2, "gene": "A"}]
    frame = association(matrix, phenotype, markers)
    assert frame.iloc[0]["logistic_status"] == "monomorphic"
    assert frame.iloc[0]["chi2_status"] == "monomorphic"
    assert frame.iloc[1]["logistic_status"] == "numerical_failure"
    assert frame["logistic_p"].isna().all()
    with pytest.raises(ValueError, match="metadata"):
        association(matrix, phenotype, markers[:1])


def test_unconverged_model_is_reported(monkeypatch):
    import gene_lab.association as module
    class BadFit:
        mle_retvals = {"converged": False}
        pvalues = [1, 0.01]
    class FakeModel:
        def fit(self, **kwargs):
            return BadFit()
    monkeypatch.setattr(module.sm, "Logit", lambda *args: FakeModel())
    frame = association(np.array([[0], [1], [1], [2]]), np.array([0, 1, 0, 1]), [{"variant_id": "x", "chrom": "1", "pos": 1}])
    assert frame.iloc[0]["logistic_status"] == "numerical_failure"


def test_qq_quantiles_zero_p_and_missing():
    expected, observed = qq_values([0.5, 0.01, np.nan, 0])
    assert expected == pytest.approx(-np.log10([1 / 6, 0.5, 5 / 6]))
    assert np.isfinite(observed).all() and observed[1] == 2
    assert len(qq_values([])[0]) == 0
    with pytest.raises(ValueError, match="Invalid p"):
        qq_values([1.5])


def test_simulation_reproducible_and_invalid_design(config, tmp_path):
    one = genotypes(GENES, config, tmp_path / "one")
    two = genotypes(GENES, config, tmp_path / "two")
    assert np.array_equal(one[0], two[0]) and np.array_equal(one[1], two[1])
    with pytest.raises(ValueError, match="coordinates"):
        genotypes([], config, tmp_path)
    with pytest.raises(ValueError, match="20 samples"):
        genotypes(GENES, {**config, "samples": 10}, tmp_path)


def test_plot_artifacts_with_missing_fits(config, tmp_path):
    import pandas as pd
    frame = pd.DataFrame({"variant_id": ["A", "B"], "chrom": ["X", "2"], "pos": [10, 20], "logistic_p": [np.nan, 0.02]})
    make_plots(frame, tmp_path, 0.05, "A")
    assert (tmp_path / "manhattan.png").stat().st_size > 1000
    assert (tmp_path / "qq.png").stat().st_size > 1000
