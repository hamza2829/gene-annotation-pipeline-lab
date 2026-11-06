import warnings

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import ConvergenceWarning, PerfectSeparationError, PerfectSeparationWarning


def validate_inputs(matrix, phenotype):
    matrix, phenotype = np.asarray(matrix), np.asarray(phenotype)
    if matrix.ndim != 2 or phenotype.ndim != 1 or matrix.shape[0] != len(phenotype) or matrix.shape[1] == 0:
        raise ValueError("Genotypes must be samples by variants and match the phenotype length")
    if not np.isin(matrix, [0, 1, 2]).all() or not np.isin(phenotype, [0, 1]).all() or len(np.unique(phenotype)) != 2:
        raise ValueError("Use complete 0/1/2 synthetic genotypes and both binary phenotype classes")


def allele_test(genotype, phenotype):
    case, control = genotype[phenotype == 1], genotype[phenotype == 0]
    table = np.array([[case.sum(), 2 * len(case) - case.sum()], [control.sum(), 2 * len(control) - control.sum()]])
    if (table.sum(axis=0) == 0).any():
        return np.nan, np.nan, "monomorphic"
    value, pvalue, _, expected = chi2_contingency(table, correction=False)
    return float(value), float(pvalue), "low_expected_counts" if (expected < 5).any() else "ok"


def adjust(pvalues):
    values = np.asarray(pvalues, dtype=float)
    valid = np.isfinite(values)
    if ((values[valid] < 0) | (values[valid] > 1)).any():
        raise ValueError("P-values must be between zero and one, or missing")
    bonf, bh = np.full(len(values), np.nan), np.full(len(values), np.nan)
    indices = np.flatnonzero(valid)
    order = indices[np.argsort(values[indices], kind="stable")]
    bonf[indices] = np.minimum(values[indices] * len(values), 1)
    if len(order):
        ranked = values[order] * len(values) / np.arange(1, len(order) + 1)
        bh[order] = np.minimum(np.minimum.accumulate(ranked[::-1])[::-1], 1)
    return bonf, bh


def association(matrix, phenotype, markers):
    matrix, phenotype = np.asarray(matrix), np.asarray(phenotype)
    validate_inputs(matrix, phenotype)
    if len(markers) != matrix.shape[1]:
        raise ValueError("Marker metadata length does not match the genotype matrix")
    rows = []
    for index, marker in enumerate(markers):
        genotype = matrix[:, index].astype(float)
        statistic, chi_p, chi_status = allele_test(genotype, phenotype)
        beta, pvalue, status = np.nan, np.nan, "monomorphic"
        if np.ptp(genotype) != 0:
            try:
                with warnings.catch_warnings():
                    for category in [ConvergenceWarning, PerfectSeparationWarning, RuntimeWarning]:
                        warnings.simplefilter("error", category)
                    result = sm.Logit(phenotype, sm.add_constant(genotype)).fit(disp=False, maxiter=100)
                if not result.mle_retvals["converged"] or not np.isfinite(result.pvalues[1]):
                    raise ValueError("Logistic regression did not converge")
                beta, pvalue, status = float(result.params[1]), float(result.pvalues[1]), "ok"
            except (PerfectSeparationError, PerfectSeparationWarning, ConvergenceWarning, RuntimeWarning, ValueError, np.linalg.LinAlgError):
                status = "numerical_failure"
        rows.append({**marker, "beta": beta, "logistic_p": pvalue, "logistic_status": status, "chi2": statistic, "chi2_p": chi_p, "chi2_status": chi_status})
    frame = pd.DataFrame(rows)
    for name in ["logistic", "chi2"]:
        frame[f"{name}_bonferroni"], frame[f"{name}_bh"] = adjust(frame[f"{name}_p"])
    return frame
