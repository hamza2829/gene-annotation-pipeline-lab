import json
import math

from .association import association
from .expression import summarize_expression
from .plots import make_plots
from .providers import fetch_genes
from .report import build_report
from .synthetic import genotypes
from .variants import annotate, read_vcf


def fetch_stage(symbols, client, config, store):
    genes = fetch_genes(symbols, client, config, store)
    store.metadata("run", {"api_mode": "recorded public responses" if client.offline else "live public API with cache", "seed": config["seed"], "assembly": config["assembly"]})
    if len(genes) != len(symbols):
        raise ValueError("Some Ensembl gene lookups failed; see provider_status.csv for every finding")
    return genes


def annotate_stage(client, config, store):
    variants, issues = read_vcf(config["vcf"], config["assembly"])
    genes = store.read("genes").to_dict("records")
    rows = annotate(variants, genes, client, config, store)
    store.write("vcf_issues", issues)
    return rows


def associate_stage(config, store):
    genes = store.read("genes").to_dict("records")
    matrix, phenotype, markers = genotypes(genes, config, store.directory)
    frame = association(matrix, phenotype, markers)
    store.write("association", frame.to_dict("records"))
    planted_id = markers[config["planted_index"]]["variant_id"]
    planted = frame[frame["variant_id"] == planted_id].iloc[0]
    detected = bool(planted["logistic_bonferroni"] < config["alpha"] and planted["chi2_bonferroni"] < config["alpha"])
    store.metadata("association-summary", {"samples": int(matrix.shape[0]), "markers": int(matrix.shape[1]), "cases": int(phenotype.sum()),
        "controls": int(len(phenotype) - phenotype.sum()), "planted_id": planted_id, "planted_detected": detected,
        "logistic_p": float(planted["logistic_p"]) if math.isfinite(planted["logistic_p"]) else None,
        "logistic_bonferroni": float(planted["logistic_bonferroni"]) if math.isfinite(planted["logistic_bonferroni"]) else None,
        "chi2_p": float(planted["chi2_p"]) if math.isfinite(planted["chi2_p"]) else None, "failed_logistic_fits": int((frame["logistic_status"] != "ok").sum())})
    make_plots(frame, store.directory, config["alpha"], planted_id)
    if not detected:
        raise ValueError("The planted marker did not pass both Bonferroni tests for this simulation configuration")
    return frame


def report_stage(symbols, config, store):
    summarize_expression(config["expression"], symbols, store)
    return build_report(store)


def verify_outputs(directory):
    from pathlib import Path
    import pandas as pd
    path = Path(directory)
    summary = json.loads((path / "association-summary.json").read_text())
    counts = {name: len(pd.read_csv(path / f"{name}.csv")) for name in ["genes", "variants", "association", "transcripts", "proteins", "ncrna"]}
    if counts["genes"] != 10 or counts["variants"] != 200 or counts["association"] != 300 or not summary["planted_detected"]:
        raise ValueError("The full demo output does not match its declared dimensions or planted effect")
    return {"counts": counts, "association": summary}
