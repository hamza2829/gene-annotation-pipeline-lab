import json
from pathlib import Path

import numpy as np
import pandas as pd

from .variants import Variant


def make_vcf(genes, sequences, path):
    rows = []
    for gene in genes:
        sequence = sequences[gene["symbol"]].upper()
        if len(sequence) != int(gene["end"]) - int(gene["start"]) + 1:
            raise ValueError("Reference sequence length does not match the gene region")
        positions = np.linspace(5, len(sequence) - 6, 20, dtype=int)
        for index, offset in enumerate(positions):
            ref = sequence[offset:offset + (3 if index >= 16 else 1)]
            alt = next(base for base in "ACGT" if base != ref[0]) if index < 12 else ref + "AC" if index < 16 else ref[0]
            rows.append(Variant(str(gene["chrom"]), int(gene["start"]) + int(offset), f"SYN{len(rows) + 1:04}", ref, alt))
    lines = ["##fileformat=VCFv4.2", "##reference=GRCh38", "##source=Synthetic_variants_using_public_reference_sequence", '##INFO=<ID=SYNTHETIC,Number=1,Type=Integer,Description="Artificial variant marker">', "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    lines += [f"{row.chrom}\t{row.pos}\t{row.identifier}\t{row.ref}\t{row.alt}\t.\tPASS\tSYNTHETIC=1" for row in rows]
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return rows


def make_expression(symbols, path, seed=2829):
    rng = np.random.default_rng(seed)
    counts = rng.poisson(np.linspace(1, 8, len(symbols))[:, None], size=(len(symbols), 24))
    counts[rng.random(counts.shape) < 0.2] = 0
    pd.DataFrame(counts, index=pd.Index(symbols, name="gene"), columns=[f"CELL{index + 1:03}" for index in range(24)]).to_csv(path)


def genotypes(genes, config, directory):
    samples, count = int(config["samples"]), int(config["markers"])
    if samples < 20 or count < 2 or not 0 <= config["planted_index"] < count:
        raise ValueError("Simulation requires at least 20 samples, two markers and a valid planted index")
    if not genes:
        raise ValueError("Gene coordinates are required for marker positions")
    rng = np.random.default_rng(config["seed"])
    frequencies = rng.uniform(0.1, 0.45, count)
    matrix = rng.binomial(2, frequencies, size=(samples, count))
    linear = -1.3 + config["effect_beta"] * matrix[:, config["planted_index"]]
    phenotype = rng.binomial(1, 1 / (1 + np.exp(-linear)))
    markers = []
    for index in range(count):
        gene = genes[index % len(genes)]
        fraction = (index // len(genes) + 1) / (count // len(genes) + 2)
        markers.append({"variant_id": f"SIM{index + 1:04}", "chrom": str(gene["chrom"]), "pos": int(gene["start"] + fraction * (gene["end"] - gene["start"])), "gene": gene["symbol"]})
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(directory / "synthetic-genotypes.npz", genotypes=matrix, phenotype=phenotype)
    pd.DataFrame(markers).to_csv(directory / "simulation-markers.csv", index=False)
    pd.DataFrame({"sample_id": [f"SYNTHETIC{index + 1:04}" for index in range(samples)], "synthetic_case": phenotype}).to_csv(directory / "synthetic-samples.csv", index=False)
    (directory / "simulation-design.json").write_text(json.dumps({"synthetic_only": True, "samples": samples, "markers": count,
        "planted_id": markers[config["planted_index"]]["variant_id"], "effect_beta": config["effect_beta"], "seed": config["seed"]}, indent=2) + "\n")
    return matrix, phenotype, markers
