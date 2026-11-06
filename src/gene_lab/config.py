from pathlib import Path
import re

import yaml


def load_config(path):
    path = Path(path).resolve()
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise ValueError(f"{path}: expected a YAML mapping")
    needed = {"assembly", "species", "ensembl_url", "uniprot_url", "rnacentral_url", "recordings", "vcf", "expression", "output", "cache", "timeout", "attempts", "backoff", "max_backoff", "min_interval", "cache_ttl", "vep_batch", "seed", "samples", "markers", "planted_index", "effect_beta", "alpha"}
    if set(config) != needed:
        raise ValueError(f"{path}: missing or unknown config keys: {sorted(needed ^ set(config))}")
    if config["assembly"] != "GRCh38" or config["species"] != "homo_sapiens":
        raise ValueError("This lab supports homo_sapiens GRCh38 only")
    for key in ["timeout", "attempts", "cache_ttl", "samples", "markers", "vep_batch"]:
        if not isinstance(config[key], (int, float)) or config[key] <= 0:
            raise ValueError(f"{key} must be positive")
    if config["attempts"] > 10 or not 0 < config["alpha"] < 1 or not 0 <= config["planted_index"] < config["markers"]:
        raise ValueError("Invalid attempts, alpha or planted index")
    for key in ["backoff", "max_backoff", "min_interval"]:
        if config[key] < 0:
            raise ValueError(f"{key} must be nonnegative")
    for key in ["recordings", "vcf", "expression", "output", "cache"]:
        config[key] = (path.parent / config[key]).resolve()
    return config


def read_genes(path):
    genes = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        symbol = line.strip().upper()
        if not symbol or symbol.startswith("#"):
            continue
        if not re.fullmatch(r"[A-Z0-9][A-Z0-9_.-]{0,39}", symbol):
            raise ValueError(f"{path}: invalid gene symbol {symbol!r}")
        if symbol not in genes:
            genes.append(symbol)
    if not genes:
        raise ValueError(f"{path}: gene list is empty")
    return genes
