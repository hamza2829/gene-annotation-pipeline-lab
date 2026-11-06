"""Capture genuine public responses on the runner before freezing the demo fixtures."""
import json
from pathlib import Path

from gene_lab.config import load_config, read_genes
from gene_lab.http import CaptureSession, HttpClient
from gene_lab.pipeline import fetch_stage
from gene_lab.storage import Store
from gene_lab.synthetic import make_expression, make_vcf
from gene_lab.variants import annotate

config = load_config("config.yaml")
symbols = read_genes("data/genes.txt")
capture = CaptureSession(config["recordings"])
client = HttpClient(config, capture)
store = Store("artifacts/recording-run")
try:
    genes = fetch_stage(symbols, client, config, store)
    sequences = {}
    for gene in genes:
        region = f"{gene['chrom']}:{gene['start']}..{gene['end']}:1"
        result = client.json("GET", config["ensembl_url"] + f"/sequence/region/{config['species']}/{region}")
        sequences[gene["symbol"]] = result["seq"]
    variants = make_vcf(genes, sequences, config["vcf"])
    make_expression(symbols, config["expression"], config["seed"])
    annotations = annotate(variants, genes, client, config, store)
    assert len(annotations) == 200 and sum(row["vep_status"] == "ok" for row in annotations) > 0
    assert len(store.read("proteins")) > 0 and len(store.read("ncrna")) > 0, "Public protein and RNAcentral results are required for the demo"
    manifest = []
    for path in sorted(Path(config["recordings"]).glob("*.json")):
        record = json.loads(path.read_text())
        manifest.append({"file": path.name, "method": record["method"], "url": record["url"], "retrieved_at": record["retrieved_at"], "wire_sha256": record["wire_sha256"]})
    (Path(config["recordings"]) / "index.json").write_text(json.dumps({"source": "Actual public HTTP responses captured in GitHub Actions", "records": manifest}, indent=2) + "\n")
finally:
    client.close()
    store.close()
