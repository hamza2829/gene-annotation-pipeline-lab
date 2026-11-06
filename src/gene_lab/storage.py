import json
from pathlib import Path
import sqlite3

import pandas as pd

SCHEMAS = {
    "genes": ["symbol", "ensembl_id", "chrom", "start", "end", "strand", "biotype", "assembly", "description"],
    "transcripts": ["symbol", "transcript_id", "start", "end", "biotype"],
    "proteins": ["symbol", "accession", "name", "length", "function"],
    "ncrna": ["symbol", "entry_id", "description", "mapping"],
    "provider_status": ["symbol", "provider", "status", "count", "message"],
    "variants": ["variant_id", "chrom", "pos", "ref", "alt", "change", "overlapping_genes", "consequence", "vep_status"],
    "vcf_issues": ["line", "message"],
    "expression": ["gene", "cells", "mean", "median", "detected_fraction", "status"],
}


class Store:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.connection = sqlite3.connect(self.directory / "pipeline.sqlite")

    def write(self, name, rows):
        if name not in SCHEMAS and name != "association":
            raise ValueError(f"Unknown output table: {name}")
        frame = pd.DataFrame(rows, columns=SCHEMAS.get(name))
        frame.to_sql(name, self.connection, if_exists="replace", index=False)
        frame.to_csv(self.directory / f"{name}.csv", index=False)
        return frame

    def read(self, name):
        if name not in SCHEMAS and name != "association":
            raise ValueError(f"Unknown output table: {name}")
        try:
            return pd.read_sql_query(f'SELECT * FROM "{name}"', self.connection)
        except pd.errors.DatabaseError as error:
            raise ValueError(f"Missing pipeline stage: {name}; run the corresponding command first") from error

    def metadata(self, name, data):
        (self.directory / f"{name}.json").write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    def close(self):
        self.connection.close()
