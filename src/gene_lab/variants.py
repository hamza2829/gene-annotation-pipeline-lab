from dataclasses import dataclass
import gzip
from pathlib import Path
import re

from .http import ApiError


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int
    identifier: str
    ref: str
    alt: str

    @property
    def key(self):
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"

    @property
    def change(self):
        if len(self.ref) == len(self.alt) == 1:
            return "SNV"
        if len(self.alt) > len(self.ref):
            return "insertion"
        return "deletion" if len(self.alt) < len(self.ref) else "MNV"

    def vep_input(self):
        return f"{self.chrom} {self.pos} {self.identifier} {self.ref} {self.alt} . . ."


def read_vcf(path, assembly="GRCh38"):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    variants, issues, seen, header = [], [], set(), False
    with opener(path, "rt", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith("##reference=") and line.split("=", 1)[1] != assembly:
                raise ValueError(f"{path}:{line_number}: VCF assembly does not match {assembly}")
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                if line.split("\t")[:8] != ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]:
                    raise ValueError(f"{path}:{line_number}: invalid VCF column header")
                header = True
                continue
            if not header:
                raise ValueError(f"{path}:{line_number}: missing #CHROM header")
            fields = line.split("\t")
            try:
                if len(fields) < 8 or int(fields[1]) < 1:
                    raise ValueError
                chrom, pos, identifier, ref, alts = fields[:5]
                chrom, ref = chrom.removeprefix("chr"), ref.upper()
                if not chrom or not re.fullmatch(r"[ACGTN]+", ref):
                    raise ValueError
                for alt in alts.split(","):
                    alt = alt.upper()
                    if not re.fullmatch(r"[ACGTN]+", alt):
                        issues.append({"line": line_number, "message": f"Unsupported or missing ALT {alt!r}; allele skipped"})
                        continue
                    if alt == ref:
                        issues.append({"line": line_number, "message": "REF equals ALT; allele skipped"})
                        continue
                    variant = Variant(chrom, int(pos), identifier, ref, alt)
                    if variant.key in seen:
                        raise ValueError("duplicate allele")
                    seen.add(variant.key)
                    variants.append(variant)
            except (ValueError, IndexError) as error:
                raise ValueError(f"{path}:{line_number}: invalid VCF row ({error})") from error
    if not header:
        raise ValueError(f"{path}: missing #CHROM header")
    return variants, issues


def overlaps(variant, genes):
    end = variant.pos + len(variant.ref) - 1
    return sorted(gene["symbol"] for gene in genes if str(gene["chrom"]) == variant.chrom and int(gene["start"]) <= end and int(gene["end"]) >= variant.pos)


def annotate(variants, genes, client, config, store):
    answers, failures = {}, set()
    size = int(config["vep_batch"])
    for offset in range(0, len(variants), size):
        batch = variants[offset:offset + size]
        try:
            result = client.json("POST", config["ensembl_url"] + f"/vep/{config['species']}/region", body={"variants": [row.vep_input() for row in batch]})
            if not isinstance(result, list):
                raise ApiError("VEP response is not an array")
            for item in result:
                key = " ".join(str(item.get("input", "")).split())
                answers[key] = item.get("most_severe_consequence")
        except ApiError:
            failures.update(row.key for row in batch)
    rows = []
    for variant in variants:
        consequence = answers.get(variant.vep_input())
        rows.append({"variant_id": variant.identifier, "chrom": variant.chrom, "pos": variant.pos, "ref": variant.ref, "alt": variant.alt,
            "change": variant.change, "overlapping_genes": ";".join(overlaps(variant, genes)), "consequence": consequence or "unavailable",
            "vep_status": "ok" if consequence else "request_failed" if variant.key in failures else "not_returned"})
    store.write("variants", rows)
    return rows
