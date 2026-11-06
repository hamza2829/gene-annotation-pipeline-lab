from urllib.parse import quote

from .http import ApiError


def fetch_genes(symbols, client, config, store):
    genes, transcripts, proteins, ncrna, statuses, raw = [], [], [], [], [], {}
    for symbol in symbols:
        raw[symbol] = {}
        gene = None
        try:
            lookup = f"{config['ensembl_url']}/lookup/symbol/{config['species']}/{quote(symbol)}"
            try:
                gene = client.json("GET", lookup, {"expand": 1})
            except ApiError:
                gene = client.json("GET", lookup)
                if not isinstance(gene, dict) or "id" not in gene:
                    raise ApiError("Basic Ensembl lookup omitted a gene id")
                transcript_rows = client.json("GET", f"{config['ensembl_url']}/overlap/id/{gene['id']}", {"feature": "transcript"})
                if not isinstance(transcript_rows, list):
                    raise ApiError("Ensembl transcript overlap is not an array")
                gene["Transcript"] = transcript_rows
            if not isinstance(gene, dict) or not all(key in gene for key in ["id", "seq_region_name", "start", "end", "assembly_name"]):
                raise ApiError("Ensembl lookup omitted required gene coordinates")
            if gene["assembly_name"] != config["assembly"]:
                raise ApiError(f"Assembly mismatch: {gene['assembly_name']} versus {config['assembly']}")
            genes.append({"symbol": symbol, "ensembl_id": gene["id"], "chrom": str(gene["seq_region_name"]), "start": gene["start"],
                "end": gene["end"], "strand": gene.get("strand", 1), "biotype": gene.get("biotype", "unknown"),
                "assembly": gene["assembly_name"], "description": gene.get("description", "")})
            for transcript in gene.get("Transcript", []):
                transcripts.append({"symbol": symbol, "transcript_id": transcript["id"], "start": transcript["start"], "end": transcript["end"], "biotype": transcript.get("biotype", "unknown")})
            raw[symbol]["ensembl"] = gene
            statuses.append({"symbol": symbol, "provider": "Ensembl", "status": "ok", "count": len(gene.get("Transcript", [])), "message": "transcripts returned"})
        except ApiError as error:
            gene = None
            statuses.append({"symbol": symbol, "provider": "Ensembl", "status": "error", "count": 0, "message": str(error)})
        try:
            data = client.json("GET", config["uniprot_url"] + "/uniprotkb/search", {"query": f"gene_exact:{symbol} AND organism_id:9606 AND reviewed:true", "format": "json", "size": 3})
            if not isinstance(data, dict) or not isinstance(data.get("results"), list):
                raise ApiError("UniProt response has no results array")
            raw[symbol]["uniprot"] = data
            for protein in data["results"]:
                description = protein.get("proteinDescription", {})
                name = description.get("recommendedName", {}).get("fullName", {}).get("value", "Unnamed protein")
                functions = [text["value"] for comment in protein.get("comments", []) if comment.get("commentType") == "FUNCTION" for text in comment.get("texts", [])]
                proteins.append({"symbol": symbol, "accession": protein["primaryAccession"], "name": name,
                    "length": protein.get("sequence", {}).get("length"), "function": " ".join(functions)})
            statuses.append({"symbol": symbol, "provider": "UniProt", "status": "ok" if data["results"] else "empty", "count": len(data["results"]), "message": "reviewed human exact gene-symbol search, up to three proteins"})
        except ApiError as error:
            statuses.append({"symbol": symbol, "provider": "UniProt", "status": "error", "count": 0, "message": str(error)})
        if gene is not None and "assembly_name" in gene and gene.get("assembly_name") == config["assembly"]:
            try:
                region = f"{gene['seq_region_name']}:{gene['start']}-{gene['end']}"
                data = client.json("GET", f"{config['rnacentral_url']}/overlap/region/homo_sapiens/{region}", {"format": "json"})
                entries = data if isinstance(data, list) else data.get("results", data.get("features", []))
                if not isinstance(entries, list) or any(not isinstance(entry, dict) for entry in entries):
                    raise ApiError("RNAcentral response has no overlap array")
                raw[symbol]["rnacentral"] = data
                for index, entry in enumerate(entries):
                    mapping = f"{entry.get('seq_region_name', gene['seq_region_name'])}:{entry.get('start', gene['start'])}-{entry.get('end', gene['end'])}"
                    ncrna.append({"symbol": symbol, "entry_id": str(entry.get("ID", entry.get("id", entry.get("rnacentral_id", f"overlap-{index + 1}")))),
                        "description": str(entry.get("description", entry.get("biotype", "RNAcentral coordinate overlap"))), "mapping": mapping})
                statuses.append({"symbol": symbol, "provider": "RNAcentral", "status": "ok" if entries else "empty", "count": len(entries), "message": "coordinate overlap, not gene-symbol equivalence"})
            except ApiError as error:
                statuses.append({"symbol": symbol, "provider": "RNAcentral", "status": "error", "count": 0, "message": str(error)})
        else:
            statuses.append({"symbol": symbol, "provider": "RNAcentral", "status": "unavailable", "count": 0, "message": "No usable GRCh38 gene coordinates"})
    for name, rows in [("genes", genes), ("transcripts", transcripts), ("proteins", proteins), ("ncrna", ncrna), ("provider_status", statuses)]:
        store.write(name, rows)
    store.metadata("raw_api_bundles", raw)
    return genes
