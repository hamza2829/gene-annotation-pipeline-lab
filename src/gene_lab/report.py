import base64
from html import escape
import json


def table(frame, columns):
    return frame[columns].to_html(index=False, escape=True, border=0, na_rep="unavailable", float_format=lambda value: f"{value:.4g}")


def build_report(store):
    frames = {name: store.read(name) for name in ["genes", "transcripts", "proteins", "ncrna", "provider_status", "variants", "association", "expression"]}
    simulation = json.loads((store.directory / "association-summary.json").read_text())
    provenance = json.loads((store.directory / "run.json").read_text())
    cards = []
    for _, gene in frames["genes"].iterrows():
        symbol = gene["symbol"]
        proteins = frames["proteins"][frames["proteins"]["symbol"] == symbol]
        ncrna = frames["ncrna"][frames["ncrna"]["symbol"] == symbol]
        transcripts = frames["transcripts"][frames["transcripts"]["symbol"] == symbol]
        cards.append(f'<article class="gene-card"><h3>{escape(symbol)}</h3><p class="muted">{escape(gene["ensembl_id"])} / {escape(gene["biotype"])}</p>'
            f'<p>GRCh38 chromosome {escape(str(gene["chrom"]))}: {gene["start"]:,} to {gene["end"]:,}</p><p>{escape(str(gene["description"]))}</p>'
            f'<p><strong>{len(transcripts)}</strong> transcripts; <strong>{len(proteins)}</strong> reviewed protein results; <strong>{len(ncrna)}</strong> RNAcentral overlaps.</p>'
            '<details><summary>Protein and non-coding results</summary>' + table(proteins, ["accession", "name", "length"])
            + (f"<p>Coordinate matches, not gene-symbol equivalence. Showing up to 10 of {len(ncrna)} overlap features; complete results are in ncrna.csv.</p>" + table(ncrna.head(10), ["entry_id", "description", "mapping"]) if len(ncrna) else "<p>No RNAcentral overlap returned.</p>") + "</details></article>")

    def picture(name, alt):
        data = base64.b64encode((store.directory / name).read_bytes()).decode()
        return f'<img src="data:image/png;base64,{data}" alt="{escape(alt)}">'

    html = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Gene annotation pipeline lab</title><style>
*{{box-sizing:border-box}}body{{margin:0;background:#f1f5f9;color:#20374f;font:15px/1.6 Arial,sans-serif}}header{{background:#15364b;color:white;padding:36px max(28px,calc((100vw - 1200px)/2))}}h1{{font-size:32px;letter-spacing:-.8px;margin:10px 0}}header p{{color:#d0e2eb}}.eyebrow{{font-size:11px;letter-spacing:2px;color:#a1c7ce}}main{{max-width:1256px;margin:auto;padding:28px}}section{{background:white;padding:28px;border:1px solid #dce5ee;border-radius:10px;margin:22px 0}}h2{{font-size:23px;margin:0 0 14px}}h3{{margin:0;font-size:23px}}.muted{{color:#687e90;font-size:12px}}.notice{{background:#e4eff5;padding:16px;border-radius:6px}}.cards{{display:grid;grid-template-columns:1fr 1fr;gap:18px}}.gene-card{{padding:22px;border:1px solid #dce5ee;border-radius:8px;overflow:auto}}.gene-card p{{font-size:12px}}.scroll{{max-height:430px;overflow:auto}}table{{border-collapse:collapse;width:100%;font-size:12px}}th,td{{text-align:left;vertical-align:top;padding:10px;border-bottom:1px solid #e3eaf1}}th{{background:#eff4f8;position:sticky;top:0}}img{{max-width:100%;height:auto;display:block;margin:20px auto}}.plots{{display:grid;grid-template-columns:1.5fr 1fr;gap:15px;align-items:center}}summary{{cursor:pointer;font-size:12px;color:#226e83}}footer{{padding:20px 28px;color:#687e90;font-size:12px}}@media(max-width:800px){{.cards,.plots{{grid-template-columns:1fr}}main{{padding:15px}}section{{padding:20px}}}}
</style></head><body><header><div class="eyebrow">COMPUTATIONAL RESEARCH LAB / SYNTHETIC DATA</div><h1>Gene annotation pipeline</h1><p>Public gene annotations. Synthetic variants, genotypes and expression counts.</p></header><main>
<div class="notice">API mode: {escape(provenance['api_mode'])}. Reference assembly: GRCh38. All {simulation['samples']} simulated sample records are artificial. Consequence labels describe predicted sequence effects and are not clinical classifications.</div>
<section><h2>Genes of interest</h2><p>Public Ensembl genes and expanded transcripts, reviewed human UniProt proteins, and RNAcentral coordinate overlaps.</p><div class="cards">{''.join(cards)}</div></section>
<section><h2>Provider availability</h2><p>Empty results and failed requests are kept separate. Recorded mode uses captured public responses, not invented API annotations.</p><div class="scroll">{table(frames['provider_status'], ['symbol','provider','status','count','message'])}</div></section>
<section><h2>Synthetic VCF annotation</h2><p>{len(frames['variants'])} allele records. Local overlaps use the inclusive VCF REF span, with insertions assigned to their anchor base. Missing VEP responses are labeled unavailable.</p><div class="scroll">{table(frames['variants'], list(frames['variants'].columns))}</div></section>
<section><h2>Association simulation</h2><p>{simulation['samples']} artificial samples and {simulation['markers']} separate synthetic markers. Planted marker: <strong>{escape(simulation['planted_id'])}</strong>. Detected after logistic Bonferroni correction: <strong>{simulation['planted_detected']}</strong>.</p>
<p>A GWAS tests associations between many genetic variants and a trait. This small experiment fits one additive logistic model per marker and also tests allele counts. Association does not establish causation. This is a sparse simulation, not a genome-wide real cohort analysis. It omits ancestry structure, imputation, relatedness, covariates, linkage disequilibrium modeling and real quality control. The allelic chi square test treats allele counts as a simplified contingency table.</p>
<div class="plots">{picture('manhattan.png','Manhattan plot of the synthetic logistic association scan')}{picture('qq.png','QQ plot of synthetic logistic Wald p-values')}</div>
<p>Bonferroni and Benjamini Hochberg corrections are applied separately to each test family. All configured markers count toward the correction, including failed fits. Numerical failures stay missing rather than being assigned a successful p-value.</p><details><summary>All association results</summary><div class="scroll">{table(frames['association'],list(frames['association'].columns))}</div></details></section>
<section><h2>Toy single-cell example</h2><p>A synthetic genes-by-cells count matrix. Means and detection fractions are descriptive only. There is no normalization, clustering, batch correction or cell-type inference.</p>{picture('expression.png','Mean counts and nonzero detection in the toy synthetic expression matrix')}{table(frames['expression'],list(frames['expression'].columns))}</section>
<section><h2>Provenance and limits</h2><p>Seed {provenance['seed']}. {len(frames['genes'])} public gene regions. API recordings retain source URL, retrieval time and response hash. The association markers are a separate experiment from the 200 VCF variants. No real cohort, clinical phenotype or person data is used.</p></section>
</main><footer>Gene annotation pipeline lab. Synthetic teaching example with public reference annotations.</footer></body></html>'''
    path = store.directory / "report.html"
    path.write_text(html, encoding="utf-8")
    return path
