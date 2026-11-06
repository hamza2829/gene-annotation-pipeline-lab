import json
from pathlib import Path

from playwright.sync_api import sync_playwright

from gene_lab.pipeline import verify_outputs

root = Path("artifacts")
output = root / "container" / "run"
proof = verify_outputs(output)
with sync_playwright() as playwright:
    browser = playwright.chromium.launch()
    page = browser.new_page(viewport={"width": 1440, "height": 1100})
    errors = []
    page.on("pageerror", lambda error: errors.append(str(error)))
    page.goto((output / "report.html").resolve().as_uri())
    assert page.locator(".gene-card").count() == 10
    assert page.locator("img").count() == 3
    assert page.evaluate("Array.from(document.images).every(image => image.complete && image.naturalWidth > 0)")
    page.screenshot(path=str(root / "report-screen.png"), full_page=True)
    assert not errors, errors
    browser.close()
proof.update({"docker_network": "none", "browser_errors": errors, "gene_cards": 10, "embedded_plots": 3})
(root / "pipeline-proof.json").write_text(json.dumps(proof, indent=2) + "\n")
print(json.dumps(proof, indent=2))
