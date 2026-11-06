import json
import os
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET

root = Path("artifacts")
cases = list(ET.parse(root / "tests.xml").getroot().iter("testcase"))
skipped = [case for case in cases if case.find("skipped") is not None]
passed = [case for case in cases if all(case.find(tag) is None for tag in ["skipped", "failure", "error"])]
assert len(passed) >= 40 and len(skipped) == 1
assert all("test_live" in case.attrib["classname"] for case in skipped)
assert len(passed) + len(skipped) == len(cases)
coverage = ET.parse(root / "coverage.xml").getroot()
covered, valid = int(coverage.attrib["lines-covered"]), int(coverage.attrib["lines-valid"])
assert covered / valid >= 0.85
pipeline = json.loads((root / "pipeline-proof.json").read_text())
metrics = {"passed_tests": len(passed), "live_skipped": len(skipped), "failures": 0, "line_coverage_percent": round(100 * covered / valid, 2),
    "covered_lines": covered, "measured_lines": valid, "pipeline": pipeline, "commit": os.environ["GITHUB_SHA"],
    "run_url": f'https://github.com/{os.environ["GITHUB_REPOSITORY"]}/actions/runs/{os.environ["GITHUB_RUN_ID"]}'}
(root / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
print(json.dumps(metrics, indent=2))
with open(os.environ["GITHUB_STEP_SUMMARY"], "a", encoding="utf-8") as summary:
    summary.write(f'| Metric | Measured in CI |\n| --- | ---: |\n| Passed pytest cases | {len(passed)} |\n| Live tests skipped | {len(skipped)} |\n| Line coverage | {metrics["line_coverage_percent"]}% ({covered}/{valid}) |\n')
for name in subprocess.run(["git", "ls-files", "-z"], check=True, capture_output=True, text=True).stdout.split("\0"):
    if name and Path(name).suffix not in {".png", ".npz"}:
        assert chr(0x2014) not in Path(name).read_text(encoding="utf-8"), name
