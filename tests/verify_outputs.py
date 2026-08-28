#!/usr/bin/env python3
"""Small dependency-free semantic checks for the xamdst integration fixture."""
import gzip
import json
import pathlib
import sys

out = pathlib.Path(sys.argv[1])
expected = {
    "coverage.report",
    "coverage.report.json",
    "cumu.plot",
    "insert.plot",
    "chromosome.report",
    "region.tsv.gz",
    "depth.tsv.gz",
    "uncover.bed",
}
actual = {p.name for p in out.iterdir() if p.is_file()}
missing = expected - actual
if missing:
    raise SystemExit(f"missing outputs: {sorted(missing)}")
legacy = {"depth_distribution.plot", "insertsize.plot", "chromosomes.report"}
if actual & legacy:
    raise SystemExit(f"legacy outputs still present: {sorted(actual & legacy)}")

with (out / "coverage.report.json").open(encoding="utf-8") as handle:
    report = json.load(handle)
assert report["schema_version"] == "3.1"
assert report["version"] == "3.1.0"
assert report["depth_output"] is True
assert report["target"]["region_length"] == 10
assert report["flank"]["region_length"] == 4

with gzip.open(out / "depth.tsv.gz", "rt", encoding="utf-8") as handle:
    lines = list(handle)
assert lines[0].startswith("#Chr\tPos\tRaw depth\tRmdup depth\tCoverage with deletions")
rows = [line.rstrip("\n").split("\t") for line in lines if not line.startswith("#")]
assert len(rows) == 10, len(rows)
assert rows[0][0] == "chr1" and rows[0][1] == "3"
# SAM positions are 1-based: the first target base is covered by r1 only;
# the next two bases overlap r1 and duplicate r2, while rmdup keeps r1.
assert rows[0][2:] == ["1", "1", "1"], rows[0]
assert rows[1][2:] == ["2", "1", "2"], rows[1]
# The N operation advances reference coordinates and must not add depth.
assert rows[8][1] == "11"
assert rows[9][1] == "12"

with gzip.open(out / "region.tsv.gz", "rt", encoding="utf-8") as handle:
    region_rows = [line for line in handle if not line.startswith("#")]
assert len(region_rows) == 1
