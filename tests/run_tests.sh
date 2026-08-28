#!/bin/sh
set -eu

if [ "$#" -ne 1 ]; then
    echo "usage: $0 path/to/xamdst" >&2
    exit 2
fi
binary=$1
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
fixture="$root/tests/fixtures"
tmp=$(mktemp -d "${TMPDIR:-/tmp}/xamdst-test.XXXXXX")
trap 'rm -rf "$tmp"' EXIT HUP INT TERM

hash_file() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    else
        shasum -a 256 "$1" | awk '{print $1}'
    fi
}

"$binary" -p "$fixture/target.bed" -o "$tmp/out" --flank 2 "$fixture/mini.sam"
python3 "$fixture/../verify_outputs.py" "$tmp/out"
gzip -t "$tmp/out/depth.tsv.gz" "$tmp/out/region.tsv.gz"
before=$(hash_file "$tmp/out/coverage.report")
if "$binary" -p "$fixture/target.bed" -o "$tmp/out" --flank 2 "$fixture/unsorted.sam" >/dev/null 2>&1; then
    echo "unsorted replacement unexpectedly succeeded" >&2
    exit 1
fi
after=$(hash_file "$tmp/out/coverage.report")
[ "$before" = "$after" ]

"$binary" -p "$fixture/target.bed" -o "$tmp/bamout" --flank 2 \
    --threads 2 --compute-threads 4 --bamout "$tmp/bamout/target.bam" "$fixture/mini.sam"
python3 - "$tmp/bamout/target.bam" <<'PY'
import gzip
import pathlib
import sys
path = pathlib.Path(sys.argv[1])
assert path.read_bytes()[:2] == b"\x1f\x8b"
with gzip.open(path, "rb") as handle:
    assert handle.read(4) == b"BAM\1"
PY
"$binary" -p "$fixture/target.bed" -o "$tmp/bamout-serial" --flank 2 \
    --bamout "$tmp/bamout-serial/target.bam" "$fixture/mini.sam"
cmp "$tmp/bamout/target.bam" "$tmp/bamout-serial/target.bam"

"$binary" -p "$fixture/target.bed" -o "$tmp/merged" --flank 2 \
    "$fixture/mini.sam" "$fixture/mini.sam"
python3 - "$tmp/out/depth.tsv.gz" "$tmp/merged/depth.tsv.gz" <<'PY'
import gzip
import sys

def rows(path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line for line in handle if not line.startswith("#")]

one = rows(sys.argv[1])
two = rows(sys.argv[2])
assert len(one) == len(two)
for left, right in zip(one, two):
    l = left.rstrip().split("\t")
    r = right.rstrip().split("\t")
    assert l[:2] == r[:2]
    assert [int(v) for v in r[2:]] == [2 * int(v) for v in l[2:]], (l, r)
PY
if "$binary" -p "$fixture/target.bed" -o "$tmp/header-mismatch" \
    "$fixture/mini.sam" "$fixture/mismatch.sam" >/dev/null 2>&1; then
    echo "incompatible headers unexpectedly succeeded" >&2
    exit 1
fi

"$binary" -p "$fixture/target.bed" -o "$tmp/summary" --flank 2 --summary-only \
    --compute-threads 4 "$fixture/mini.sam"
printf 'keep-old-depth\n' > "$tmp/summary/depth.tsv.gz"
if "$binary" -p "$fixture/target.bed" -o "$tmp/summary" --flank 2 --summary-only \
    "$fixture/unsorted.sam" >/dev/null 2>&1; then
    echo "summary failure unexpectedly succeeded" >&2
    exit 1
fi
grep -q '^keep-old-depth$' "$tmp/summary/depth.tsv.gz"
"$binary" -p "$fixture/target.bed" -o "$tmp/summary" --flank 2 --summary-only \
    --compute-threads 2 "$fixture/mini.sam"
python3 - "$tmp/summary/coverage.report.json" "$tmp/summary" <<'PY'
import json
import pathlib
import sys
report = json.load(open(sys.argv[1], encoding="utf-8"))
assert report["schema_version"] == "3.1"
assert report["depth_output"] is False
assert not (pathlib.Path(sys.argv[2]) / "depth.tsv.gz").exists()
PY

"$binary" -p "$fixture/target.bed" -o "$tmp/fragment" --flank 0 --fragment-mode \
    "$fixture/fragment.sam"
python3 - "$tmp/fragment/depth.tsv.gz" <<'PY'
import gzip
import sys
rows = [line.rstrip().split("\t") for line in gzip.open(sys.argv[1], "rt")
        if not line.startswith("#")]
assert [int(row[2]) for row in rows[:6]] == [1, 1, 1, 1, 1, 1], rows
PY
"$binary" -p "$fixture/target.bed" -o "$tmp/fragment-parallel" --flank 0 \
    --fragment-mode --compute-threads 4 "$fixture/fragment.sam"
cmp "$tmp/fragment/depth.tsv.gz" "$tmp/fragment-parallel/depth.tsv.gz"

"$binary" -p "$fixture/target.bed" -o "$tmp/fragment-conflict" --flank 0 \
    --fragment-mode "$fixture/fragment_conflict.sam"
python3 - "$tmp/fragment-conflict/depth.tsv.gz" "$tmp/fragment-conflict/coverage.report.json" <<'PY'
import gzip
import json
import sys
rows = [line.rstrip().split("\t") for line in gzip.open(sys.argv[1], "rt")
        if not line.startswith("#")]
# The reciprocal mate position is deliberately inconsistent, so the overlap
# is retained and the pair is diagnosed rather than clipped.
assert rows[2][2:] == ["2", "2", "2"], rows
report = json.load(open(sys.argv[2], encoding="utf-8"))
assert report["fragment"]["paired"] == 0
assert report["fragment"]["ambiguous"] == 1
assert report["fragment"]["unmatched"] == 2
PY

"$binary" -p "$fixture/deletion.bed" -o "$tmp/deletion" "$fixture/deletion.sam"
python3 - "$tmp/deletion/depth.tsv.gz" <<'PY'
import gzip
import sys
rows = [line.rstrip().split("\t") for line in gzip.open(sys.argv[1], "rt")
        if not line.startswith("#")]
assert rows[5][2:5] == ["0", "0", "1"], rows[5]
PY

"$binary" -p "$fixture/deletion.bed" -o "$tmp/fragment-deletion" --fragment-mode \
    "$fixture/fragment_deletion.sam"
python3 - "$tmp/fragment-deletion/depth.tsv.gz" <<'PY'
import gzip
import sys
rows = [line.rstrip().split("\t") for line in gzip.open(sys.argv[1], "rt")
        if not line.startswith("#")]
# The M/D overlap is counted once in both the raw and D-inclusive unions.
assert all(row[2:] == ["1", "1", "1"] for row in rows), rows
PY

"$binary" -p "$fixture/empty.bed" -o "$tmp/empty" --summary-only \
    --compute-threads 2 "$fixture/mini.sam"
python3 - "$tmp/empty/coverage.report.json" <<'PY'
import json
import sys
report = json.load(open(sys.argv[1], encoding="utf-8"))
assert report["target"]["region_length"] == 0
assert report["target"]["average_depth"] == 0
PY

"$binary" -p "$fixture/overlap.bed" -o "$tmp/overlap" --flank 2 "$fixture/mini.sam"
python3 - "$tmp/overlap/region.tsv.gz" <<'PY'
import gzip
import sys
rows = [line.rstrip().split("\t") for line in gzip.open(sys.argv[1], "rt")
        if not line.startswith("#")]
assert len(rows) == 1 and rows[0][1:3] == ["2", "12"], rows
PY
printf 'chr1\t2\t3\nchr1\t15\t16\n' > "$tmp/large-flank.bed"
"$binary" -p "$tmp/large-flank.bed" -o "$tmp/large-flank" --flank 20 \
    --summary-only "$fixture/mini.sam"
python3 - "$tmp/large-flank/coverage.report.json" <<'PY'
import json
import sys
report = json.load(open(sys.argv[1], encoding="utf-8"))
# Both expanded windows cover the full chromosome; flank must still be the
# chromosome minus the merged target union, not the full expansion itself.
assert report["target"]["region_length"] == 2
assert report["flank"]["region_length"] == 18
PY
printf 'chr1\t3\t12\n' > "$tmp/one-based.bed"
"$binary" -1 -p "$tmp/one-based.bed" -o "$tmp/one-based" --flank 2 "$fixture/mini.sam"
python3 - "$tmp/one-based/region.tsv.gz" <<'PY'
import gzip
import sys
rows = [line.rstrip().split("\t") for line in gzip.open(sys.argv[1], "rt")
        if not line.startswith("#")]
assert len(rows) == 1 and rows[0][1:3] == ["2", "12"], rows
PY
"$binary" -p "$fixture/target.bed" -o "$tmp/maxdepth" --flank 2 --maxdepth 1 "$fixture/mini.sam"
python3 - "$tmp/maxdepth/cumu.plot" <<'PY'
import sys
for line in open(sys.argv[1], encoding="utf-8"):
    if line.strip():
        assert int(line.split("\t", 1)[0]) <= 1, line
PY
quoted="$tmp/input\"quoted.sam"
cp "$fixture/mini.sam" "$quoted"
"$binary" -p "$fixture/target.bed" -o "$tmp/quoted" "$quoted"
python3 - "$tmp/quoted/coverage.report.json" <<'PY'
import json
import sys
with open(sys.argv[1], encoding="utf-8") as handle:
    report = json.load(handle)
assert 'input"quoted.sam' in report["files"][0]
PY

if "$binary" -p "$fixture/bad.bed" -o "$tmp/bad" "$fixture/mini.sam" >/dev/null 2>&1; then
    echo "malformed BED unexpectedly succeeded" >&2
    exit 1
fi
if "$binary" -p "$fixture/target.bed" -o "$tmp/unsorted" "$fixture/unsorted.sam" >/dev/null 2>&1; then
    echo "unsorted input unexpectedly succeeded" >&2
    exit 1
fi
if "$binary" -p "$fixture/target.bed" -o "$tmp/invalid" --mapthres 256 "$fixture/mini.sam" >/dev/null 2>&1; then
    echo "invalid mapQ unexpectedly succeeded" >&2
    exit 1
fi
if "$binary" -p "$fixture/target.bed" -o "$tmp/invalid-compute" --compute-threads 0 "$fixture/mini.sam" >/dev/null 2>&1; then
    echo "invalid compute thread count unexpectedly succeeded" >&2
    exit 1
fi
if "$binary" -p "$fixture/target.bed" -o "$tmp/invalid-compute2" --compute-threads 257 "$fixture/mini.sam" >/dev/null 2>&1; then
    echo "oversized compute thread count unexpectedly succeeded" >&2
    exit 1
fi
if "$binary" -p "$fixture/target.bed" -o "$tmp/invalid-ratio" --depthratio nan "$fixture/mini.sam" >/dev/null 2>&1; then
    echo "NaN depth ratio unexpectedly succeeded" >&2
    exit 1
fi
if "$binary" --compute-threads 4 -p "$fixture/target.bed" -o "$tmp/invalid-cigar" \
    "$fixture/bad_cigar.sam" >/dev/null 2>&1; then
    echo "invalid CIGAR unexpectedly succeeded" >&2
    exit 1
fi
if "$binary" -p "$fixture/target.bed" -o "$tmp/invalid-list" --cutoffdepth 1, "$fixture/mini.sam" >/dev/null 2>&1; then
    echo "empty cutoff item unexpectedly succeeded" >&2
    exit 1
fi
if "$binary" -p "$fixture/target.bed" -o "$tmp/invalid-stdin" - - </dev/null >/dev/null 2>&1; then
    echo "duplicate stdin unexpectedly succeeded" >&2
    exit 1
fi

python3 "$root/tests/oracle.py" "$binary"

echo "xamdst integration tests passed"
