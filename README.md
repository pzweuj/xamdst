# xamdst

xamdst is a streaming depth and coverage report generator for coordinate-sorted
BAM, CRAM and SAM files over BED target regions. Version 3.1.0 keeps stdin,
multiple inputs, CRAM references, HTSlib I/O threads, JSON and target-BAM export,
while making interval semantics and filtering explicit.

## Dependencies

- HTSlib >= 1.13
- zlib
- a POSIX C99 toolchain (Linux or macOS)
- Python 3 (for `make test` and the oracle/benchmark helpers)

Debian/Ubuntu:

```bash
sudo apt-get install build-essential pkg-config libhts-dev zlib1g-dev python3 time
```

macOS:

```bash
brew install htslib
```

## Build

```bash
make
make test
```

If `pkg-config` cannot find HTSlib, provide its installation/source directory:
`make HTSLIB_DIR=/path/to/htslib`.

The project does not provide a native Windows build; use Docker, WSL or a
POSIX-compatible environment.

## Usage

```bash
xamdst -p targets.bed -o results input.bam
xamdst -p targets.bed -o results --threads 8 input.bam
xamdst -p targets.bed -o results input.cram --reference reference.fa
samtools view input.bam -u | xamdst -p targets.bed -o results -
```

All input files must have the same reference dictionary and coordinate sort
order. Multiple files are merged by `(reference, position)` before counting;
the exported target BAM is coordinate ordered as well. Only one `-` stdin input
is allowed. The BAM output path must be a `.bam` file distinct from every input
and from the report files.

## Options

| Option | Default | Description |
| --- | --- | --- |
| `-p, --bed FILE` | required | target regions; 0-based, half-open BED |
| `-o, --outdir DIR` | required | output directory, created if needed |
| `-T, --reference FILE` | — | reference FASTA, required for CRAM |
| `-f, --flank N` | `200` | outside flank size; target itself is excluded |
| `-q, --mapthres N` | `20` | clean-depth mapQ cutoff, `0..255` |
| `--maxdepth N` | `0` | cap rows written to `cumu.plot`; `0` means no cap |
| `--cutoffdepth D1,D2,...` | — | up to 10 positive `>=` depth cutoffs |
| `--depthratio R1,R2,...` | `0.2,0.5` | up to 10 ratios in `(0,1]`; report depth `> R × average` |
| `--isize N` | `2000` | include positive insert sizes below N |
| `--uncover N` | `5` | write positions with raw depth `< N` to `uncover.bed` |
| `--bamout FILE` | — | export primary mapped reads intersecting target |
| `--threads N` | `0` | one shared HTSlib compression/decompression thread pool across inputs and BAMout |
| `-1` | off | interpret input coordinates as 1-based inclusive |
| `-h, --help` | — | show help |
| `--compute-threads N` | `1` | CIGAR/coverage workers (`1..256`) |
| `--fragment-mode` | off | subtract actual overlaps of matched proper pairs |
| `--summary-only` | off | omit `depth.tsv.gz` while retaining other reports |
| `-v, --version` | — | show `3.1.0` |

Blank/comment-only BED files are accepted and produce empty, zero-depth target
reports; malformed intervals are rejected with their line number.

## Depth definitions

- `raw`: primary mapped `M`, `=`, and `X` bases, including duplicate, QC-fail
  and low-mapQ reads.
- `rmdup`: raw depth excluding duplicate and QC-fail reads and requiring mapQ
  at least `--mapthres`.
- `coverage_with_deletions`: raw depth plus `D` CIGAR operations. `N` advances the reference
  coordinate but does not add depth.

Secondary, supplementary and unmapped records do not contribute coverage.
Target read counts and `--bamout` include a record when one of its `M/=/X/D`
segments intersects a target interval.
Read-level counters are based on primary records; unmapped primaries contribute
to read totals but never to depth.

## Outputs

The default output directory contains exactly these eight files (summary-only
mode intentionally contains seven):

| File | Content |
| --- | --- |
| `coverage.report` | text summary of read, target and flank statistics |
| `coverage.report.json` | same summary with `schema_version: "3.1"`; `depth_output` records whether depth was emitted |
| `cumu.plot` | target raw-depth distribution and cumulative counts |
| `insert.plot` | insert-size distribution and cumulative counts |
| `chromosome.report` | per-reference target coverage |
| `region.tsv.gz` | merged target interval raw mean/median/coverage plus explicitly named D-inclusive columns |
| `depth.tsv.gz` | chromosome, 1-based position, raw/rmdup/coverage-with-deletions depth |
| `uncover.bed` | 0-based half-open contiguous low-coverage spans |

Outputs are written to temporary files and atomically installed only after a
successful run. Existing results remain intact if input parsing or writing
fails.

## 3.0 → 3.1 migration

The three historical names `depth_distribution.plot`, `insertsize.plot` and
`chromosomes.report` are no longer produced. Raw `M/=/X` depth is now the
default series for `target_data*`, `average_depth`, `cumu.plot`, chromosome
columns, coverage percentages, ratios and `uncover.bed`. D-inclusive values
are available as `coverage_with_deletions`/`average_depth_with_deletions` and
explicit report columns. The flank section excludes target bases, CIGAR `N`
is handled correctly, and malformed parameters/BED records are errors instead
of undefined behavior. `--summary-only` removes a previous `depth.tsv.gz`
transactionally; a failed run restores it.

| 3.0 field | 3.1 field |
| --- | --- |
| `target_data_mb`, `average_depth`, `coverage` | raw `M/=/X` values |
| implicit `coverage` with `D` | `target_data_with_deletions_mb`, `average_depth_with_deletions`, `coverage_with_deletions` |
| `Coverage(FIX)` region column | named D-inclusive mean/median/coverage columns |
| per-input `--threads` pools | one shared I/O budget (implementation may use HTSlib pool) |

For high-depth panels, use `--compute-threads 4` and `--summary-only` when a
per-base depth file is not needed.

`--compute-threads` validates and segments CIGARs in a bounded persistent worker
pool. Results are always reduced in input coordinate order, so changing the
worker count does not change any report or BAMout bytes. `--threads` is kept
separate and is shared by all HTSlib input/output handles.

Run `make test` for the fixture-based integration suite. Set
`XAMDST_BENCH_INPUT` and `XAMDST_BENCH_BED` before `make benchmark` to measure a
representative production workload (`/usr/bin/time` is required by the
benchmark target; Debian/Ubuntu provide it in the `time` package).

For a reproducible high-depth panel workload, generate data outside the
repository and pass its paths to the benchmark:

```bash
python3 tests/generate_panel.py --outdir /tmp/xamdst-panel --depth 500
XAMDST_BENCH_INPUT=/tmp/xamdst-panel/panel.sam \
XAMDST_BENCH_BED=/tmp/xamdst-panel/panel.bed make benchmark
```
