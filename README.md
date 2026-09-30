# xamdst

xamdst is a streaming depth and coverage report generator for coordinate-sorted
BAM, CRAM and SAM files over BED target regions. Version 3.2.0 adds optional RNA
QC (unique-mapping depth, splicing, annotation distribution and strand inference),
and keeps stdin,
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

## Recommend

```
docker pull ghcr.io/pzweuj/mapping:2026Aug
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
| `--rna` | off | unique-mapping RNA depth and QC; incompatible with `--fragment-mode` |
| `--annotation FILE` | — | optional GTF/GFF2, plain or gzip; requires `--rna` |
| `--rna-dedup-mem N` | `67108864` | RNA multimap sorting buffer budget in bytes, minimum 128 |
| `--max-region-mem N` | `0` | per-region coverage buffer limit in bytes; 0 means unlimited |
| `-v, --version` | — | show `3.2.0` |

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
| `coverage.report.json` | same summary with `schema_version: "3.2"`; `depth_output` records whether depth was emitted |
| `cumu.plot` | target raw-depth distribution and cumulative counts |
| `insert.plot` | insert-size distribution and cumulative counts |
| `chromosome.report` | per-reference target coverage |
| `region.tsv.gz` | merged target interval raw mean/median/coverage plus explicitly named D-inclusive columns |
| `depth.tsv.gz` | chromosome, 1-based position, raw/rmdup/coverage-with-deletions depth |
| `uncover.bed` | 0-based half-open contiguous low-coverage spans |

Outputs are written to temporary files and atomically installed only after a
successful run. Existing results remain intact if input parsing or writing
fails.

## RNA usage and counting

```bash
xamdst --rna -p targets.bed -o rna-results input.bam
xamdst --rna --annotation genes.gtf.gz -p targets.bed -o rna-results \
  --compute-threads 4 --summary-only input.bam
```

The annotation is optional. RNA mode produces ten report files, or nine with
`--summary-only`: the DNA files above plus `splice.tsv.gz` and `distribution.tsv`.
The text report adds `[RNA]` lines and JSON adds an `rna` object only in RNA mode.

RNA reads are **read ends**: each mate is counted separately. `total.raw_reads`
and the existing total statistics still count primary records (neither secondary
nor supplementary), including unmapped primary records. Multiple primary
alignments can therefore inflate the total primary-record counts without
inflating the RNA multimapping read count.

Only mapped primary records with `NH=1` contribute RNA depth, target/flank read
counts and target BAM export. A missing NH uses the compatibility fallback of
treating the record as unique; every run with missing NH emits a warning and
reports `nh_missing_records`. NH and HI, when present on mapped primary records,
must be positive integer tags; HI must not exceed NH when NH is present. Missing
HI is accepted and never changes which alignment is primary.

For `NH>1`, `multimap_records` counts primary records, while `multimap_reads`
counts distinct `(input file ordinal, QNAME, read1/read2 flags)` keys over the
entire run, including across chromosomes. Supplying the same input twice counts
it twice. This supports both STAR's default one-primary output and
`--outSAMprimaryFlag AllBestScore`, which may emit multiple primaries; it does
not assume `HI=1` identifies a primary. See the
[STAR manual](https://github.com/alexdobin/STAR/blob/master/extras/doc-latex/STARmanual.tex).

The dedup buffer defaults to 64 MiB. At its budget it spills sorted partitions
into the output directory and performs exact external merges; there is no
external `sort` dependency. Memory also includes pointer capacity, merge buffers
and at most one oversized key. `--rna-dedup-mem` is optional and can lower the
budget. Handled failures remove these temporary files and preserve existing
reports. An abrupt process kill can leave `.xamdst-rna-dedup-*` scratch files.

RNA depth retains the raw/rmdup/deletion-inclusive definitions above after the
NH filter. Splice counts use unique primary records: `spliced_reads` counts a
read end once if its CIGAR contains `N`; the intron histogram counts every `N`
operation and its reference length. Spliced read1 records are excluded from the
insert-size histogram because their TLEN spans introns.

## RNA annotation and strand evidence

Annotation coordinates are 1-based inclusive and reference names must match the
BAM header exactly. Out-of-range or malformed features are errors; references
absent from the header are skipped with a warning. An annotation with no matching
exon/gene/transcript features is rejected. GTF `gene_id` or GFF2 `gene` groups
features; `transcript_id` or GFF2 `Transcript` is a fallback for exon-only files.
Grouped exons use `gene_id`/`transcript_id` or GFF2 `gene`/`Transcript`; an exon
with `.` attributes is retained as an exon interval but cannot extend an
inferred transcript body. GFF3 `ID=...;Parent=...` is not supported.

Exons are merged over all genes/transcripts. Gene bodies use explicit gene or
transcript extents and the span of grouped exons. Introns are the union of gene
bodies minus the global exon union. Only `M/=/X` bases on unique primary records
contribute to the exonic/intronic/intergenic distribution, across the entire
alignment rather than only the BED targets. `D`, `N`, insertions and clips add no
distribution bases. Each read with at least one aligned `M/=/X` base takes the
class with the most bases; ties prefer exon, then intron, then intergenic.
References without annotated features are intergenic. Region-length denominators
cover the full BAM reference dictionary (`annotated_span`), including references
without annotation.

Strand inference uses unique reads overlapping exons with `+` or `-` annotation.
A read overlapping both directions anywhere in its aligned bases is excluded
and counted in `ambiguous_reads`. Unstranded exon annotations add no evidence.
Read1 (or a single-end read) antisense plus read2 sense supports `fr-firststrand`;
the reverse supports `fr-secondstrand`. At least 1,000 informative read ends are
required. Fractions above 0.7 call the corresponding protocol; a firststrand
fraction in [0.4, 0.6] calls `unstranded`; the remainder is `mixed`.
Below the sample threshold the result is `insufficient_data`, while observed
fractions are still reported. Without annotation the evidence is `none` and
inference is always `insufficient_data`. XS tag counts remain available as raw
diagnostics, independently of annotation-based inference.

| RNA output/field | Meaning |
| --- | --- |
| `splice.tsv.gz` | `IntronLength`, `Count`, `Fraction` for observed CIGAR `N` lengths |
| `distribution.tsv` | class, read count/fraction, aligned base count/fraction, reference region length; comment-only data section without annotation |
| `rna.unique_mapping_reads`, `unique_mapping_rate` | unique/fallback records and fraction of unique plus distinct multimapping read ends |
| `rna.multimap_records`, `multimap_reads`, `multimap_read_fraction` | primary multimap records, exact distinct ends and their fraction of mapped RNA ends |
| `rna.nh_missing_records` | primary mapped records using the unique fallback |
| `rna.spliced_reads`, `spliced_read_fraction`, `introns` | spliced ends/fraction of unique ends; N count, mean, median, P5/P25/P75/P95 lengths |
| `rna.strand.evidence`, `effective_reads`, `ambiguous_reads` | evidence source, usable ends, and opposite-strand overlap exclusions |
| `rna.strand.annotation_read{1,2}_{antisense,sense}` | four annotation-derived strand quadrants |
| `rna.strand.xs_reads`, `q_read{1,2}_{antisense,sense}` | raw valid XS counts and quadrants |
| `rna.strand.fr_firststrand_fraction`, `fr_secondstrand_fraction`, `inference` | observed evidence fractions and inferred protocol |
| `rna.strand.single_end_observed` | whether informative annotation evidence includes single-end reads |
| `rna.distribution` | per-class reads/bases/fractions, exon/intron union lengths and header span; null without annotation |

RNA JSON fractions and TSV fractions are in [0,1]; text report rates are
percentages except the explicitly named strand fractions. Existing DNA JSON
percentage fields retain their 3.1 conventions.

## 3.1 → 3.2 migration

DNA invocations keep the eight-file set (seven with `--summary-only`), primary
record counts and depth definitions. JSON `schema_version` changes from `3.1`
to `3.2` in both modes. Consumers enforcing a schema version must accept `3.2`;
the optional `rna` object and two extra files appear only with `--rna`.
No new required command-line argument is introduced. RNA mode and annotation
must be selected explicitly. Placed unmapped records retain coordinate order
and contribute only to primary read totals. Summary-only output removes stale
depth output transactionally in both modes.
Reusing an RNA output directory in DNA mode also removes the two stale RNA
reports transactionally. All ten report names are reserved against input-path
collisions in either mode.

The existing 3.2 cumulative-plot correction changes columns 4/5 of `cumu.plot`
and `insert.plot` from strictly greater than the row value to greater than or
equal to it. Histogram bin counts and depth values are unchanged. Consumers
using the former strict convention should subtract the row's count/fraction.
The existing insert-size correction includes only positive TLEN below `--isize`
on clean, proper-pair read1 records whose mate maps to the same reference;
3.1 included any positive TLEN below the cutoff. This affects `insert.plot`
and insert-size summaries, while primary read totals are unchanged.

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
