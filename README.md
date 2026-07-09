# xamdst

A fast depth-coverage statistic tool for **BAM/CRAM/SAM** files over target regions.
A hard fork of [bamdst](https://github.com/shiquan/bamdst), refactored on
[HTSlib](https://github.com/samtools/htslib) — adding CRAM support, multi-threaded
I/O, and a JSON report.

Input files must be **coordinate-sorted**. The probe file and output directory are
mandatory.

## Features

- BAM / CRAM / SAM input, single file or multiple files merged
- Multi-threaded compression & decompression (HTSlib native)
- Per-position, per-region, per-chromosome and cumulative depth statistics
- `raw`, `rmdup` (dedup + primary + mapQ filter) and `coverage` (CIGAR-aware) depths
- JSON report for programmatic access
- `depth.tsv.gz` and `region.tsv.gz` are bgzipped and tabix-ready

## Dependencies

- **htslib**, **zlib**, **pthread**

```bash
# Debian/Ubuntu
apt-get install libhts-dev zlib1g-dev
# macOS (Homebrew)
brew install htslib
```

## Install

```bash
git clone https://github.com/pzweuj/xamdst
cd xamdst
make
```

Or use the Docker image:

```bash
docker pull ghcr.io/pzweuj/mapping:2026Jan
```

## Usage

```bash
# basic
xamdst -p probe.bed -o ./ in1.bam

# multi-threaded
xamdst -p probe.bed -o ./ --threads 8 in1.bam

# CRAM (reference required)
xamdst -p probe.bed -o ./ in1.cram -T ref.fa

# streaming from stdin
samtools view in1.bam -u | xamdst -p probe.bed -o ./ -
```

## Options

| Option | Default | Description |
| --- | --- | --- |
| `-p, --bed FILE` | — | target / probe regions (merged first). **Required** |
| `-o, --outdir DIR` | — | output directory (created if missing). **Required** |
| `-T, --reference FILE` | — | reference FASTA, **required for CRAM** |
| `-f, --flank N` | 200 | flank bp to stat around each region |
| `-q, --mapthres N` | 20 | mapQ cutoff (≥ counted as quality reads) |
| `--maxdepth N` | 0 | cap depth in the cumulative distribution (0 = no cap) |
| `--cutoffdepth D1,D2,...` | 0 | report coverage at these depths (max 10) |
| `--depthratio R1,R2,...` | 0.2,0.5 | coverage at ratios of average depth |
| `--isize N` | 2000 | insert-size cutoff for the distribution |
| `--uncover N` | 5 | positions below this depth go to `uncover.bed` |
| `--bamout FILE` | — | export target reads to this BAM |
| `--threads N` | 0 | threads for BAM/CRAM I/O (0 = single-threaded) |
| `-1` | off | BED begin position is 1-based |
| `-h, --help` | | show help |
| `-v, --version` | | show version |

## Output

Eight files are written to the output directory:

| File | Content |
| --- | --- |
| `coverage.report` | coverage & read statistics for target and flank regions |
| `coverage.report.json` | same report in JSON |
| `cumu.plot` | depth cumulative distribution |
| `insert.plot` | inferred insert-size distribution |
| `chromosome.report` | per-chromosome coverage |
| `region.tsv.gz` | mean / median depth and coverage per region (bgzip + tabix) |
| `depth.tsv.gz` | per-position raw / rmdup / coverage depth (bgzip + tabix) |
| `uncover.bed` | poorly covered or uncovered regions |

### `depth.tsv.gz` columns

| Column | Meaning |
| --- | --- |
| chromosome | chromosome name |
| position | 1-based position |
| raw depth | unfiltered depth |
| rmdup depth | deduplicated, primary, mapQ ≥ cutoff (default 20) — close to `samtools depth` |
| coverage depth | raw depth plus CIGAR deletions (≥ raw depth) |

Coverage statistics in `coverage.report` use the **coverage depth**. Pass
`--depthratio` / `--cutoffdepth` for custom thresholds.

## Notes

- All input files must use the same reference and sort order.
- The first file's header is used for output.
- `region.tsv.gz` and `depth.tsv.gz` can be indexed with `tabix -p bed` / `tabix -p vcf`.

Homepage: <https://github.com/pzweuj/xamdst>
