# xamdst tests

`run_tests.sh` exercises the SAM, BED (including overlap and 1-based input),
CIGAR, multi-input, JSON, strict-parameter, summary-only, fragment-mode and
atomic-output paths without
requiring `samtools`. CRAM tests are run in Linux CI with HTSlib tooling.
`benchmark.sh` accepts `XAMDST_BENCH_INPUT`, `XAMDST_BENCH_BED`,
`XAMDST_BENCH_OUT`, `XAMDST_BENCH_COMPUTE_THREADS` and
`XAMDST_BENCH_SUMMARY_ONLY`/`XAMDST_BENCH_BASELINE` for representative production data and reports the
median of five wall-clock runs with CPU/RSS/output-size measurements.
`oracle.py` generates a deterministic mixed-CIGAR case and checks raw, rmdup
and deletion-inclusive per-base depth. `generate_panel.py` creates a larger
coordinate-sorted high-depth SAM/BED pair for reproducible benchmarks.

`rna_oracle.py` generates STAR-default and AllBestScore fixtures, including
cross-reference primary alignments, paired read ends, missing/malformed tags,
mixed CIGAR operations and overlapping GTF/GFF2 features. A Python per-base
oracle checks depth, distribution, read counts and strand quadrants. The suite
also checks gzip/truncated annotations, exact serial/parallel/spill output,
summary file sets, multi-input identity, and rollback on parsing and disk-write
failures. Forced spill tests use the optional `--rna-dedup-mem 128` limit.
