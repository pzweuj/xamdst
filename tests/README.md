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
