# DNA performance comparison

The v3.2.0 build used for the DNA-path comparison was compared with 3.1.0 (`af559ba`) in WSL2 Ubuntu 24.04, GCC 13.3.0, using `-O2 -g -std=c99 -Wall -Wextra -Wpedantic -Wformat=2 -Werror`. Each main case used five million synthetic 100 bp records, fixed CPUs 12–15, and interleaved baseline/candidate runs. `panel` uses a 40 kb reference and two targeted regions; `broad` uses a 2 Mb reference and two regions.

| Dataset | Threads | Summary only | Wall change | CPU change |
| --- | ---: | --- | ---: | ---: |
| panel | 1 | False | +48.42% | -2.91% |
| panel | 1 | True | -9.17% | -8.70% |
| panel | 4 | False | +2.05% | +0.64% |
| panel | 4 | True | -6.39% | -9.36% |
| broad | 1 | False | -1.81% | -2.53% |
| broad | 1 | True | -6.26% | -6.60% |
| broad | 4 | False | +5.34% | +6.27% |
| broad | 4 | True | +2.55% | +2.04% |

The wide-reference parallel cases were repeated for ten more rounds with balanced trial order. Median CPU changes were:

| Dataset | Threads | Summary only | Baseline CPU | v3.2.0 CPU | Change |
| --- | ---: | --- | ---: | ---: | ---: |
| broad | 4 | False | 4.141s | 4.234s | +2.24% |
| broad | 4 | True | 3.892s | 4.004s | +2.87% |

All depth rows, read counts, summary fields and output files matched after normalizing version fields and the documented `>` to `>=` cumulative-plot change. Serial results were faster; panel parallel CPU time was within 1%; broad-reference parallel CPU time remained about 2–3% higher. The parallel reducer stages CIGAR overlaps as per-record alignment deltas and then applies them in coordinate order; that path is a plausible contributor to the small broad-reference cost. No representative production alignment file was available for this comparison.

Raw measurements: [eight-case comparison](dna-performance.json) and [ten-round broad parallel repeat](dna-performance-broad-repeat.json).
