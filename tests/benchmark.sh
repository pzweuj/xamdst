#!/bin/sh
set -eu

if [ "$#" -ne 1 ]; then
    echo "usage: $0 path/to/xamdst" >&2
    exit 2
fi
binary=$1
if [ ! -x /usr/bin/time ]; then
    echo "benchmark requires /usr/bin/time (install the 'time' package)" >&2
    exit 2
fi
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
input=${XAMDST_BENCH_INPUT:-"$root/tests/fixtures/mini.sam"}
bed=${XAMDST_BENCH_BED:-"$root/tests/fixtures/target.bed"}
compute_threads=${XAMDST_BENCH_COMPUTE_THREADS:-1}
summary_only=${XAMDST_BENCH_SUMMARY_ONLY:-0}
baseline=${XAMDST_BENCH_BASELINE:-}
if [ -n "${XAMDST_BENCH_OUT:-}" ]; then
    out=$XAMDST_BENCH_OUT
    cleanup=0
else
    out=$(mktemp -d "${TMPDIR:-/tmp}/xamdst-bench.XXXXXX")
    cleanup=1
fi
trap '[ "$cleanup" -eq 1 ] && rm -rf "$out"' EXIT HUP INT TERM

mkdir -p "$out"

echo "benchmark input: $input"
echo "benchmark bed:   $bed"
echo "compute threads: $compute_threads"
echo "summary only:    $summary_only"
if [ -n "$baseline" ]; then
    echo "baseline binary: $baseline"
fi
: > "$out/timing.log"
for run in 1 2 3 4 5; do
    run_out="$out/run-$run"
    rm -rf "$run_out"
    summary_arg=
    if [ "$summary_only" = 1 ]; then summary_arg=--summary-only; fi
    /usr/bin/time -f "run=$run wall=%e user=%U sys=%S rss_kb=%M" \
        "$binary" --compute-threads "$compute_threads" $summary_arg -p "$bed" -o "$run_out" "$input" >/dev/null 2>>"$out/timing.log"
    du -sk "$run_out" | awk -v run="$run" '{print "run=" run " output_kb=" $1}' >>"$out/timing.log"
done
cat "$out/timing.log"
median_metric() {
    metric=$1
    awk -v metric="$metric" '
        /^run=[0-9]+/ {
            for (i = 1; i <= NF; ++i) {
                split($i, pair, "=")
                if (pair[1] == metric) print pair[2]
            }
        }
    ' "$out/timing.log" | sort -n | sed -n '3p'
}
for metric in wall user sys rss_kb output_kb; do
    value=$(median_metric "$metric")
    if [ -n "$value" ]; then
        echo "median_${metric}=$value"
    fi
done

if [ -n "$baseline" ] && [ -x "$baseline" ]; then
    echo "running five baseline iterations"
    for run in 1 2 3 4 5; do
        baseline_out="$out/baseline-$run"
        rm -rf "$baseline_out"
        /usr/bin/time -f "baseline_run=$run wall=%e user=%U sys=%S rss_kb=%M" \
            "$baseline" -p "$bed" -o "$baseline_out" "$input" >/dev/null 2>>"$out/timing.log"
        du -sk "$baseline_out" | awk -v run="$run" '{print "baseline_run=" run " output_kb=" $1}' >>"$out/timing.log"
    done
    baseline_median_metric() {
        metric=$1
        awk -v metric="$metric" '
            /^baseline_run=[0-9]+/ {
                for (i = 1; i <= NF; ++i) {
                    split($i, pair, "=")
                    if (pair[1] == metric) print pair[2]
                }
            }
        ' "$out/timing.log" | sort -n | sed -n '3p'
    }
    for metric in wall user sys rss_kb output_kb; do
        value=$(baseline_median_metric "$metric")
        if [ -n "$value" ]; then
            echo "baseline_median_${metric}=$value"
        fi
    done
fi
