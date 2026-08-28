#!/usr/bin/env python3
"""Generate a deterministic high-depth panel for local benchmarks.

The generated SAM is coordinate sorted and deliberately modest by default so
it is practical on a laptop.  Increase ``--depth``/``--length`` for a
production-like stress run; no generated data is committed to the repository.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--length", type=int, default=20000)
    parser.add_argument("--depth", type=int, default=500)
    parser.add_argument("--read-length", type=int, default=100)
    args = parser.parse_args()
    if args.length <= 0 or args.depth <= 0 or args.read_length <= 0:
        parser.error("length, depth and read-length must be positive")
    if args.read_length > args.length:
        parser.error("read-length must not exceed chromosome length")

    args.outdir.mkdir(parents=True, exist_ok=True)
    sam_path = args.outdir / "panel.sam"
    bed_path = args.outdir / "panel.bed"
    # A 20% stride keeps the file size manageable while producing the
    # requested depth over the target intervals.
    stride = max(1, args.read_length // 5)
    starts = range(0, args.length - args.read_length + 1, stride)
    with sam_path.open("w", encoding="ascii", newline="\n") as sam:
        sam.write("@HD\tVN:1.6\tSO:coordinate\n")
        sam.write(f"@SQ\tSN:chrPanel\tLN:{args.length}\n")
        sequence = "A" * args.read_length
        quality = "I" * args.read_length
        read_index = 0
        for start in starts:
            for replicate in range(args.depth):
                sam.write(
                    f"panel{read_index:09d}\t0\tchrPanel\t{start + 1}\t60\t"
                    f"{args.read_length}M\t*\t0\t0\t{sequence}\t{quality}\n"
                )
                read_index += 1
    # Two separated targets exercise both the window lookup and flank
    # subtraction paths while remaining easy to inspect.
    first_start = args.length // 20
    first_end = max(first_start + 1, args.length // 5)
    second_start = min(args.length - 1, max(first_end + 1, args.length * 3 // 5))
    second_end = min(args.length, max(second_start + 1, second_start + args.length // 10))
    bed_path.write_text(
        f"chrPanel\t{first_start}\t{first_end}\n"
        f"chrPanel\t{second_start}\t{second_end}\n",
        encoding="ascii",
    )
    print(sam_path)
    print(bed_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
