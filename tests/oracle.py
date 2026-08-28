#!/usr/bin/env python3
"""Slow, dependency-free oracle for the three 3.1 depth series.

The fixture is deterministic (the seed is fixed) but exercises a mixture of
M/= /X/D/N/I/S CIGAR operations and duplicate, QC-fail and low-mapQ flags.
It intentionally uses a small chromosome so that a straightforward
per-base implementation can serve as a regression oracle for the optimized
target-local difference arrays.
"""

from __future__ import annotations

import gzip
import random
import subprocess
import sys
import tempfile
from pathlib import Path


CHROMOSOME = "chr1"
CHROMOSOME_LENGTH = 160
TARGET_START = 8
TARGET_END = 145
MAPQ_CUTOFF = 20


def parse_cigar(cigar: str):
    number = ""
    for char in cigar:
        if char.isdigit():
            number += char
            continue
        if not number:
            raise ValueError(f"invalid CIGAR: {cigar}")
        yield int(number), char
        number = ""
    if number:
        raise ValueError(f"invalid CIGAR: {cigar}")


def cigar_lengths(cigar: str):
    query = 0
    reference = 0
    for length, op in parse_cigar(cigar):
        if op in "M=X":
            query += length
            reference += length
        elif op in "IS":
            query += length
        elif op in "DN":
            reference += length
        elif op in "HP":
            pass
        else:
            raise ValueError(f"unsupported CIGAR operation: {op}")
    return query, reference


def expected_depth(records):
    raw = [0] * (TARGET_END - TARGET_START)
    rmdup = [0] * (TARGET_END - TARGET_START)
    coverage = [0] * (TARGET_END - TARGET_START)
    for start, mapq, flag, cigar in records:
        duplicate = bool(flag & 0x400)
        qc_fail = bool(flag & 0x200)
        clean = not duplicate and not qc_fail and mapq >= MAPQ_CUTOFF
        reference = start
        for length, op in parse_cigar(cigar):
            if op in "M=X":
                for position in range(reference, reference + length):
                    if TARGET_START <= position < TARGET_END:
                        index = position - TARGET_START
                        raw[index] += 1
                        if clean:
                            rmdup[index] += 1
                        coverage[index] += 1
                reference += length
            elif op == "D":
                for position in range(reference, reference + length):
                    if TARGET_START <= position < TARGET_END:
                        coverage[position - TARGET_START] += 1
                reference += length
            elif op == "N":
                reference += length
            elif op in "ISHP":
                pass
            else:
                raise ValueError(op)
    return raw, rmdup, coverage


def make_records():
    rng = random.Random(310)
    patterns = [
        "5M",
        "3=2X",
        "2S4M1I3M",
        "4M2D3M",
        "3M5N4M",
        "1S2M1I2=1X1S",
    ]
    records = []
    for index in range(36):
        start = 2 + index * 4
        cigar = patterns[index % len(patterns)]
        query_length, reference_length = cigar_lengths(cigar)
        if start + reference_length >= CHROMOSOME_LENGTH:
            break
        mapq = (10, 20, 35, 60)[index % 4]
        flag = (0x400 if index % 7 == 0 else 0) | (0x200 if index % 11 == 0 else 0)
        sequence = "".join(rng.choice("ACGT") for _ in range(query_length))
        quality = "I" * query_length
        records.append((start, mapq, flag, cigar, sequence, quality))
    return records


def write_sam(path: Path, records):
    with path.open("w", encoding="ascii", newline="\n") as handle:
        handle.write("@HD\tVN:1.6\tSO:coordinate\n")
        handle.write(f"@SQ\tSN:{CHROMOSOME}\tLN:{CHROMOSOME_LENGTH}\n")
        for index, (start, mapq, flag, cigar, sequence, quality) in enumerate(records):
            handle.write(
                f"oracle{index:03d}\t{flag}\t{CHROMOSOME}\t{start + 1}\t{mapq}\t{cigar}"
                f"\t*\t0\t0\t{sequence}\t{quality}\n"
            )


def read_depth(path: Path):
    rows = []
    with gzip.open(path, "rt", encoding="ascii") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            rows.append(tuple(int(value) for value in fields[2:5]))
    return rows


def main(binary: str) -> int:
    records = make_records()
    expected = expected_depth([(r[0], r[1], r[2], r[3]) for r in records])
    expected_rows = list(zip(*expected))
    with tempfile.TemporaryDirectory(prefix="xamdst-oracle-") as temporary:
        root = Path(temporary)
        sam = root / "oracle.sam"
        bed = root / "oracle.bed"
        output = root / "output"
        write_sam(sam, records)
        bed.write_text(f"{CHROMOSOME}\t{TARGET_START}\t{TARGET_END}\n", encoding="ascii")
        subprocess.run(
            [
                binary,
                "--compute-threads",
                "4",
                "--mapthres",
                str(MAPQ_CUTOFF),
                "--flank",
                "0",
                "-p",
                str(bed),
                "-o",
                str(output),
                str(sam),
            ],
            check=True,
        )
        actual = read_depth(output / "depth.tsv.gz")
        if actual != expected_rows:
            for index, (want, got) in enumerate(zip(expected_rows, actual)):
                if want != got:
                    raise AssertionError(
                        f"depth mismatch at target offset {index}: expected {want}, got {got}"
                    )
            raise AssertionError(f"depth row count mismatch: expected {len(expected_rows)}, got {len(actual)}")
    print("xamdst oracle test passed")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(f"usage: {sys.argv[0]} path/to/xamdst")
    raise SystemExit(main(sys.argv[1]))
