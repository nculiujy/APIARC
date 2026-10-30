#!/usr/bin/env python3
"""Derive a TSS/gene-body BED from a GENCODE GTF (shared helper).

The enhancer module needs the genomic coordinates of genes. The previous version
passed a promoter-peak BED as ``genefile``, so distances were measured to peaks
rather than to transcription start sites. This helper produces the correct file
once per assembly.
"""

from __future__ import annotations

import argparse
import gzip
import re


ATTRIBUTE = re.compile(r'(\w+) "([^"]+)"')


def open_maybe_gzip(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def main() -> None:
    parser = argparse.ArgumentParser(description="GENCODE GTF -> TSS BED")
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--feature", default="gene")
    parser.add_argument("--extend", type=int, default=0,
                        help="extend the TSS window by N bp on each side")
    args = parser.parse_args()

    seen: set[str] = set()
    with open_maybe_gzip(args.gtf) as handle, open(args.out, "w") as out:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != args.feature:
                continue
            attributes = dict(ATTRIBUTE.findall(fields[8]))
            gene_id = attributes.get("gene_id", "")
            symbol = attributes.get("gene_name") or attributes.get("gene_id", "")
            if not gene_id or gene_id in seen:
                continue
            seen.add(gene_id)

            start, end = int(fields[3]), int(fields[4])
            if fields[6] == "+":
                tss_start, tss_end = start - 1, start
            else:
                tss_start, tss_end = end - 1, end
            if args.extend:
                tss_start = max(0, tss_start - args.extend)
                tss_end = tss_end + args.extend
            out.write(
                f"{fields[0]}\t{tss_start}\t{tss_end}\t{symbol}\t{gene_id}\t{fields[6]}\n"
            )

    print(f"[tss] wrote {len(seen)} genes to {args.out}")


if __name__ == "__main__":
    main()
