#!/usr/bin/env python3
import argparse
import gzip
import json
import sys
from typing import Dict, Iterable, List, Optional, Set, Tuple

def open_maybe_gz(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")

def ensure_chr(chrom: str) -> str:
    if chrom.startswith("chr"):
        return chrom
    return f"chr{chrom}"

def parse_gtf_attrs(attr_str: str) -> Dict[str, str]:
    """
    Ensembl GTF attributes look like:
      gene_id "ENSG..."; gene_version "1"; gene_name "BRCA1"; ...
    """
    d: Dict[str, str] = {}
    for part in attr_str.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        # split at first space
        if " " not in part:
            continue
        k, v = part.split(" ", 1)
        v = v.strip().strip('"')
        d[k] = v
    return d

def load_genes_json(path: str) -> Set[str]:
    """
    Accepts common shapes:
      - ["BRCA1", "BRCA2", ...]
      - {"genes": ["BRCA1", ...]}
      - {"BRCA1": {...}, "BRCA2": {...}}  (keys are symbols)
      - [{"gene_symbol":"BRCA1"}, ...]
    Returns a set of gene symbols.
    """
    with open(path) as fh:
        data = json.load(fh)

    genes: Set[str] = set()

    if isinstance(data, list):
        for item in data:
            if isinstance(item, str):
                genes.add(item.strip())
            elif isinstance(item, dict):
                for k in ("gene_symbol", "symbol", "gene", "name", "gene_name"):
                    if k in item and isinstance(item[k], str):
                        genes.add(item[k].strip())
                        break

    elif isinstance(data, dict):
        if "genes" in data and isinstance(data["genes"], list):
            for item in data["genes"]:
                if isinstance(item, str):
                    genes.add(item.strip())
                elif isinstance(item, dict):
                    for k in ("gene_symbol", "symbol", "gene", "name", "gene_name"):
                        if k in item and isinstance(item[k], str):
                            genes.add(item[k].strip())
                            break
        else:
            # assume keys are gene symbols
            for k in data.keys():
                if isinstance(k, str):
                    genes.add(k.strip())

    genes = {g for g in genes if g and not g.startswith("#")}
    return genes

def merge_intervals(intervals: List[Tuple[int,int]]) -> List[Tuple[int,int]]:
    """Merge 0-based half-open intervals."""
    if not intervals:
        return []
    intervals.sort(key=lambda x: (x[0], x[1]))
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        ps, pe = merged[-1]
        if s <= pe:  # overlap/adjacent
            merged[-1] = (ps, max(pe, e))
        else:
            merged.append((s, e))
    return merged

def main():
    ap = argparse.ArgumentParser(description="Build BED from Ensembl GTF + genes.json")
    ap.add_argument("--gtf", required=True, help="Ensembl GTF (.gtf or .gtf.gz)")
    ap.add_argument("--genes_json", required=True, help="genes.json containing gene symbols")
    ap.add_argument("--feature", choices=["gene", "exon"], default="exon",
                    help="Which GTF feature to use (default: exon)")
    ap.add_argument("--chr_prefix", action="store_true",
                    help="Force chr prefix (e.g., 17 -> chr17)")
    ap.add_argument("--out", required=True, help="Output BED file (BED4)")
    args = ap.parse_args()

    gene_set = load_genes_json(args.genes_json)
    if not gene_set:
        print("ERROR: gene set is empty after parsing genes.json", file=sys.stderr)
        sys.exit(2)

    # Collect intervals per (gene_symbol, chrom)
    # In most cases each gene is on one chrom; we still key by chrom to be safe.
    per_gene_chr: Dict[Tuple[str,str], List[Tuple[int,int]]] = {}

    with open_maybe_gz(args.gtf) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attrs = c
            if feature != args.feature:
                continue

            a = parse_gtf_attrs(attrs)
            gene_name = a.get("gene_name") or a.get("gene_symbol")
            if not gene_name or gene_name not in gene_set:
                continue

            # BED is 0-based start, 1-based end
            s0 = int(start) - 1
            e1 = int(end)

            if args.chr_prefix:
                chrom = ensure_chr(chrom)

            key = (gene_name, chrom)
            per_gene_chr.setdefault(key, []).append((s0, e1))

    if not per_gene_chr:
        print("WARNING: No matching genes found in GTF for provided genes.json", file=sys.stderr)

    # Write merged intervals as BED4: chrom start end gene
    # If feature=exon and you want per-exon (not merged), remove merge_intervals().
    with open(args.out, "w") as out:
        for (gene, chrom), ivs in sorted(per_gene_chr.items(), key=lambda x: (x[0][1], x[0][0])):
            for s0, e1 in merge_intervals(ivs):
                out.write(f"{chrom}\t{s0}\t{e1}\t{gene}\n")

if __name__ == "__main__":
    main()
