#!/usr/bin/env python3
import argparse, sys


def ensure_chr(chrom: str) -> str:
    return chrom if chrom.startswith("chr") else f"chr{chrom}"

def load_bed_gene_map(bed_path: str):
    # map (chrom, pos1based) -> gene
    m = {}
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            chrom, start0, end, gene = line.strip().split("\t")[:4]
            start0 = int(start0)
            pos1 = start0 + 1
            m[(chrom, pos1)] = gene
    return m

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--assay", required=True, help="WES or WGS")
    ap.add_argument("--positions_bed", required=True, help="Same BED used for samtools depth")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    gene_map = load_bed_gene_map(args.positions_bed)

    with open(args.out, "w") as out:
        out.write("sample_id\tassay_type\tgene_symbol\tchrom\tposition\tdp\tge10\tge20\tge30\n")
        for line in sys.stdin:
            if not line.strip():
                continue
            chrom, pos, dp = line.strip().split("\t")[:3]
            pos_i = int(pos)
            dp_i = int(float(dp))
            gene = gene_map.get((chrom, pos_i), "NA")
            ge10 = 1 if dp_i >= 10 else 0
            ge20 = 1 if dp_i >= 20 else 0
            ge30 = 1 if dp_i >= 30 else 0
            chrom = ensure_chr(chrom)
            out.write(f"{args.sample}\t{args.assay}\t{gene}\t{chrom}\t{pos_i}\t{dp_i}\t{ge10}\t{ge20}\t{ge30}\n")

if __name__ == "__main__":
    main()
