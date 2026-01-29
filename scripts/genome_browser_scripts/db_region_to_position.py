#!/usr/bin/env python3
import argparse

def ensure_chr(chrom: str) -> str:
    return chrom if chrom.startswith("chr") else f"chr{chrom}"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regions_bed", required=True, help="BED: chrom start end geneSymbol")
    ap.add_argument("--step", type=int, default=200, help="Sampling step bp")
    ap.add_argument("--out", required=True, help="Output BED: chrom start end geneSymbol")
    args = ap.parse_args()

    with open(args.out, "w") as out:
        with open(args.regions_bed) as fh:
            for line in fh:
                if not line.strip():
                    continue
                chrom, start0, end, gene = line.strip().split("\t")[:4]
                start0 = int(start0); end = int(end)

                chrom = ensure_chr(chrom)

                positions = list(range(start0, end, args.step))
                if not positions or positions[-1] != end - 1:
                    positions.append(end - 1)

                for p0 in positions:
                    out.write(f"{chrom}\t{p0}\t{p0+1}\t{gene}\n")

if __name__ == "__main__":
    main()
