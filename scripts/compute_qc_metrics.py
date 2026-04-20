#!/usr/bin/env python3
"""
Aggregate per-sample QC metrics into a machine-readable JSON for DB ingestion.

Reuses Sarek outputs wherever possible:
  --samtools-stats        samtools stats *.md.cram.stats   (from Sarek)
  --mosdepth-summary      mosdepth *.md.mosdepth.summary.txt (from Sarek)
  --mosdepth-global-dist  mosdepth *.md.mosdepth.global.dist.txt (from Sarek)
  --bcftools-stats        bcftools stats on consensus VCF  (run by QC module)
  --callers               optional caller composition      (run by QC module)

Outputs:
  <outdir>/<sample>.qc.json   structured QC report
  <outdir>/<sample>.qc.tsv    single-row summary for aggregation
"""

import argparse
import json
import re
from pathlib import Path

QC_VERSION = "1.0"

THRESHOLDS = {
    "WES": {
        "min_mapped_pct": 95.0,
        "max_duplicate_pct": 30.0,
        "min_mean_coverage": 80.0,
        "min_pct_20x": 95.0,
    },
    "WGS": {
        "min_mapped_pct": 95.0,
        "max_duplicate_pct": 20.0,
        "min_mean_coverage": 30.0,
        "min_pct_20x": 90.0,
    },
}
DEFAULT_THRESHOLDS = THRESHOLDS["WES"]


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_samtools_stats(path):
    """Extract alignment metrics from Sarek's samtools stats SN section.

    Format:  SN\\tkey:\\tvalue\\t# optional comment
    """
    raw = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith("SN\t"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 3:
                continue
            key = cols[1].rstrip(":").strip()
            val = cols[2].split("#")[0].strip()
            raw[key] = val

    total = int(raw.get("raw total sequences", 0))
    mapped = int(raw.get("reads mapped", 0))
    duplicated = int(raw.get("reads duplicated", 0))
    properly_paired = int(raw.get("reads properly paired", 0))
    properly_paired_pct = float(raw.get("percentage of properly paired reads (%)", 0))

    return {
        "total_reads": total,
        "mapped_reads": mapped,
        "mapped_pct": round(100.0 * mapped / total, 2) if total > 0 else 0.0,
        "duplicate_reads": duplicated,
        "duplicate_pct": round(100.0 * duplicated / total, 2) if total > 0 else 0.0,
        "properly_paired_reads": properly_paired,
        "properly_paired_pct": round(properly_paired_pct, 2),
    }


def parse_markdup_metrics(path):
    """Extract alignment metrics from Picard MarkDuplicates metrics.

    Fallback when samtools stats is unavailable.
    Properly paired info is not available from this source.
    """
    header = None
    data = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("#") or line.startswith("##") or not line:
                continue
            if header is None:
                if "UNPAIRED_READS_EXAMINED" in line:
                    header = line.split("\t")
                continue
            data = line.split("\t")
            break

    if header is None or data is None:
        return {
            "total_reads": 0, "mapped_reads": 0, "mapped_pct": 0.0,
            "duplicate_reads": 0, "duplicate_pct": 0.0,
            "properly_paired_reads": 0, "properly_paired_pct": 0.0,
        }

    def col(name):
        try:
            return data[header.index(name)]
        except (ValueError, IndexError):
            return "0"

    unpaired = int(col("UNPAIRED_READS_EXAMINED"))
    pairs = int(col("READ_PAIRS_EXAMINED"))
    unmapped = int(col("UNMAPPED_READS"))
    dup_unpaired = int(col("UNPAIRED_READ_DUPLICATES"))
    dup_pairs = int(col("READ_PAIR_DUPLICATES"))

    total = unpaired + (pairs * 2) + unmapped
    mapped = unpaired + (pairs * 2)
    duplicated = dup_unpaired + (dup_pairs * 2)

    try:
        dup_pct_raw = col("PERCENT_DUPLICATION")
        dup_pct = round(float(dup_pct_raw) * 100.0, 2) if dup_pct_raw else 0.0
    except ValueError:
        dup_pct = round(100.0 * duplicated / total, 2) if total > 0 else 0.0

    return {
        "total_reads": total,
        "mapped_reads": mapped,
        "mapped_pct": round(100.0 * mapped / total, 2) if total > 0 else 0.0,
        "duplicate_reads": duplicated,
        "duplicate_pct": dup_pct,
        "properly_paired_reads": 0,
        "properly_paired_pct": 0.0,
    }


def parse_mosdepth_summary(path):
    """Extract mean coverage from Sarek's mosdepth summary.

    Prefers ``total_region`` (present when --by is used).
    Falls back to ``total`` (genome-wide).
    """
    mean_cov = 0.0
    fallback_cov = 0.0
    with open(path) as fh:
        for line in fh:
            cols = line.strip().split("\t")
            if len(cols) < 4:
                continue
            if cols[0] == "total_region":
                mean_cov = float(cols[3])
            elif cols[0] == "total":
                fallback_cov = float(cols[3])
    cov = mean_cov if mean_cov > 0 else fallback_cov
    return round(cov, 2)


def parse_mosdepth_global_dist(path):
    """Derive coverage threshold percentages from Sarek's mosdepth global.dist.txt.

    Each line: region\\tdepth\\tfraction_at_or_above.
    We extract the ``total`` row at depths 10, 20, 30, 50.
    """
    dist = {}
    with open(path) as fh:
        for line in fh:
            cols = line.strip().split("\t")
            if len(cols) < 3 or cols[0] != "total":
                continue
            try:
                depth = int(cols[1])
                frac = float(cols[2])
            except ValueError:
                continue
            dist[depth] = frac

    return {
        t: round(dist.get(t, 0.0) * 100.0, 2)
        for t in (10, 20, 30, 50)
    }


def parse_bcftools_stats(path):
    """Extract variant counts, Ti/Tv, and het/hom from bcftools stats."""
    metrics = {
        "pass_variants": 0, "snvs": 0, "indels": 0,
        "titv": 0.0, "het_hom_ratio": 0.0,
    }
    n_het, n_hom = 0, 0

    with open(path) as fh:
        for line in fh:
            if line.startswith("SN\t"):
                cols = line.strip().split("\t")
                if len(cols) < 4:
                    continue
                key = cols[2].rstrip(":").strip()
                val = cols[3].strip()
                if key == "number of records":
                    metrics["pass_variants"] = int(val)
                elif key == "number of SNPs":
                    metrics["snvs"] = int(val)
                elif key == "number of indels":
                    metrics["indels"] = int(val)
            elif line.startswith("TSTV\t"):
                cols = line.strip().split("\t")
                if len(cols) >= 5:
                    try:
                        metrics["titv"] = round(float(cols[4]), 4)
                    except ValueError:
                        pass
            elif line.startswith("PSC\t"):
                cols = line.strip().split("\t")
                if len(cols) >= 6:
                    try:
                        n_hom = int(cols[4])  # nNonRefHom
                        n_het = int(cols[5])  # nHet
                    except (ValueError, IndexError):
                        pass

    metrics["het_hom_ratio"] = round(n_het / n_hom, 4) if n_hom > 0 else 0.0
    return metrics


def parse_caller_composition(path):
    """Count caller annotations from an optional CALLERS query file."""
    counts = {}
    if path is None or not Path(path).exists() or Path(path).stat().st_size == 0:
        return counts
    with open(path) as fh:
        for line in fh:
            val = line.strip()
            if val and val != ".":
                counts[val] = counts.get(val, 0) + 1
    return counts


# ---------------------------------------------------------------------------
# QC evaluation
# ---------------------------------------------------------------------------

def evaluate_qc(alignment, coverage, thresholds):
    """Compare metrics against thresholds; return status, recommendation, flags."""
    flags = []

    if alignment.get("mapped_pct", 100) < thresholds["min_mapped_pct"]:
        flags.append(f"LOW_MAPPED_PCT:{alignment['mapped_pct']}")

    if alignment.get("duplicate_pct", 0) > thresholds["max_duplicate_pct"]:
        flags.append(f"HIGH_DUPLICATE_PCT:{alignment['duplicate_pct']}")

    if coverage["mean_coverage"] < thresholds["min_mean_coverage"]:
        flags.append(f"LOW_MEAN_COVERAGE:{coverage['mean_coverage']}")

    if coverage["pct_20x"] < thresholds["min_pct_20x"]:
        flags.append(f"LOW_20X_COVERAGE:{coverage['pct_20x']}")

    qc_status = "FAIL" if flags else "PASS"
    recommendation = "HOLD" if flags else "READY"
    return qc_status, recommendation, flags


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samtools-stats", default=None,
                    help="samtools stats output (from Sarek, preferred)")
    ap.add_argument("--markdup-metrics", default=None,
                    help="Picard MarkDuplicates metrics (from Sarek, fallback)")
    ap.add_argument("--mosdepth-summary", required=True,
                    help="mosdepth summary (from Sarek)")
    ap.add_argument("--mosdepth-global-dist", required=True,
                    help="mosdepth global.dist.txt (from Sarek)")
    ap.add_argument("--bcftools-stats", required=True,
                    help="bcftools stats on consensus VCF")
    ap.add_argument("--callers", default=None,
                    help="optional caller composition file")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--assay", default="NA")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    assay_upper = args.assay.upper()
    thresholds = THRESHOLDS.get(assay_upper, DEFAULT_THRESHOLDS)

    if args.samtools_stats:
        alignment = parse_samtools_stats(args.samtools_stats)
    elif args.markdup_metrics:
        alignment = parse_markdup_metrics(args.markdup_metrics)
    else:
        raise SystemExit("ERROR: provide --samtools-stats or --markdup-metrics")
    mean_cov = parse_mosdepth_summary(args.mosdepth_summary)
    cov_pcts = parse_mosdepth_global_dist(args.mosdepth_global_dist)

    coverage = {
        "scope": "genome_wide",
        "mean_coverage": mean_cov,
        "pct_10x": cov_pcts.get(10, 0.0),
        "pct_20x": cov_pcts.get(20, 0.0),
        "pct_30x": cov_pcts.get(30, 0.0),
        "pct_50x": cov_pcts.get(50, 0.0),
    }

    variant_summary = parse_bcftools_stats(args.bcftools_stats)
    caller_comp = parse_caller_composition(args.callers)

    qc_status, recommendation, flags = evaluate_qc(alignment, coverage, thresholds)

    result = {
        "sample": args.sample,
        "assay": args.assay,
        "qc_status": qc_status,
        "qc_version": QC_VERSION,
        "alignment": {
            "total_reads": alignment["total_reads"],
            "mapped_reads": alignment["mapped_reads"],
            "mapped_pct": alignment["mapped_pct"],
            "duplicate_reads": alignment["duplicate_reads"],
            "duplicate_pct": alignment["duplicate_pct"],
            "properly_paired_reads": alignment["properly_paired_reads"],
            "properly_paired_pct": alignment["properly_paired_pct"],
        },
        "coverage": coverage,
        "variant_summary": {
            "pass_variants": variant_summary["pass_variants"],
            "snvs": variant_summary["snvs"],
            "indels": variant_summary["indels"],
            "titv": variant_summary["titv"],
            "het_hom_ratio": variant_summary["het_hom_ratio"],
        },
        "thresholds": thresholds,
        "flags": flags,
        "db_ingestion_recommendation": recommendation,
    }

    if caller_comp:
        result["variant_summary"]["caller_composition"] = caller_comp

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # --- JSON ---
    json_path = outdir / f"{args.sample}.qc.json"
    with open(json_path, "w") as fh:
        json.dump(result, fh, indent=2)

    # --- single-row TSV ---
    tsv_path = outdir / f"{args.sample}.qc.tsv"
    cols = [
        "sample", "assay", "qc_status", "db_ingestion_recommendation",
        "total_reads", "mapped_pct", "duplicate_pct",
        "mean_coverage", "pct_20x",
        "pass_variants", "snvs", "indels", "titv", "het_hom_ratio",
    ]
    vals = [
        args.sample, args.assay, qc_status, recommendation,
        str(alignment["total_reads"]),
        str(alignment["mapped_pct"]),
        str(alignment["duplicate_pct"]),
        str(coverage["mean_coverage"]),
        str(coverage["pct_20x"]),
        str(variant_summary["pass_variants"]),
        str(variant_summary["snvs"]),
        str(variant_summary["indels"]),
        str(variant_summary["titv"]),
        str(variant_summary["het_hom_ratio"]),
    ]
    with open(tsv_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        fh.write("\t".join(vals) + "\n")

    print(f"QC export: {json_path} — {qc_status} / {recommendation}")


if __name__ == "__main__":
    main()
