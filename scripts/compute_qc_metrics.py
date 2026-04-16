#!/usr/bin/env python3
"""
Aggregate per-sample QC metrics into a machine-readable JSON for DB ingestion.

Inputs (standard bioinformatics tool outputs):
  --flagstat              samtools flagstat output
  --mosdepth-summary      mosdepth *.mosdepth.summary.txt
  --mosdepth-thresholds   mosdepth *.thresholds.bed.gz
  --bcftools-stats        bcftools stats output
  --callers               optional caller composition (one CALLERS value per line)
  --sample                sample ID
  --assay                 assay type (WES, WGS, etc.)
  --outdir                output directory

Outputs:
  <outdir>/<sample>.qc.json   structured QC report
  <outdir>/<sample>.qc.tsv    single-row summary for aggregation
"""

import argparse
import gzip
import json
import re
import sys
from pathlib import Path

QC_VERSION = "1.0"

THRESHOLDS = {
    "WES": {
        "min_mapped_pct": 95.0,
        "max_duplicate_pct": 30.0,
        "min_target_mean_coverage": 80.0,
        "min_targets_20x_pct": 95.0,
    },
    "WGS": {
        "min_mapped_pct": 95.0,
        "max_duplicate_pct": 20.0,
        "min_target_mean_coverage": 30.0,
        "min_targets_20x_pct": 90.0,
    },
}
DEFAULT_THRESHOLDS = THRESHOLDS["WES"]


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_flagstat(path):
    """Extract alignment metrics from samtools flagstat output."""
    metrics = {
        "total_reads": 0, "mapped_reads": 0, "mapped_pct": 0.0,
        "duplicate_reads": 0, "duplicate_pct": 0.0,
        "properly_paired_reads": 0, "properly_paired_pct": 0.0,
    }
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            m = re.match(r"(\d+)\s+\+\s+\d+\s+(.*)", line)
            if not m:
                continue
            count = int(m.group(1))
            desc = m.group(2)

            if "in total" in desc:
                metrics["total_reads"] = count
            elif desc.startswith("duplicates"):
                metrics["duplicate_reads"] = count
            elif desc.startswith("mapped") and "primary" not in desc:
                metrics["mapped_reads"] = count
                pct = re.search(r"\(([0-9.]+)%", desc)
                metrics["mapped_pct"] = float(pct.group(1)) if pct else 0.0
            elif "properly paired" in desc:
                metrics["properly_paired_reads"] = count
                pct = re.search(r"\(([0-9.]+)%", desc)
                metrics["properly_paired_pct"] = float(pct.group(1)) if pct else 0.0

    total = metrics["total_reads"]
    dup = metrics["duplicate_reads"]
    metrics["duplicate_pct"] = round(100.0 * dup / total, 2) if total > 0 else 0.0
    return metrics


def parse_mosdepth_summary(path):
    """Extract mean coverage from mosdepth summary.

    Prefers ``total_region`` (present when --by is used, i.e. WES with BED).
    Falls back to ``total`` (genome-wide, i.e. WGS without BED).
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
    return {"target_mean_coverage": round(cov, 2)}


def parse_mosdepth_thresholds(path):
    """Compute % of target bases at >=Nx from mosdepth thresholds BED."""
    open_fn = gzip.open if str(path).endswith(".gz") else open

    threshold_cols = []
    total_bases = {}
    total_region_size = 0

    with open_fn(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                header = line.lstrip("#").split("\t")
                for i, label in enumerate(header):
                    if re.match(r"^\d+[Xx]$", label):
                        threshold_cols.append((i, label.lower()))
                continue
            cols = line.split("\t")
            if not threshold_cols or len(cols) < 3:
                continue
            start, end = int(cols[1]), int(cols[2])
            total_region_size += end - start
            for col_idx, label in threshold_cols:
                if col_idx < len(cols):
                    total_bases[label] = total_bases.get(label, 0) + int(cols[col_idx])

    result = {}
    for label, bases in total_bases.items():
        pct = round(100.0 * bases / total_region_size, 2) if total_region_size > 0 else 0.0
        result[f"targets_{label}_pct"] = pct
    return result


def parse_mosdepth_global_dist(path):
    """Derive coverage threshold percentages from mosdepth global.dist.txt.

    Each line is: region \\t depth \\t fraction_at_or_above.
    We look for the ``total`` region at depths 10, 20, 30, 50.
    Used as a fallback for WGS when no --by BED thresholds are available.
    """
    dist = {}
    if path is None or not Path(path).exists():
        return {}
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

    result = {}
    for threshold in (10, 20, 30, 50):
        frac = dist.get(threshold, 0.0)
        result[f"targets_{threshold}x_pct"] = round(frac * 100.0, 2)
    return result


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

    if coverage.get("target_mean_coverage", 0) < thresholds["min_target_mean_coverage"]:
        flags.append(f"LOW_MEAN_COVERAGE:{coverage['target_mean_coverage']}")

    pct_20x = coverage.get("targets_20x_pct", 0)
    if pct_20x < thresholds["min_targets_20x_pct"]:
        flags.append(f"LOW_20X_COVERAGE:{pct_20x}")

    qc_status = "FAIL" if flags else "PASS"
    recommendation = "HOLD" if flags else "READY"
    return qc_status, recommendation, flags


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--flagstat", required=True)
    ap.add_argument("--mosdepth-summary", required=True)
    ap.add_argument("--mosdepth-thresholds", required=True)
    ap.add_argument("--mosdepth-global-dist", default=None,
                    help="mosdepth global.dist.txt (fallback for WGS without BED)")
    ap.add_argument("--bcftools-stats", required=True)
    ap.add_argument("--callers", default=None)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--assay", default="NA")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    assay_upper = args.assay.upper()
    thresholds = THRESHOLDS.get(assay_upper, DEFAULT_THRESHOLDS)

    alignment = parse_flagstat(args.flagstat)
    cov_summary = parse_mosdepth_summary(args.mosdepth_summary)
    cov_thresholds = parse_mosdepth_thresholds(args.mosdepth_thresholds)
    if cov_thresholds:
        coverage_scope = "on_target"
    elif args.mosdepth_global_dist:
        cov_thresholds = parse_mosdepth_global_dist(args.mosdepth_global_dist)
        coverage_scope = "genome_wide"
    else:
        coverage_scope = "genome_wide"
    coverage = {**cov_summary, **cov_thresholds}
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
        "coverage": {
            "scope": coverage_scope,
            "mean_coverage": coverage.get("target_mean_coverage", 0.0),
            "pct_10x": coverage.get("targets_10x_pct", 0.0),
            "pct_20x": coverage.get("targets_20x_pct", 0.0),
            "pct_30x": coverage.get("targets_30x_pct", 0.0),
            "pct_50x": coverage.get("targets_50x_pct", 0.0),
        },
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

    # --- single-row TSV (easy to aggregate later) ---
    tsv_path = outdir / f"{args.sample}.qc.tsv"
    cols = [
        "sample", "assay", "qc_status", "db_ingestion_recommendation",
        "total_reads", "mapped_pct", "duplicate_pct",
        "target_mean_coverage", "targets_20x_pct",
        "pass_variants", "snvs", "indels", "titv", "het_hom_ratio",
    ]
    vals = [
        args.sample, args.assay, qc_status, recommendation,
        str(alignment["total_reads"]),
        str(alignment["mapped_pct"]),
        str(alignment["duplicate_pct"]),
        str(coverage.get("target_mean_coverage", 0.0)),
        str(coverage.get("targets_20x_pct", 0.0)),
        str(variant_summary["pass_variants"]),
        str(variant_summary["snvs"]),
        str(variant_summary["indels"]),
        str(variant_summary["titv"]),
        str(variant_summary["het_hom_ratio"]),
    ]
    with open(tsv_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        fh.write("\t".join(vals) + "\n")

    print(f"QC export complete: {json_path} — status={qc_status}, recommendation={recommendation}")


if __name__ == "__main__":
    main()
