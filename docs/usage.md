# Usage Guide

## Pipeline Overview

Bionl_Lean_Call is a flexible pipeline that can run in three modes:
1. **Mode 1**: Full pipeline starting from FASTQ files (runs variant calling + post-processing)
2. **Mode 2**: Post-processing from variant calling output directory (automatic file discovery)
3. **Mode 3**: Post-processing with custom file paths (explicit VCF/BAM specification)

---

## Input Options

### Mode 1: Full Pipeline (FASTQ Input)

Run the complete variant calling workflow followed by VEP annotation and reporting.

#### Required Parameters
- `--input`: Path to CSV samplesheet with FASTQ file paths
- `--outdir`: Directory where results will be saved

#### Sample Sheet Format

Create a CSV file with the following columns:

```csv
patient,sex,status,sample,fastq_1,fastq_2,lane
Patient001,XX,0,Sample_A,data/sample_A_R1.fastq.gz,data/sample_A_R2.fastq.gz,1
Patient002,XY,0,Sample_B,data/sample_B_R1.fastq.gz,data/sample_B_R2.fastq.gz,1
```

**Column Descriptions:**
- `patient`: Patient identifier (can be the same for multiple samples)
- `sex`: Biological sex - use `XX` (female), `XY` (male), or `0` (unknown)
- `status`: Disease status - `0` (unaffected) or `1` (affected)
- `sample`: Unique sample identifier
- `fastq_1`: Full path to R1 FASTQ file
- `fastq_2`: Full path to R2 FASTQ file
- `lane`: Sequencing lane number

#### Command Example

```bash
nextflow run main.nf \
  --input samples.csv \
  --outdir results
```

---

### Mode 2: Post-Processing from Variant Calling Output Directory

Use this mode when you have an existing variant calling results directory. The pipeline automatically discovers VCF and BAM files following the standard variant calling output structure.

#### Required Parameters
- `--variant_calling_outdir`: Path to variant calling results directory
- `--outdir`: Directory where post-processing results will be saved

#### How It Works

The pipeline automatically searches for files in the variant calling output structure:
- **VCF files**: `{variant_calling_outdir}/variant_calling/*/*/*.vcf.gz`
- **BAM files**: `{variant_calling_outdir}/preprocessing/mapped/*/*.sorted.bam`
- **BAI files**: `{variant_calling_outdir}/preprocessing/mapped/*/*.sorted.bam.bai`

Sample names are extracted from the directory structure.

#### Command Example

```bash
nextflow run main.nf \
  --variant_calling_outdir /path/to/variant_calling/results \
  --outdir results
```

**Note:** This mode assumes standard variant calling output structure. If your files are in custom locations, use Mode 3 instead.

---

### Mode 3: Post-Processing with Custom File Paths

Use this mode for maximum control when files are not in standard variant calling structure, or when you want to process specific samples with explicit paths.

#### Required Parameters
- `--post_samplesheet`: Path to CSV file with explicit VCF/BAM/BAI paths
- `--outdir`: Directory where results will be saved

#### Post Sample Sheet Format

Create a CSV file with the following columns:

```csv
sample,vcf,bam,bai
HG003,/path/to/HG003.deepvariant.vcf.gz,/path/to/HG003.sorted.bam,/path/to/HG003.sorted.bam.bai
HG004,/path/to/HG004.deepvariant.vcf.gz,/path/to/HG004.sorted.bam,/path/to/HG004.sorted.bam.bai
```

**Column Descriptions:**
- `sample`: Sample identifier (will be used in output file names)
- `vcf`: Full path to variant call file (gzipped VCF)
- `bam`: Full path to aligned BAM file
- `bai`: Full path to BAM index file

**Important Notes:**
- All paths must be absolute
- VCF files should be gzipped (`.vcf.gz`)
- BAM files must have corresponding index files
- Files must exist and be readable

#### Command Example

```bash
nextflow run main.nf \
  --post_samplesheet data/post_samplesheet.csv \
  --outdir results
```

---

## Additional Parameters

### Common Options

- `--genome`: Reference genome (default: `GRCh38`)
- `--target_bed`: BED file with target regions for coverage analysis
- `--vep_cache`: Path to VEP cache directory
- `--vep_plugins`: Path to VEP plugins directory

### Advanced Options

- `--min_coverage`: Minimum coverage threshold for reporting (default: 20)
- `--min_vaf`: Minimum variant allele frequency (default: 0.2)

### Resource Configuration

- `--max_cpus`: Maximum CPUs to use
- `--max_memory`: Maximum memory to allocate
- `--max_time`: Maximum run time

---

## Output Structure

Results are organized by sample:

```
outdir/
├── SAMPLE1/
│   ├── vcf/
│   │   ├── SAMPLE1.filtered.vcf.gz
│   │   └── SAMPLE1.annotated.vcf.gz
│   ├── qc/
│   │   ├── coverage_summary.txt
│   │   ├── coverage_per_exon.tsv
│   │   ├── gaps.tsv
│   │   └── sample_summary.json
│   └── reports/
│       ├── SAMPLE1_acmg_report.xlsx
│       └── SAMPLE1_report.html
├── SAMPLE2/
│   └── ...
└── pipeline_report.html
```

### Output Files Description

**VCF Directory:**
- Filtered and annotated variant calls
- Includes functional annotations (VEP)
- ACMG SF gene variants

**QC Directory:**
- Coverage statistics at 20x, 30x, 50x, 100x thresholds
- Per-exon coverage metrics
- Coverage gaps in target regions
- Sample quality summary (JSON)

**Reports Directory:**
- Excel report with variant classifications
- HTML report with interactive visualizations
- ACMG secondary findings

---

## Complete Usage Examples

### Example 1: Run full pipeline with FASTQ files

```bash
nextflow run main.nf \
  --samplesheet samples.csv \
  --outdir results \
  --genome GRCh38 \
  --target_bed regions.bed
```

### Example 2: Post-process from variant calling output

```bash
nextflow run main.nf \
  --variant_calling_outdir /data/variant_calling_results \
  --outdir lean_call_results \
  --target_bed ACMG_genes.bed
```

### Example 3: Post-process with explicit file paths

```bash
nextflow run main.nf \
  --post_samplesheet post_samples.csv \
  --outdir results \
  --vep_cache /ref/vep_cache \
  --vep_plugins /ref/vep_plugins
```

### Example 4: Resume a failed run

```bash
nextflow run main.nf \
  --samplesheet samples.csv \
  --outdir results \
  -resume
```

---
# Validated Configuration (GIAB Benchmarking)

The workflow was benchmarked using Genome in a Bottle (GIAB) datasets
and achieved ≥95% precision within ACMG SF v3.3 regions using the
following validated configuration.

## Pipeline Components

- **Aligner:** BWA-MEM2
- **Variant Callers:** DeepVariant + GATK HaplotypeCaller
- **Variant Strategy:** Consensus integration (DV + HC)
- **Reference Genome:** GRCh38
- **Target Regions:** ACMG SF v3.3 BED

## Region & Variant Filtering

- Variants restricted to ACMG SF regions
- Only `PASS` variants retained
- GIAB high-confidence regions used during benchmarking

## Recommended Quality Thresholds

- **Minimum Depth (DP):** ≥ 20
- **Minimum Genotype Quality (GQ):** ≥ 20
- **Variant Status:** PASS only

## Benchmarking Notes

- Validation performed using GIAB HG003 and HG002 samples.
- Precision metrics were calculated after ACMG region filtering.
- Performance may vary if caller combinations or filtering thresholds are modified.

---

## Tips and Best Practices

### Choosing the Right Input Mode

- Use **FASTQ input** when:
  - Starting from raw sequencing data
  - You need complete control over alignment and variant calling
  - Running samples for the first time

- Use **post-processing mode** when:
  - You already have variant calling results
  - You want to re-annotate variants without re-running alignment
  - You need to update reports or change filtering criteria
  - You have BAM/VCF files from compatible pipelines

### File Path Guidelines

- Always use **absolute paths** in sample sheets for reliability
- Ensure VCF files are **indexed** (`.tbi` for gzipped VCFs)
- Verify BAM files have corresponding **`.bai` index files**
- Keep index files in the same directory as the data files

### Performance Optimization

- For large cohorts, consider running samples in batches
- Use `--max_cpus` and `--max_memory` to control resource usage
- Enable `-resume` to restart from the last completed step if interrupted

### Quality Control

- Review coverage statistics in the QC directory
- Check for coverage gaps in critical regions
- Verify sex determination matches expected values
- Review Ti/Tv ratios for data quality indicators

---

## Troubleshooting

### Common Issues

**Issue: Sample not found in variant calling output**
- Solution: Check that sample names match directory names in variant calling output
- Verify the variant calling pipeline completed successfully

**Issue: VCF or BAM file not found**
- Solution: Verify all paths in post_samplesheet.csv are absolute paths
- Ensure files exist and are accessible

**Issue: VEP annotation fails**
- Solution: Check that VEP cache and plugins are correctly specified
- Verify cache version matches genome assembly

**Issue: Out of memory errors**
- Solution: Increase `--max_memory` parameter
- Consider processing fewer samples at once

---

## Getting Help

If you encounter issues:

1. Check the Nextflow log files in the `work/` directory
2. Review `.command.log` files in failed task directories
3. Contact your bioinformatics support team
4. Email: support@example.com

For bug reports or feature requests, please include:
- Pipeline version
- Command used
- Error messages
- Relevant log files

