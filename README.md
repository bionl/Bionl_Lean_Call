# Bionl_Lean_Call

Germline variant calling and reporting for ACMG SF genes, including QC and functional annotation—designed for research, not diagnostic use.

---

## Inputs

The pipeline can run in three modes:

### Mode 1: Full Pipeline (from FASTQ files)

Run the complete variant calling pipeline followed by post-processing.

**Required parameters:**
- `--input`: Samplesheet CSV
- `--outdir`: Output directory

**Example samplesheet:**

```csv
patient,sex,status,sample,fastq_1,fastq_2,lane
Patient001,XX,0,Sample_A,data/sample_A_R1.fastq.gz,data/sample_A_R2.fastq.gz,1
Patient002,XY,0,Sample_B,data/sample_B_R1.fastq.gz,data/sample_B_R2.fastq.gz,1
Patient003,XX,1,Sample_C,data/sample_C_R1.fastq.gz,data/sample_C_R2.fastq.gz,1
```

The CSV should include:
- **patient**: Patient identifier
- **sex**: Sex (XX for female, XY for male, or 0 if unknown)
- **status**: Affected status (0 for unaffected, 1 for affected)
- **sample**: Sample identifier (must be unique)
- **fastq_1**: Path to forward reads file
- **fastq_2**: Path to reverse reads file
- **lane**: Sequencing lane number

### Mode 2: Post-Processing from Variant Calling Output Directory

Use existing variant calling results (automatic file discovery).

**Required parameters:**
- `--variant_calling_outdir`: Path to variant calling results directory
- `--outdir`: Output directory

The pipeline automatically finds VCF and BAM files in the variant calling output structure:
- VCF files: `{variant_calling_outdir}/variant_calling/*/*/*.vcf.gz`
- BAM files: `{variant_calling_outdir}/preprocessing/mapped/*/*.sorted.bam`

### Mode 3: Post-Processing with Custom File Paths

Use a CSV to specify exact paths to VCF and BAM files.

**Required parameters:**
- `--post_samplesheet`: CSV with sample, VCF, BAM, and BAI paths
- `--outdir`: Output directory

**Example post_samplesheet.csv:**

```csv
sample,vcf,bam,bai
HG003,/path/to/HG003.deepvariant.vcf.gz,/path/to/HG003.sorted.bam,/path/to/HG003.sorted.bam.bai
HG004,/path/to/HG004.deepvariant.vcf.gz,/path/to/HG004.sorted.bam,/path/to/HG004.sorted.bam.bai
```

The CSV should include:
- **sample**: Sample identifier
- **vcf**: Full path to variant call file
- **bam**: Full path to aligned BAM file
- **bai**: Full path to BAM index file

---

## Outputs

Results are organized by sample in the output directory:

```
outdir/
└── {SAMPLE}/
    ├── vcf/          # Variant files (filtered and annotated)
    ├── qc/           # Quality metrics and coverage statistics
    └── reports/      # Excel and HTML reports
```

**What you'll find:**

- **VCF files**: Identified genetic variants with annotations
- **QC metrics**: 
  - Coverage analysis (20x, 30x, 50x, 100x thresholds)
  - Coverage gaps identification
  - Alignment statistics (mapped reads, duplicates)
  - Read balance metrics (R1/R2 ratios, strand balance)
  - Sex determination from sequencing data
  - Variant statistics (transition/transversion ratios)
- **Reports**: Easy-to-read Excel summaries with variant classifications and clinical annotations

---

## How to Run

### Using the Platform UI

1. Navigate to the pipeline in your platform interface
2. Choose your input mode:
   - **Mode 1**: Upload `--input` samplesheet (FASTQ files)
   - **Mode 2**: Specify `--variant_calling_outdir` (existing variant calling results)
   - **Mode 3**: Upload `--post_samplesheet` (custom VCF/BAM paths)
3. Specify your **output directory** path
4. Click **Run**

### Command Line Examples

**Mode 1 - Full pipeline from FASTQ files:**
```bash
nextflow run main.nf \
  --input samples.csv \
  --outdir results
```

**Mode 2 - Post-process existing variant calling output:**
```bash
nextflow run main.nf \
  --variant_calling_outdir /path/to/variant_calling/results \
  --outdir results
```

**Mode 3 - Post-process with custom file paths:**
```bash
nextflow run main.nf \
  --post_samplesheet data/post_samplesheet.csv \
  --outdir results
```

**Post-processing from Sarek output directory:**
```bash
nextflow run main.nf \
  --sarek_outdir /path/to/sarek/results \
  --outdir results
```

**Post-processing with explicit file paths:**
```bash
nextflow run main.nf \
  --post_samplesheet data/post_samplesheet.csv \
  --outdir results
```

That's it! The pipeline handles everything else automatically.

---
# Validation & Benchmarking

This workflow was benchmarked using Genome in a Bottle (GIAB) reference samples
to validate performance within clinically relevant ACMG Secondary Findings regions.

## Benchmarking Overview

Datasets:
- HG003 (GIAB)
- HG002 (GIAB)

Reference:
- GRCh38

Target Regions:
- ACMG SF v3.3 BED

Pipeline Configuration:
- Aligner: BWA-MEM
- Variant Callers: DeepVariant + GATK HaplotypeCaller
- Strategy: Consensus variant integration (DV+HC)
- Filtering:
  - ACMG region restriction
  - Depth and genotype quality thresholds
  - PASS variants only

## Performance Summary (ACMG Regions)

| Metric | Result |
|---|---|
| Precision | 94.8% |
| Recall | 99.5% |
| F1 Score | 97.1% |

Benchmarking was performed against GIAB high-confidence truth sets using hap.py.

## Reproducibility

To reproduce validated results, use the configuration shown in `USAGE.md`.

## Full Validation Documentation

See detailed validation report: docs/validation/GIAB_validation_report.md

## Support

For questions or assistance:
- Contact your bioinformatics team
- Email: khatib@bionl.ai

**Pipeline version**: 1.0.0
