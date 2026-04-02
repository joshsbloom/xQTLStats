# Fluconazole xQTL Example

Bulk-segregant QTL mapping of fluconazole resistance in yeast biparental crosses. This example covers the full workflow: FASTQ alignment, allele counting, and statistical analysis.

## Prerequisites

- [Nextflow](https://www.nextflow.io/) (>= 24.x)
- [Conda](https://docs.conda.io/) or [Mamba](https://mamba.readthedocs.io/)
- FASTQ files downloaded from BaseSpace (see below)

## 1. Set up the conda environment

### From the environment file (recommended)

```bash
mamba env create -f environment.yml
conda activate genomics
```

### On Hoffman2

```bash
# Request an interactive node first
qrsh -l highp,h_rt=48:00:00,h_data=8G -pe shared 8

# Load mamba and create the environment
module load mamba
mamba env create -f environment.yml
conda activate genomics
```

### Install the xQTLStats R package (inside the environment)

```r
devtools::install_github("joshsbloom/xQTLStats", ref = "main")
```

The pipeline assumes the conda environment is named `genomics`. If yours has a different name, either pass it on the command line with `--conda_env myenv` or change the `conda_env` parameter in `nextflow.config`.

## 2. Download FASTQ files

FASTQs are downloaded from Illumina BaseSpace. Authenticate and download to your project's fastq directory:

```bash
bs auth
bs list projects

# Replace the project ID and output path for your setup
bs download project -i 451502321 \
    -o /path/to/your/project/fastq/ \
    --extension fastq.gz
```

BaseSpace creates a subdirectory per sample (e.g., `fastq/Sample001/Sample001_S1_L001_R1_001.fastq.gz`). The pipeline handles this naming convention automatically.

## 3. Run the alignment + counting pipeline

The Nextflow pipeline (`main.nf`) runs per sample in parallel: `bwa mem` alignment, `sambamba` sort, and `gatk ASEReadCounter`.

### Local

```bash
conda activate genomics
./run_pipeline.sh local
```

### Hoffman2

```bash
./run_pipeline.sh hoffman2
```

This submits a Nextflow master job via `qsub`, which then submits individual per-sample SGE jobs (up to 96 concurrent). Monitor with `qstat -u $USER`.

### Customization

Edit `run_pipeline.sh` to adjust paths for your setup:
- `FASTQ_DIR` — location of downloaded FASTQ files
- `OUTDIR` — where BAMs and count files are written
- `GENOME` — path to indexed sacCer3 reference FASTA
- `VCF` — path to the VCF used for allele counting

Edit `nextflow.config` to adjust resource limits:
- `executor.cpus` / `executor.memory` — local parallelism budget
- `executor.queueSize` — max concurrent SGE jobs on Hoffman2

### Pipeline management

```bash
# Resume after interruption (picks up where it left off)
./run_pipeline.sh local    # -resume is on by default

# Clean cached work files
nextflow clean -f

# Full reset (removes all cached results and outputs)
rm -rf work/ .nextflow .nextflow.log*
```

### Output

- `{outdir}/bam_files/` — sorted BAM files per sample
- `{outdir}/count_files/` — GATK ASEReadCounter output per sample
- `{outdir}/pipeline_info/` — Nextflow timeline and resource reports

## 4. Statistical analysis

After the pipeline completes, run the R analysis script:

```r
source("xQTL_analysis.R")
```

This script:
1. Loads the genotype reference (`Z_1011_filtered.qs`) and sample key file
2. Builds phased allele count tables per sample (parallelized with `pbmclapply`)
3. Computes smoothed allele frequency differences per chromosome
4. Calculates contrast statistics between drug and control conditions
5. Generates genome-wide plots of allele frequency and significance

Adjust `ncores` at the top of the script to match your machine.

## File overview

| File | Description |
|---|---|
| `environment.yml` | Conda environment specification |
| `flatkey.csv` | Sample key: sample IDs, parents, conditions |
| `main.nf` | Nextflow pipeline (align + sort + count) |
| `nextflow.config` | Nextflow config with local and hoffman2 profiles |
| `run_pipeline.sh` | Driver script for local and hoffman2 execution |
| `xQTL_analysis.R` | R analysis: phasing, AFD smoothing, contrasts |
| `orig/` | Original serial scripts (setup, alignment, analysis) for reference |
