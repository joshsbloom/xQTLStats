#!/bin/bash
set -euo pipefail

# ── Usage ──────────────────────────────────────────────────────────
# Local:    ./run_pipeline.sh local
# Hoffman2: ./run_pipeline.sh hoffman2
#
# Override conda env name: --conda_env myenv
# Override work dir:       WORKDIR=/path/to/work ./run_pipeline.sh local
# Create env from scratch: mamba env create -f environment.yml
# ───────────────────────────────────────────────────────────────────

MODE=${1:-local}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE="${SCRIPT_DIR}/main.nf"

# ── Shared ─────────────────────────────────────────────────────────
KEYFILE="${SCRIPT_DIR}/flatkey.csv"

if [ "$MODE" == "local" ]; then
    # ── Local paths ────────────────────────────────────────────────
    FASTQ_DIR="/data2/xQTL/FluconazoleV1_041825/fastq"
    OUTDIR="/data2/xQTL/FluconazoleV1_041825"
    WORKDIR="${WORKDIR:-${OUTDIR}/work}"
    GENOME="/media/hoffman2/PUBLIC_SHARED/yeast/reference/1011/sacCer3.fasta"
    VCF="/media/hoffman2/PUBLIC_SHARED/yeast/reference/3000/VCFs_SNPs/1011_allhet_SNPs.vcf.gz"

    nextflow run "${PIPELINE}" \
        -profile local \
        --keyfile "${KEYFILE}" \
        --fastq_dir "${FASTQ_DIR}" \
        --outdir "${OUTDIR}" \
        --workdir "${WORKDIR}" \
        --genome "${GENOME}" \
        --vcf "${VCF}" \
        -resume

elif [ "$MODE" == "hoffman2" ]; then
    # ── Hoffman2 paths ─────────────────────────────────────────────
    # Change username to your hoffman2 username
    USERNAME="${USER}"
    PROJECT_ROOT="/u/project/kruglyak/${USERNAME}/projects/FluconazoleV1_041825"
    FASTQ_DIR="${PROJECT_ROOT}/fastq"
    OUTDIR="${PROJECT_ROOT}"
    WORKDIR="${WORKDIR:-/u/home/j/jsbloom/scratch/xQTL_work/FluconazoleV1_041825}"
    GENOME="/u/project/kruglyak/PUBLIC_SHARED/yeast/reference/1011/sacCer3.fasta"
    VCF="/u/project/kruglyak/PUBLIC_SHARED/yeast/reference/3000/VCFs_SNPs/1011_allhet_SNPs.vcf.gz"

    # Create logs dir before submitting
    mkdir -p "${OUTDIR}/logs"

    # Submit the nextflow master process as an SGE job.
    # Nextflow runs on this node and submits per-sample tasks via SGE,
    # so all 96 samples can run in parallel across the cluster.
    # 4G/1 core is enough for the master; individual tasks request their own resources.
    qsub -cwd \
         -o "${OUTDIR}/logs" \
         -e "${OUTDIR}/logs" \
         -l h_data=4G,h_rt=48:00:00,highp \
         -pe shared 2 \
         -N xQTL_nf \
         <<QSUB
#!/bin/bash
#\$ -V
#\$ -cwd

nextflow run ${PIPELINE} \\
    -profile hoffman2 \\
    --keyfile ${KEYFILE} \\
    --fastq_dir ${FASTQ_DIR} \\
    --outdir ${OUTDIR} \\
    --workdir ${WORKDIR} \\
    --genome ${GENOME} \\
    --vcf ${VCF} \\
    -resume
QSUB

    echo "Submitted nextflow master job to SGE. Monitor with: qstat -u ${USERNAME}"
    echo "Logs will be in: ${OUTDIR}/logs/"
    echo "Per-sample tasks are submitted as separate SGE jobs (up to 96 concurrent)."

else
    echo "Usage: $0 {local|hoffman2}"
    exit 1
fi
