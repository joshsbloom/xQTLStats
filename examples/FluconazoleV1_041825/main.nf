#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.keyfile   = null
params.fastq_dir = null
params.genome    = null
params.vcf       = null
params.outdir    = null
params.workdir   = null
params.conda_env = 'genomics'

// ── Validate required parameters ────────────────────────────
if (!params.keyfile)   error "Missing required param: --keyfile"
if (!params.fastq_dir) error "Missing required param: --fastq_dir"
if (!params.genome)    error "Missing required param: --genome"
if (!params.vcf)       error "Missing required param: --vcf"
if (!params.outdir)    error "Missing required param: --outdir"

// ── Parse key file and resolve FASTQ pairs into the channel ──
// Row index (1-based) is used as the Illumina S-number in FASTQ filenames
Channel
    .fromPath(params.keyfile)
    .splitCsv(header: true)
    .toList()
    .flatMap { rows ->
        rows.withIndex().collect { row, idx ->
            def sid = row.'Sample ID'
            def s = idx + 1
            // BaseSpace puts FASTQs in a sample-named subdirectory
            def r1 = file("${params.fastq_dir}/${sid}*/${sid}_S${s}_L001_R1_001.fastq.gz")
            def r2 = file("${params.fastq_dir}/${sid}*/${sid}_S${s}_L001_R2_001.fastq.gz")
            // file() with glob returns a list; fall back to flat directory
            if (!r1) r1 = file("${params.fastq_dir}/${sid}_S${s}_L001_R1_001.fastq.gz")
            if (!r2) r2 = file("${params.fastq_dir}/${sid}_S${s}_L001_R2_001.fastq.gz")
            def f1 = r1 instanceof List ? r1[0] : r1
            def f2 = r2 instanceof List ? r2[0] : r2
            if (!f1.exists()) error "FASTQ R1 not found for sample ${sid} (S${s}): looked in ${params.fastq_dir}/${sid}*/ and ${params.fastq_dir}/"
            if (!f2.exists()) error "FASTQ R2 not found for sample ${sid} (S${s}): looked in ${params.fastq_dir}/${sid}*/ and ${params.fastq_dir}/"
            tuple(sid, f1, f2)
        }
    }
    .set { samples_ch }

process ALIGN_SORT_COUNT {
    tag "${sample_id}"
    cpus 2
    memory '4 GB'
    time '1h'

    publishDir "${params.outdir}/bam_files",   mode: 'copy', pattern: '*.bam*'
    publishDir "${params.outdir}/count_files",  mode: 'copy', pattern: '*.txt'

    input:
    tuple val(sample_id), path(r1), path(r2)
    path genome
    path genome_indices
    path vcf
    path vcf_index

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam
    tuple val(sample_id), path("${sample_id}.txt"), emit: counts

    script:
    """
    # Align and sort
    bwa mem -t ${task.cpus} \\
        -R '@RG\\tID:${sample_id}\\tSM:${sample_id}' \\
        ${genome} ${r1} ${r2} | \\
    sambamba view --nthreads=${task.cpus} --sam-input --format=bam --with-header /dev/stdin | \\
    sambamba sort --nthreads=${task.cpus} --tmpdir=\${TMPDIR:-/tmp} --out=${sample_id}.bam /dev/stdin

    # Count allele-specific reads
    gatk ASEReadCounter \\
        --min-mapping-quality 20 \\
        -I ${sample_id}.bam \\
        -V ${vcf} \\
        -R ${genome} \\
        -O ${sample_id}.txt
    """
}

workflow {
    // Collect reference files as value channels
    genome_file    = file(params.genome)
    // Include both sacCer3.fasta.* (bwa indices, .fai) and sacCer3.dict (GATK)
    genome_dir     = file(params.genome).parent
    genome_base    = file(params.genome).baseName  // "sacCer3" from "sacCer3.fasta"
    genome_indices = Channel.fromPath(["${params.genome}.*", "${genome_dir}/${genome_base}.dict"]).collect()
    vcf_file       = file(params.vcf)
    vcf_index      = file("${params.vcf}.*")

    ALIGN_SORT_COUNT(samples_ch, genome_file, genome_indices, vcf_file, vcf_index)
}
