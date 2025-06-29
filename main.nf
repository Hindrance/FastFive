#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Parameters
params.fast5 = null
params.outdir = "results"
params.help = false

// Help message
def helpMessage() {
    log.info"""
    Usage:
    nextflow run FastFive.nf --fast5 <path_to_fast5> --reference <path_to_reference_genome>

    Required arguments:
    --fast5         Path to input FAST5 file

    Optional arguments:
    --outdir        Output directory (default: results)
    --help          Show this help message
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// Validate required parameters
if (!params.fast5 ) {
    log.error "Both --fast5 and --reference arguments are required"
    helpMessage()
    exit 1
}

// Input channels
ch_fast5 = Channel.fromPath(params.fast5, checkIfExists: true)
ch_reference = Channel.fromPath(params.reference, checkIfExists: true)

workflow {
    // Basecalling with Guppy
    BASECALLING(ch_fast5)
    
    // Quality control with NanoPlot
    QUALITY_CONTROL(BASECALLING.out.fastq)
    
    // Alignment with minimap2
    ALIGNMENT(BASECALLING.out.fastq, ch_reference)
    
    // quick SAMtools indexing
    INDEXING(ALIGNMENT.out.sam)
}

// Processes
process BASECALLING {
    container 'genomicpariscentre/guppy-gpu:6.5.7'
    publishDir "${params.outdir}/basecalling", mode: 'copy'
    
    input:
    path fast5
    
    output:
    path "pass/*.fastq*", emit: pass_fastq
    path "fail/*.fastq*", emit: fail_fastq, optional: true
    path "sequencing_summary.txt", emit: txt
    path "guppy_basecaller_log*", emit: log
    path "*.fastq", emit: fastq
    
    script:
    """
    guppy_basecaller \\
        --input_path ${fast5} \\
        --save_path . \\
        --flowcell FLO-MIN106 \\
        --kit SQK-LSK109 \\
        --num_callers 1
    
    # Concatenate all fastq files
    cat pass/*.fastq > combined_reads.fastq
    """
}


process QUALITY_CONTROL {
    container 'quay.io/biocontainers/nanoplot:1.41.0--pyhdfd78af_0'
    publishDir "${params.outdir}/qc", mode: 'copy'
    
    input:
    path fastq
    
    output:
    path "NanoPlot-report.html"
    path "NanoStats.txt"
    path "*.png"
    
    script:
    """
    NanoPlot \\
        --fastq ${fastq} \\
        --plots kde hex dot \\
        --format png \\
        --N50 \\
        --dpi 300 \\
        --store \\
        --raw \\
        --tsv_stats \\
        --info_in_report \\
        --prefix "" \\
        --threads ${task.cpus}
    """
}


process ALIGNMENT {
    container 'quay.io/biocontainers/minimap2:2.24--h7132678_1'
    publishDir "${params.outdir}/alignment", mode: 'copy'
    
    input:
    path fastq
    path reference
    
    output:
    path "aligned.sam", emit: sam
    
    script:
    """
    # Index reference if needed
    # minimap2 -d reference.mmi ${reference}
    # path "aligned_reads.bam.bai", emit: bai
    # path "aligned_reads.bam", emit: bam
    
    # Align reads
    minimap2 -ax map-ont ${reference} ${fastq} > aligned.sam \\
    
    #    samtools view -bS - | \\
    #    samtools sort -o aligned_reads.bam -
    
    # Index BAM file
    #samtools index aligned_reads.bam
    """
}

process INDEXING {
    container 'mgibio/samtools:v1.21-noble'
    publishDir "${params.outdir}/bam", mode: 'copy'
    
    input:
    path sam
    
    output:
    path "aligned_reads.bam", emit: bam
    path "aligned_reads.bai", emit: bai
    
    script:
    """
      # Run samtools for indexing sorting reads
      samtools view -bS ${sam} | samtools sort -o aligned_reads.bam
      samtools index -M aligned_reads.bam aligned_reads.bai
    """
}

