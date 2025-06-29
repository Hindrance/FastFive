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
      
      # Align reads
      minimap2 \\
        -ax map-ont \\
        ${reference} \\
        ${fastq} > \\
        aligned.sam
    """
}

