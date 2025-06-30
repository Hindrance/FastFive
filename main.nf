#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Parameters
params.fast5 = null
params.reference = null
params.outdir = "results"
params.help = false

// Help message
def helpMessage() {
    log.info"""
    Usage:
    nextflow run FastFive.nf --fast5 <path_to_fast5> --reference <path_to_reference_genome>

    Required arguments:
    --fast5         Path to input FAST5 file
    --reference         Path to input genome fasta file
    
    
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
if (!params.fast5 | !params.reference) {
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
    SAM2BAM(ALIGNMENT.out.sam)    
    
    // clair3 variant calling
    VARIANT_CALLING(SAM2BAM.out.bam, SAM2BAM.out.bai, ch_reference)
    
     // Methylation calling with Nanopolish
    METHYLATION_CALLING(SAM2BAM.out.bam, SAM2BAM.out.bai, ch_reference, ch_fast5, BASECALLING.out.fastq)
    
    // Reports and fun plot times with R
    ANALYSIS_REPORTS(METHYLATION_CALLING.out.methyl_tsv, ch_reference, ch_fast5)
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
        --prefix ""
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
    # Align reads
    minimap2 -ax map-ont ${reference} ${fastq} > aligned.sam \\
    

    """
}

process SAM2BAM {
    container 'mgibio/samtools:v1.21-noble'
    publishDir "${params.outdir}/bam", mode: 'copy'
    
    input:
    path sam
    
    output:
    path "aligned_reads.bam", emit: bam
    path "aligned_reads.bam.bai", emit: bai
    
    script:
    """
      # Run samtools for indexing sorting reads
      samtools view -bS ${sam} | samtools sort -o aligned_reads.bam
      samtools index aligned_reads.bam aligned_reads.bam.bai
    """
}

process VARIANT_CALLING {
    container 'hkubal/clair3:latest'
    publishDir "${params.outdir}/variants", mode: 'copy'
    
    input:
    path bam
    path bai
    path reference
    
    output:
    path "merge_output.vcf.gz"
    path "run_clair3.log"

    
    
    script:
    """
    # quickly index the genome FA
    samtools faidx ${reference}
    
    # Run Clair3 for variant calling
    # apparently Clair3 hates white space...
    /opt/bin/run_clair3.sh \\
        --bam_fn=${bam} \\
        --ref_fn=${reference} \\
        --threads=${task.cpus} \\
        --platform=ont \\
        --model_path=/opt/models/ont_guppy5 \\
        --output=. \\
        --sample_name=sample001
    """
}


process METHYLATION_CALLING {
    container 'quay.io/biocontainers/nanopolish:0.14.0--hee927d3_5'
    publishDir "${params.outdir}/methylation", mode: 'copy'
   
    input:
    path bam
    path bai
    path reference
    path fast5
    path fastq
    
    output:
    path "methylation_calls.tsv", emit: methyl_tsv
    
    script:
    """
    # Index the fast5 files
    nanopolish index -d . ${fastq}
    
    # Call methylation
    nanopolish call-methylation \\
        --reads=${fastq} \\
        --bam=${bam} \\
        --genome=${reference} \\
        --t=1 > methylation_calls.tsv
    """
}

process ANALYSIS_REPORTS {
    container 'r-base:latest'
    publishDir "${params.outdir}/analysis_reports", mode: 'copy'
   
    input:
    path methyl_tsv
    path fast5
    path reference
    
    output:
    path "*.png"
    
    script:
    """
    # quick R script for fun...    
    cat > analyze_data.R << 'EOF'
    #!/usr/bin/env R
    quickPlot = function(methyl_path){
      genome_df = data.frame(
        "chromosome" = as.character(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22, "x", "y")),
        "min" = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
        "max" = c(247000000,242951149,199501827,191273063,180857866,170899992,158821424,146274826,140273252,135374737,
        134452384,132349534,114142980,106368585,100338915,88827254,78774742,76117153,63811651,62435964,46944323,49691432,154913754,57772954
        )
      )

      #x = read.csv("methylation_calls.tsv", sep="\t")
      x = read.csv(methyl_path, sep="\t")
      x_meth = x[x[,"log_lik_ratio"] > 0,]

      bw = 500000

      png("genome_pileup_plot.png", width=1200, height=1000)
      par(cex=1.5)
      plot(c(0, max(genome_df[,"max"])), c(1,24), pch=NA, xlab="bases", ylab="chromosome", yaxt="n")
      axis(2, at=1:24, label=c("x", "y", as.character(22:1)))#

      for(chr in 1:nrow(genome_df)){
        lines(c(genome_df[chr,"min"], genome_df[chr,"max"]), rep(25-chr,2), lwd=1)
        
        if(chr %in% unique(as.character(x[,"chromosome"]))){
        
          freq_bins = seq(genome_df[chr,"min"], genome_df[chr,"max"], bw)

          freq_table = sapply(1:(length(freq_bins)-1), function(i){sum(x[,4] > freq_bins[i] & x[,4] <= freq_bins[i+1])})
          freq_table2 = (freq_table / max(freq_table))*0.8
          
          
          
          for(i in 1:length(freq_bins)){
            rect(freq_bins[i], (25-chr)+0, freq_bins[i+1], (25-chr)+freq_table2[i+1], col="grey", lty=0)
          }
          
          
          freq_table_meth = sapply(1:(length(freq_bins)-1), function(i){sum(x_meth[,4] > freq_bins[i] & x_meth[,4] <= freq_bins[i+1])})
          freq_table_meth2 = (freq_table_meth / max(freq_table))*0.8
          
          for(i in 1:length(freq_bins)){
            rect(freq_bins[i], (25-chr)+0, freq_bins[i+1], (25-chr)+freq_table_meth2[i+1], col="red", lty=0)
          }
        }
      }
      dev.off()
    }
    quickPlot("${methyl_tsv}")
EOF
    Rscript analyze_data.R
    """
}



