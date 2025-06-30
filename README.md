# Nanopore Analysis Pipeline

# FastFive
Putting the fast in FAST5 analysis, from alignment, to base, variant, and 5MC calling, and results.

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

## Introduction

This pipeline will perform a reasonable first-pass analysis of Oxford Nanopore FAST5 files including:

- **Basecalling** with Guppy
- **Quality Control** with NanoPlot
- **Read Alignment** with minimap2 followed by Samtools, for sorting the sam -> bam
- **Variant Calling** with Clair3
- **Methylation Analysis** with Nanopolish
- **Reporting and plots** with R

## Guide

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.0`)

2. Install [`Docker`](https://docs.docker.com/engine/installation/) also [`NVIDIA-Container Toolkit`](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html)


## Pipeline Summary

The pipeline performs the following steps:

1. **Basecalling** ([`Guppy`](https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revak_14dec2018))
2. **Quality Control** ([`NanoPlot`](https://github.com/wdecoster/NanoPlot))
3. **Read Alignment** ([`minimap2`](https://github.com/lh3/minimap2))
4. **Sam to Bam sorting and indexing** ([`samtools`](https://www.htslib.org/))
5. **Clair3 Variant calls** ([`Clair3`](https://github.com/HKU-BAL/Clair3/))
6. **Nanopolish for methylation calls** ([`Nanopolish`](https://github.com/jts/nanopolish))
7. **Reports and output using R** ([`R`](https://cran.r-project.org/))



## Usage

### Typical Command

```bash
nextflow run Hindrance/FastFive \
    --fast5 '/path/to/fast5/files/*.fast5' \
    --reference '/path/to/reference.fa' \
    --outdir './results' \
    -profile docker
```

### Example for Test Data
```bash
cd FastFive
unzip test_files.zip
nextflow run main.nf \
    --fast5 "test_files/fast5_files" \
    --reference "test_files/genomes/GRCh38/chr20_chr22_combined.fa" \
    --outdir "test_results" \
    -profile docker
```

### Input Requirements

- **FAST5 files**: Oxford Nanopore raw signal files
- **Reference genome**: FASTA format reference genome

### Core Nextflow Arguments

| Argument | Description |
|----------|-------------|
| `-profile` | Configuration profile (docker only currently) |
| `-resume` | Resume pipeline from last successful step |

### Custom Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--fast5` | Path to FAST5 input file(s) | None |
| `--reference` | Path to reference genome | None |
| `--outdir` | Output directory | `./results` |
| `--flowcell` | Flowcell type for basecalling | `FLO-MIN106` |
| `--kit` | Sequencing kit for basecalling | `SQK-LSK109` |
| `--sample_name` | Sample name for analysis | `sample001` |

## Output

The pipeline generates the following outputs in the specified output directory:

```
results/
├── basecalling/         # FASTQ files and sequencing summary
├── qc/                  # Initial basecall quality control reports and plots
├── alignment/           # SAM alignment files
├── bam/                 # BAM files and indices
├── variants/            # VCF files with variant calls
├── methylation/         # Methylation calling results
└── analysis_reports/    # Pipeline execution and summary reports?
```


## The future?

Future directions for this are:

1. Comprehensive summary report file (HTML for example the one that Claude threw together in misc_files)
2. Alternatively implement solutions offered in R shiny: ([`SCI-VCF`](https://github.com/HimanshuLab/SCI-VCF))
3. fix Test profile runs
4. Better error handling, error reporting.
5. Move to AWS!

