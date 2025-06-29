# Nanopore Analysis Pipeline

# FastFive
Putting the fast in FAST5 analysis, from alignment, to base, variant, and 5MC calling, and results.

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

This pipeline performs comprehensive analysis of Oxford Nanopore FAST5 files including:

- **Basecalling** with Guppy
- **Quality Control** with NanoPlot
- **Read Alignment** with minimap2
- **Variant Calling** with Clair3
- **Methylation Analysis** with Nanopolish

## Guide

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.0`)

2. Install [`Docker`](https://docs.docker.com/engine/installation/) also (https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html)




