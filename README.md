# Disentangling-sRNA-Seq

## Introduction
Assigning the correct organism to each small RNA read is particularly challenging relative to longer mRNA sequences (standard RNA-Seq protocols). Due to their short length (~20-25 nt) random mappings can occur in large genomes just by chance.
We show that increasing the length of sRNA sequences even by a few nucleotides can greatly decrease random mappings. We propose using de novo or genome-guided assembly of sRNA-Seq data, to increase sRNA length and help disambiguate the organism of origin.
When mapping reads to a mixed reference, containing both organisms, we take advantage of uniquely-mapping reads as support counts that guide the assignment of multi-mapping reads. For those reads where no guidance is possible, we take the conservative approach of keeping them as ambiguous.
We provide a pipeline called diSrna to help separate small RNA sequencing data (sRNA-Seq) from libraries containing information from two interacting organisms.

## Installation
Clone the repository. 

    git clone git@github.com:ObedRamirez/Disentangling-sRNA-Seq.git

Make sure that your system supplies the following dependencies for diSrna.
- OS: Linux, Mac OS
- miniconda
- bowtie-1.2.2
- bowtie2-2.4.2
- fastp
- tally
- samtools-1.5
- seqtk-trinity
- trinityrnaseq-v2.9.0
- pandoc
- Rscript
- The following R packages installed:
  - rmarkdown
  - rjson
  - ggplot2
  - RColorBrewer
  - kableExtra
  - gridExtra
  - edgeR
  - DT

Before you start:
1. You need to modify the diSrna script to tell it where the binary files of your tools are located. 
Example:
```
bin_dir = "/home/user/bin/"
bowtie_build_bin = bin_dir + "bowtie-1.2.2/bowtie-build"
```

If you have registered it in your path, you just have to name the tool:
    bowtie_build_bin = "bowtie-build"

2. The program assumes that: 
- Host and Par genomes (ej. M_musculus.fa/H_ polygyrus.fa) are located in the "genomes" directory (ej. /home/user/my_diSra_project/genomes).
- Your host libraries (fastq.gz files) are located in the "raw_data/host" directory (ej. /home/user/my_diSra_project/raw_data/host). 
- Your hostPar libraries (fastq.gz files) are located in the "raw_data/hostPar" directory (ej. /home/user/my_diSra_project/raw_data/hostPar).

## Use:
1. Activate conda environment
    conda activate snakemake-tutorial
2. make Bowtie1/Bowtie2 indexes
    snakemake -s ~/path_to/diSrna/diSrna --cores 6 bowtie_index_all
3. run diSrna pipe-line
    snakemake -s ~/path_to/diSrna/diSrna --cores 6 contigs_report
4. (Optional) Clean all directories created by diSrna (except: genomes/bowtieIndex, and dea_contigs) 
    snakemake -s ~/path_to/diSrna/diSrna clean_tmp

## Release description
We redesigned our workflow and implemented it using the Snakemake framework.
We updated most of our scripts to adapt them to the new workflow, but we didn't modify the core functions.
For this release, we chose Trinity de-novo assembly only.

## Citation
    If you use diSrna in your work, please cite the following:
Disentangling sRNA-Seq data to study RNA communication between species. 2020. JR Bermúdez-Barrientos, O Ramírez-Sánchez, FWN Chow, AH Buck, C Abreu-Goodger. Nucleic Acids Research. https://doi.org/10.1093/nar/gkz1198
