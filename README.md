# Disentangling-sRNA-Seq


We provide a series of scripts to help separate small RNA sequencing data (sRNA-Seq) from libraries containing information from two interacting organisms.

Assigning the correct organism to each small RNA read is particularly challenging relative to longer mRNA sequences (standard RNA-Seq protocols). Due to their short length (~20-25 nt) random mappings can occur in large genomes just by chance. 

We show that increasing the length of sRNA sequences even by a few nucleotides can greatly decrease random mappings. We propose using *de novo* or genome-guided assembly of sRNA-Seq data, to increase sRNA length and help disambiguate the organism of origin.

When mapping reads to a mixed reference, containing both organisms, we take advantage of uniquely-mapping reads as support counts that guide the assignment of multi-mapping reads. For those reads where no guidance is possible, we take the conservative approach of keeping them as ambiguous.

WARNING: The scripts contained in this repository are not meant to be a completely automated and user-friendly pipeline: they they are provided "as is" with no warranty of any kind. We used these scripts for all the analyses performed in the following publication:

Disentangling sRNA-Seq data to study RNA communication between species. 2019. JR Bermúdez-Barrientos, O Ramírez-Sánchez, FWN Chow, AH Buck, C Abreu-Goodger. bioRxiv. https://doi.org/10.1101/508937
