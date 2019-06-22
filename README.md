# Disentangling-sRNA-Seq


Disentangling-sRNA-Seq consists of a series of scripts to separate small RNA sequencing data (sRNA-Seq) from libraries containing mixed information from two organisms.

Small RNA sequences organism assignment are particularly challenging relative to longer mRNA sequences (standard RNA-Seq protocols). Due to their short length (~20 nt) random mappings may occur in large genomes just by chance. 

In [1], we show that increasing sRNA sequences even by a single nucleotide decreases the chance for random mappings. We propose the usage of small RNA de novo or genome-guided assembly to increase sRNA length and by doing so disambiguate organism origin.

Unique mapping reads to a mixed reference, containing both organisms, serve as support counts that guide the assignment of multi-mapping reads. We take additional caution for those reads where no guidance is possible and keep them as ambiguous.

WARNING: Disentangling-sRNA-Seq is not meant to be a complete automated and user-friendly pipeline, it consists of a series of scripts that may be useful for the sRNA-mediated species communication community.  

Disentangling sRNA-Seq data to study RNA communication between species
bioRxiv 2019
JR Bermúdez-Barrientos, O Ramírez-Sánchez, FWN Chow, AH Buck, C Abreu-Goodger
doi: https://doi.org/10.1101/508937