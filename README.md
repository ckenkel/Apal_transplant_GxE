# Additive effects of transplant site and host genotype influence growth, survival and microbiome diversity in Acropora palmata
This study quantifies the relative effects of genotype and outplant site on coral growth and survival. Three replicate ramets of five genotypes of common garden reared sexually produced Acropora palmata were transplanted to each of nine offshore reef sites in the Lower Florida Keys. Growth, bleaching status, and survival were monitored with 3D photogrammetry every three months for two years. Tissue samples taken at the initial time point and two years post-transplantation were used to quantify shifts in gene expression, endosymbiont community, and microbiome composition. This repository contains the bioinformatic and statistical scripts necessary to re-create analyses of survival, growth and microbiome diversity, as well as raw input data used to generate figures. Raw FASTQ reads of 16S amplicons can be obtained from NCBI’s SRA under BioProject NNN. 

Files in this repository 
-----------

1. bbduk.sh: shell script to remove 16S amplicon reads containing Illumina sequencing adapters. Run notes can be found as preamble in dadaseq2_16S_ApalmHE.R
   
2. dadaseq2_16S_ApalmHE.R: Annotated R script for quality filtering raw reads, idenetifying ASVs and assigning taxonomy. This workflow follows the DADA2 tutorial available at https://benjjneb.github.io/dada2/tutorial.html and can be run on a standard desktop/laptop computer 
	- Input Files
	  - Adapter filtered paired end FASTQ read files (R1 and R2), one per sample in a single directory. 
	- Output Files (input for the subsequent phyloseq script)
	  - Apal_seqtabNochim_9May24.rds: ASV table
	  - Apal_taxonomy_9May24.rds: Associated taxonomy for ASVs based on Silva v138.1
   
3. Phyloseq_apalm16S.R: Annotated R script for statistical analyses of microbial diversity. Follows aspects of the phyloseq and ancombc tutorials available at https://joey711.github.io/phyloseq/ and https://bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html 
	- Input Files
	  - Apal_seqtabNochim_9May24.rds
    - Apal_taxonomy_9May24.rds
	  - sampledata.csv: A csv file containing sample IDs and associated metadata (e.g. site, genet ID, etc)
    
4. ApalmMorphologicalPlasticity.R: Annotated R script for analyses of survival and growth
	- Input file: ApalmMorphologicalPlasticity_rready.csv
