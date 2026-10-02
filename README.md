# Caenophidia_W_refugium
<img width="7275" height="6213" alt="Fig2_caenophidia_v9" src="https://github.com/user-attachments/assets/81fc3f37-5866-40fa-9c03-2f44c32939e3" />

This repository contains details on the data processing and analysis steps used in each aspect of this study. Here, we assemble a dataset of nine caenophidian snake species from publically available, published genomes on NCBI, and measure repeat elements across the genome for each species. We also utilize RNAseq data for three species within this dataset to measure expression of transposable elements and make comparisons between the sexes, and between chromosome types. This workflow is a companion to the methods described in Alderman et al. 2026. XXX.

Lists and reference files can be found in the resources directory. Shell and Python scripts are in shell and python directories, respectively. R scripts are in the R directory. Note that you may need to adjust the organization of path/file locations to suit your environment. This workflow assumes that you return to the main working directory after each major section.

Feel free to email me at hbs5rf[at]virginia.edu with any questions.

Visit the [Wiki](https://github.com/mgalderman/Caenophidia_W_refugium/wiki) page to see scripts and code associated with this project! 

## Contents

* Software and dependencies
* Chromosome identification - synteny and read-depth
* Pseudoautosomal region analyses
* Repeat annotation
* Chromosome repeat density
* Sliding-window repeats
* Kimura distance landscapes
* Full-length element annotation
* Refugium and toxicity index analyses
* TE expression for dysregulation

## Software and dependencies

The analysis sections below use the following software and dependencies and assume they are on the user path:

* FastQC
* MultiQC
* NCBI BLAST
* MashMap
* SRA toolkit
* Maker
* RepeatModeler2
* Repeatmasker
* trimmomatic
* bwa
* samtools
* bedtools
* mosdepth
* GenomeTools
* HMMER
* STAR
* sed
* featureCounts
* R
  
Note, I installed a number of these programs to my conda environment.
