#!/bin/bash


#takes the fastqc_data.txt file and determines the line above the section summarizing nucleotide
#frequencies
skip_line=`less $1 | awk '$1 ~/>>Per/ && $2 ~/base/ && $3 ~/sequence/ && $4 ~/content/ {print(NR)}'`

#reads in file and the number of lines to skip and returns two plots of nucleotide frequencies
#one corresponding to cell barcodes and one to UMIs
Rscript --vanilla /project2/gilad/spott/Pipelines/dropseq_pipeline/Scripts/Nuc_freq_QC.R -f $1 -s ${skip_line}
