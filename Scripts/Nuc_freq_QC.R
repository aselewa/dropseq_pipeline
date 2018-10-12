#!/bin/rscripts

#this script takes the fastq output file from read 1 to
#plot nucleotide frequencies for cell barcodes and UMIs
#both are on read1 and follow one another immediately.
#because the number of lines will be the same accross samples the line numbers for import
#are hard coded - change if adapted for different protocol

#use optparse for management of input arguments
library(optparse)
library(dplyr)
library(tidyr)
library(ggplot2)

option_list <- list(make_option(c("-f", "--file"),
                                type="character",
                                default=NULL,
                                help="fastqc file name",
                                metavar="character"),
                    make_option(c("-s", "--start"),
                                type="integer",
                                default=NULL,
                                help="lines to skip",
                                metavar="integer"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


#interrupt execution if no file is  supplied
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("fastqc file must be supplied.n", call.=FALSE)
}


#import and clean data
fq <- read.delim(file = opt$file,
                 skip = opt$start, #this needs to change for actual file
                 nrows = 20)
colnames(fq) <- c("Base_position", "G", "A", "T", "C")

fq_long <- fq[1:20,] %>%
  gather(., nuc, prop, G:C)


#get max and min values
min <- min(fq_long$prop) - 5
max <- max(fq_long$prop) + 5

#create plot for umi nucleotide diversity
fq_long %>%
  filter(., Base_position >= 13 & Base_position <= 20) %>%
  ggplot(., aes(Base_position - 12, prop, col = nuc)) +
    geom_point() +
    geom_line() +
    ggtitle("Nucleotide frequency in UMIs") +
    ylab("nucleotide proportion") +
    xlab("nucleotide position in sequence") +
    ylim(min, max)

ggsave(filename = "Nucleotide_frequency_in_UMIs.pdf" , path= "output/qc_data/", width = 8 , height = 8)


fq_long %>%
  filter(., Base_position <= 12) %>%
  ggplot(., aes(Base_position, prop, col = nuc)) +
    geom_point() +
    geom_line() +
    ggtitle("Nucleotide frequency in cell barcodes") +
    ylab("nucleotide proportion") +
    xlab("nucleotide position in sequence") +
    ylim(min, max)

ggsave(filename = "Nucleotide_frequency_in_cell_barcode.pdf" , path= "output/qc_data/", width = 8 , height = 8)
