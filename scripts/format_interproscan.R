#!/usr/local/packages/r-4.1.2/bin/Rscript
##-------------------------------------------------------
## This script takes the output tsv file from interproscan
## and changes the transcript names to the gene names
##
## How to use: ./format_interproscan.R --interproscan /path/to/interproscan/tsv --map /path/to/rna/to/gene/name/conversion/table \
##                                     --output /path/to/reformatted/tsv/file 
##
## usage: ./format_interproscan.R [-h] --interproscan INTERPROSCAN --map MAP
##                                --output OUTPUT
##
## optional arguments:
##   -h, --help            show this help message and exit
##
## required arguments:
##   --interproscan INTERPROSCAN
##                         File path to the tsv formatted interproscan file
##   --map MAP             File path to the rna to gene conversion table
##   --output OUTPUT       Output file path to reformatted interproscan table
##-------------------------------------------------------

## Parsing Arguments
library(argparse)
parser = argparse::ArgumentParser()
requiredArgs <- parser$add_argument_group("required arguments")
requiredArgs$add_argument("--interproscan", help = "File path to the tsv formatted interproscan file", required = T)
requiredArgs$add_argument("--map", help = "File path to the rna to gene conversion table", required = T)
requiredArgs$add_argument("--output", help = "Output file path to reformatted interproscan table", required = T)
args <- parser$parse_args()


## Loading in other libraries
library(readr)
library(stringr)

interproscan_path <- args$interproscan
rna_to_gene_path <- args$map
fixed_interproscan_path <- args$output

interproscan <- readr::read_delim(interproscan_path,
                                  col_names = seq(1,15),
                                  delim = "\t",
                                  show_col_types = FALSE)

rna_to_gene <- readr::read_delim(rna_to_gene_path, delim = ';',
                                 col_names = F)


interproscan <- dplyr::inner_join(interproscan, rna_to_gene, by = "X1") %>% 
  select("X2.y", everything()) %>% 
  select(-"X1") %>% 
  dplyr::rename("X1" = "X2.y", "X2" = "X2.x")


readr::write_delim(interproscan, 
                   fixed_interproscan_path,
                   delim = "\t", col_names = F, quote = "none")
