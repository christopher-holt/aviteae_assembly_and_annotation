#!/usr/local/packages/r-4.1.2/bin/Rscript
#-------------------------------------------------------
## This script reads in a GFF file and will print out:
##      The number of genes with multiple transcripts
##      Total number of genes
##      Average number of exons per gene
##      Number of inferred exons per gene
##      Average CDS length
##      Average exon length
##
## How to run: ./annotation_stats.R --annotation /path/to/gff/file
##
## usage: ./annotation_stats.R [-h] --annotation ANNOTATION
##
## optional arguments:
##   -h, --help            show this help message and exit
##
## required arguments:
##   --annotation ANNOTATION
##                         GFF3 Annotation File
#-------------------------------------------------------

## Parsing Arguments
library(argparse)

parser <- argparse::ArgumentParser()
requiredArgs <- parser$add_argument_group("required arguments")
requiredArgs$add_argument("--annotation", help = "GFF3 Annotation File", required = T)

args <- parser$parse_args()


## Loading in libraries
library(readr)
library(dplyr)
library(tidyverse)

doesFileExist <- function(file){
  ## Check to see if a file exists, if not then quit with status = -1
  if (file.access(file) == -1){
    message(paste0(file, " does not exist"))
    q(status=1)
  }
}

format_gff_file <- function(annotations){
  annotations <- annotations %>% separate('X9', into = c("Name", "X9"), sep = ';')
  annotations$X9 <- gsub(".*:", "", annotations$X9)
  annotations$X9 <- gsub(".*=", "", annotations$X9)
  annotations$Name <- gsub(".*:", "", annotations$Name)
  annotations$Name <- gsub(".*=", "", annotations$Name)
  
  return(annotations)
}

num_rows <- function(annotations, feature){
  return(annotations %>% filter(X3 == paste0(feature)) %>% nrow(.))
}

num_feature_per_gene <- function(annotations, feature){
  (num_rows(annotations, paste0(feature)))/num_rows(annotations, 'gene')
}

average_length <- function(annotations,feature){
  return((annotations %>% filter(X3 == paste0(feature)) %>% mutate(length = X5 - X4) %>% 
            summarise(x = sum(length)))/(annotations %>% filter(X3 == paste0(feature)) %>% nrow(.)))
}

summary_statistics <- function(annotations){
  message('Number of Genes with Multiple Transcripts: ', as.data.frame(table(annotations[annotations$X3 == "mRNA",]$X10)) %>% filter(Freq >1) %>% nrow(.))
  message('total number of genes: ', num_rows(annotations, 'gene'))
  message('Avg Number Exon/gene: ', num_feature_per_gene(annotations, 'exon'))
  message('Number of Inferred Introns per Gene: ', num_feature_per_gene(annotations, "exon") - num_feature_per_gene(annotations, 'mRNA'))
  message('Avg CDS Length: ', average_length(annotations, 'CDS'))
  ## Avg Exon size
  message('Avg Exon Length: ', average_length(annotations, 'exon'))
  
}


doesFileExist(args$annotation)

ANNOTATION_FILE <- args$annotation

## Read in Annotation File
annotations <- readr::read_delim(ANNOTATION_FILE, "\t", escape_double = FALSE, col_names = FALSE,
                                 trim_ws = TRUE, show_col_types = F, comment = "#")
annotations <- format_gff_file(annotations)


mRNA_features <- annotations %>% filter(X3 == "mRNA")
mRNA_features <- unique(mRNA_features$X9)

annotations <- annotations %>% mutate(X10 = ifelse(.$X3 == "gene", .$Name, NA))
annotations <- tidyr::fill(annotations, X10)
annotations <- annotations[annotations$X10 %in% mRNA_features,]

summary_statistics(annotations)



