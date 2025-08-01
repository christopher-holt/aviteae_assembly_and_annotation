#!/usr/local/packages/r-4.1.2/bin/Rscript
#-------------------------------------------------------
## This script will generate a ggplot image
## of promer alignment coordinates
## How to run: ./plot_promer_coords.R --promer_coords /path/to/promer/coords --query Aviteae|Bmalayi \
##                                    --ref mitochondria|Wolbachia --outdir /path/to/output/directory
##
## usage: ./plot_promer_coords.R [-h] --promer_coords PROMER_COORDS
##                               --query QUERY --ref REF --outdir OUTDIR
##
## optional arguments:
##   -h, --help            show this help message and exit
##
## required arguments:
##   --promer_coords PROMER_COORDS
##                         Promer coords file
##   --query QUERY         Species used as the query, either mitochondria or
##                         Wolbachia
##   --ref REF             Species used as the reference, either Aviteae or
##                         Bmalayi
##   --outdir OUTDIR       Output Directory for Alignment Plots
##
#-------------------------------------------------------

## Load in argument parser
library(argparse)

parser <- argparse::ArgumentParser()
requiredArgs <- parser$add_argument_group("required arguments")
requiredArgs$add_argument("--promer_coords", help = "Promer coords file", required = T)
requiredArgs$add_argument("--query", help = "Species used as the query, either mitochondria or Wolbachia", required = T)
requiredArgs$add_argument("--ref", help = "Species used as the reference, either Aviteae or Bmalayi", required = T)
requiredArgs$add_argument("--outdir", help = "Output Directory for Alignment Plots", required = T)
args <- parser$parse_args()

## Loading in the rest of the packages
library(tidyverse)
library(ggplot2)


promer_plot_of_contigs<- function(df, query){
  
  
  colnames(df) <- c("QStart","QEnd","RefStart","RefEnd","QLength","RefLength","PercIdentity","QSize","RefSize","QContig","RefContig")
  
  
  if(query == "mitochondria"){
    query_length = 13723
    units = " (Kbp)"
    Mbp_factor = 1000000
    mito_scale_factor = 1000
    p <- ggplot() + 
      geom_segment(data = df,
                   mapping = aes(x = RefStart/Mbp_factor,
                                 xend = RefEnd/Mbp_factor,
                                 y = QStart/mito_scale_factor,
                                 yend = QEnd/mito_scale_factor,
                                 color=PercIdentity), linewidth = 3)
  }else if(query == "Wolbachia"){
    query_length = 1080084
    Mbp_factor = 1000000
    units = " (Mbp)"
    p <- ggplot() + 
      geom_point(data = df,
                 mapping = aes(x = RefStart/Mbp_factor,
                               y = QStart/Mbp_factor,
                               color=PercIdentity), linewidth = 3)
  }
  
  
  p <- p + 
    facet_grid(RefContig~.,scales="free",space="free") +
    theme_bw() + 
    theme(plot.title = element_text(size=14)) + 
    labs(y=paste0(unique(df$QContig), units), x = "A. viteae Contigs (Mbp)") +
    scale_color_viridis_c(option = "inferno",limits = c(58, 100)) + coord_flip() 

  print(
    p
  )
  
  return(p)
}

if(args$query == "Bmalayi"){
  ## Chr1, Chr2, Chr3, Chr4, ChrX
  chromosomes <- c("NW_025062475.1", "NW_025062476.1", "NW_025062477.1", "NW_025062478.1", "NW_025062479.1")
}else if(args$query == "Aviteae"){
  ## Chr1, Chr2, Chr3, Chr4, Chr5, ChrX
  chromosomes <- c("CM102052.1", "CM102053.1", "CM102054.1", "CM102055.1", "CM102056.1", "CM102057.1")
}


promer_coords <- as.data.frame(read.delim("/Volumes/scratch/chris.holt/aviteae_v_mito.coords", #args$promer_coords,
                                          sep = "\t", header = FALSE, stringsAsFactors = F, skip = 4))

colnames(promer_coords)<-str_replace_all(colnames(promer_coords), 'X', 'V')
promer_coords <- promer_coords[,c(1:7, 10, 11, 14, 15)]
colnames(promer_coords) <- paste0("V", seq(1, 11, by = 1))

for(chromosome in chromosomes){
  pdf(paste0(args$outdir,"/", args$query, "_vs_", args$ref, "_", chromosome, ".pdf"),
      width = 10,
      height = 10)
  print(promer_plot_of_contigs(promer_coords %>% filter(V11 == chromosome), args$query))
  dev.off()
}














