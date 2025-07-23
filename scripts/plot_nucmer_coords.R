#!/usr/local/packages/r-4.1.2/bin/Rscript
#-------------------------------------------------------
## This script will generate a ggplot image
## of nucmer alignment coordinates
## How to run: ./plot_nucmer_coords.R --nucmer_coords /path/to/nucmer/coords --outdir /path/to/output/directory
##
## usage: ./plot_nucmer_coords.R [-h] --nucmer_coords NUCMER_COORDS
##                               --outdir OUTDIR
##
## optional arguments:
##   -h, --help            show this help message and exit
## 
## required arguments:
##   --nucmer_coords NUCMER_COORDS
##                         Nucmer Coords File
##   --outdir OUTDIR       Output Directory for Alignment Plots
## 
#-------------------------------------------------------

## Loading arguments
library(argparse)

parser <- argparse::ArgumentParser()
requiredArgs <- parser$add_argument_group("required arguments")
requiredArgs$add_argument("--nucmer_coords", help = "Nucmer Coords File", required = T)
requiredArgs$add_argument("--outdir", help = "Output Directory for Alignment Plots", required = T)
args <- parser$parse_args()


## Load the rest of the packages
library(tidyverse)
library(ggplot2)


## Function to check if the file exists
doesFileExist <- function(file){
  ## Check to see if a file exists, if not then quit with status = -1
  if (file.access(file) == -1){
    message(paste0(file, " does not exist"))
    q(status=1)
  }
}

nucmer_plot_of_contigs <- function(df){ ## Bm = Query, AV = Ref
  colnames(df) <- c("QStart","QEnd","RefStart","RefEnd","QLength","RefLength","PercIdentity","QSize","RefSize","QContig","RefContig")
  Mbp_factor = 1000000
  p <- ggplot() + 
    geom_segment(data = df,
                 mapping = aes(x = RefStart/Mbp_factor,
                               xend = RefEnd/Mbp_factor,
                               y = QStart/Mbp_factor,
                               yend = QEnd/Mbp_factor,
                               color=PercIdentity), linewidth = 3) + 
    facet_grid(RefContig~.,scales="free",space="free") + 
    theme_bw() + 
    theme(plot.title = element_text(size=14)) + 
    labs(y=paste0(unique(df$QContig), " (Mbp)"), x = "A. viteae Contigs (Mbp)") + 
    scale_color_viridis_c(option = "inferno",limits = c(58, 100)) + coord_flip()
  print(
    p
  )
  
  return(p)
}


## Checking validity of files and create the output directory
doesFileExist(args$nucmer_coords)
dir.create(outdir, showWarnings = F)

## Reading in the coordinates files
Bm_vs_AV <- as.data.frame(read.delim(paste0(args$nucmer_coords), 
                                        sep = "\t", header = FALSE, stringsAsFactors = F, skip = 4))
colnames(Bm_vs_AV)<-str_replace_all(colnames(Bm_vs_AV), 'X', 'V')

## Bm Chr 1 = AV Chr 1
## Bm Chr 3 = AV Chr 2
## Bm Chr 2 = AV Chr 3
## Bm Chr 4 = AV Chr 5
## Bm Chr X = AV Chr 4 and AV Chr X
## Bm Mito = AV Mito

chromosomes <- list("NW_025062475.1" = "CM102052.1",
                    "NW_025062477.1" = "CM102053.1", 
                    "NW_025062476.1" = "CM102054.1", 
                    "NW_025062478.1" = "CM102056.1",
                    "NW_025062479.1" = c("CM102055.1", "CM102057.1"),
                    "NC_004298.1" = "CM102097.1")

for(i in 1:length(names(chromosomes))){
  
  QContig <- names(chromosomes)[i]
  RefContig <- chromosomes[i][[1]]
  
  df <- Bm_vs_AV %>% filter(V10 == QContig & V11 == RefContig)
  
  
  if(length(RefContig) > 1){
    RefContig_1 <- chromosomes[i][[1]][1]
    RefContig_2 <- chromosomes[i][[1]][2]
   
    df <- Bm_vs_AV %>% filter(V10 == QContig & (V11 == RefContig_1 | V11 == RefContig_2))
  }
  
  pdf(paste0(args$outdir, "/", QContig, "_", RefContig, ".pdf"),
      width = 10,
      height = 10)
  print(nucmer_plot_of_contigs(df))
  dev.off()
  

}


