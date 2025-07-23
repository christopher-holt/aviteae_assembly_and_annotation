#!/usr/local/packages/r-4.1.2/bin/Rscript

#-------------------------------------------------------
## This script will plot the output or
## median sequencing depth from samtools depth -aa
## using a 10kb moving window
##
## How to run: ./plot_sequencing_depth.R --depthfile /path/to/depth/file --outdir /path/to/output/directory
##
## usage: ./scripts/plot_sequencing_depth.R [-h] [--median_sequencing_depth]
##                                          --outdir OUTDIR --depthfile DEPTHFILE
##
## optional arguments:
##   -h, --help            show this help message and exit
##   --median_sequencing_depth
##                         Print out the median sequencing depth for each contig
##                         and do not plot the depth, Optional
##
## required arguments:
##  --outdir OUTDIR       Output Directory, Required
##  --depthfile DEPTHFILE
##                        Samtools Depth File, Required
#-------------------------------------------------------

## Parsing arguments
library(argparse)
### Optional arguments
parser <- argparse::ArgumentParser()
parser$add_argument("--median_sequencing_depth", action = "store_true", 
                    help = "Print out the median sequencing depth for each contig and do  not plot the depth, Optional")

## Required arguments
requiredArgs <- parser$add_argument_group("required arguments")
requiredArgs$add_argument("--outdir", help = "Output Directory", required = T) 
requiredArgs$add_argument("--depthfile", help = "Samtools Depth File", required = T)
args <- parser$parse_args()

## Plot depth over genome
library(tidyverse)

doesFileExist <- function(file){
  ## Check to see if a file exists, if not then quit with status = -1
  if (file.access(file) == -1){
    message(paste0(file, " does not exist"))
    q(status=1)
  }
}

doesFileExist(args$depthfile)

DEPTH_TYPE <- args$depthfile
outdir <- args$outdir

dir.create(paste0(outdir), showWarnings = F)

###Read in the depth file
large_contig_window = 10000
Mbp = 1000000
WormHist <- read_delim(paste0(DEPTH_TYPE), 
                       "\t", escape_double = FALSE, col_names = FALSE, 
                       trim_ws = TRUE)
colnames(WormHist) <- str_replace_all(colnames(WormHist), 'X', 'V')

colnames(WormHist) <- c("Contig", "Position", "Depth")

WormHistRelevant <- WormHist

NigonList<-unique(WormHist$Contig)
Graphing<-data.frame()

### Plot depth
for (a in 1:length(NigonList)){

  WormHistSub <- WormHistRelevant[WormHistRelevant$Contig == NigonList[a],]
  
  window = large_contig_window
  
  
  
  PartialGraphing<-data.frame(matrix(ncol=3,nrow=ceiling(length(WormHistSub$Contig)/window)),stringsAsFactors = FALSE)
  colnames(PartialGraphing)<-c("Position","Avg","Nigon")
  for (b in 1:ceiling(length(WormHistSub$Contig)/window))
  {
    if(b < ceiling(length(WormHistSub$Contig)/window)){
      avg <- mean(WormHistSub[((b*window)-(window-1)):(b*window),]$Depth)
      PartialGraphing[b,]$Position<-b*window
      PartialGraphing[b,]$Avg<-avg
      PartialGraphing[b,]$Nigon<-NigonList[a]
    } else {
      endpoint<-length(WormHistSub$Contig)
      avg<-mean(WormHistSub[((b*window)-(window-1)):endpoint,]$Depth)
      PartialGraphing[b,]$Position<-endpoint
      PartialGraphing[b,]$Avg<-avg
      PartialGraphing[b,]$Nigon<-NigonList[a]
    }
  }
  Graphing <- rbind(Graphing,PartialGraphing)
}

## Rename columns and include start position
colnames(Graphing)<-c("End_Position","Sequencing_Depth","Nigon")


Graphing <- Graphing %>% mutate(Start_Position = Graphing$End_Position - (large_contig_window-1))


Graphing <- Graphing %>% dplyr::select("Nigon", "Start_Position", "End_Position", "Sequencing_Depth")

## This prints median for Table 4.2, have to include other sequencing technologies
if(args$median_sequencing_depth){
  Graphing %>% group_by(Nigon) %>%
    summarise(median(Sequencing_Depth))
}else{
  pdf(paste0(outdir, "/", basename(DEPTH_TYPE), ".pdf"),
      width = 30)
  plot(
    ggplot(Graphing %>% filter(grepl("Chr", Graphing$Nigon))) + 
      geom_point(aes(x=Start_Position/Mbp,y=Sequencing_Depth), size = 0.5) + 
      ggtitle(paste0("Sequencing Depth for ", str_split(DEPTH_TYPE, "\\.")[[1]][1] ,sep="")) + 
      labs(x = "Start_Position (Mbp)") +
      theme_bw() + 
      theme(legend.position="none",plot.title = element_text(size=14)) +
      facet_grid(~Nigon)
    
  )
  dev.off() 
}

write_delim(Graphing, 
            file = paste0(outdir, "/", basename(DEPTH_TYPE), ".histogram.csv"), 
            delim = "\t",
            col_names = T)










