#!/bin/bash
##-------------------------------------------------------
## This script reads the formatted tsv file from interproscan (arg 1),
## and a list of gene names (arg 2) and reformats it into a geneinfo format
## The GOterm and IPRterm map files are required
##
## usage: ./generate_geneinfo.sh formatted_interproscan.tsv list_of_gene_names.txt
##-------------------------------------------------------

PROJECT_DIR=
OUTPUT_DIR=${PROJECT_DIR}/annotation_table

cd ${OUTPUT_DIR}

INTERPROSCAN_TSV=${OUTPUT_DIR}/$1
GENE_NAMES=${OUTPUT_DIR}/$2

INTERPROSCAN_FORMATTED="${INTERPROSCAN_TSV}.formatted"
GOMAP=${PROJECT_DIR}/goid_and_description.map
IPRMAP=${PROJECT_DIR}/interproid_and_description.map

GENEINFO="${INTERPROSCAN_FORMATTED}.geneinfo"



sed 's/_1*//2' ${INTERPROSCAN_TSV} | sed -e 's/-[0-9]//1' | sed -e 's/rna-//1' | sed -e 's/mrna_//1' | sed -e 's/mrna\.//1' | sed -e 's/\.t/\.g/1' > ${INTERPROSCAN_FORMATTED}


THREADS=1
## Convert interproscan to gene.info file
echo -e "
${SCRIPTS}/interproscan_to_geneinfo.R --interproscan ${INTERPROSCAN_FORMATTED}  --counts ${GENE_NAMES} --gomap ${GOMAP} --iprmap ${IPRMAP} --out ${GENEINFO}
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=50G -wd ${OUTPUT_DIR} -N interproscan_to_geneinfo








