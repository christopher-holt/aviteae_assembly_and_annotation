#!/bin/bash
#-------------------------------------------------------
## This script reads in a genome fasta file and a corresponding annotation file
## This script will run interproscan and generate a tsv formatted file
## How to run: ./generate_interproscan.sh genome_file.fasta genome_annotations.gff species_name
#-------------------------------------------------------
PROJECT_DIR=
PACKAGE_DIR=
USER_PACKAGE_DIR=

OUTPUT_DIR=${PROJECT_DIR}/annotation_statistics
REFERENCE_GENOME=${OUTPUT_DIR}/$1
NUCLEAR_GFF=${OUTPUT_DIR}/$2
SPECIES_NAME=$3

NUCLEAR_FNA="${SPECIES_NAME}.gene_sequences.fna"
NUCLEAR_AA="${SPECIES_NAME}.gene_sequences.aa"

GFFREAD_DIR=${USER_PACKAGE_DIR}/gffread-0.12.7
EMBOSS_BIN_DIR=${PACKAGE_DIR}/emboss-6.6.0/bin
INTERPROSCAN_DIR=${USER_PACKAGE_DIR}/interproscan-5.56-89.0

export PATH=${PACKAGE_DIR}/jdk-19/bin:$PATH


if [ ! -e ${HOME}/lib/libpcre.so ]; then
    ln -s /usr/lib64/libpcre.so.1.2.10 ${HOME}/lib/libpcre.so
fi
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}::${HOME}/lib

THREADS=15

## Generate FNA file
${GFFREAD_DIR}/gffread -x ${NUCLEAR_FNA} -g ${REFERENCE_GENOME} ${NUCLEAR_GFF}
## Convert FNA to Amino Acid
${EMBOSS_BIN_DIR}/transeq -sequence ${NUCLEAR_FNA} -outseq ${NUCLEAR_AA} -table 1 -clean

echo -e "
${INTERPROSCAN_DIR}/interproscan.sh -cpu ${THREADS} -i ${NUCLEAR_AA}  -d ${OUTPUT_DIR} -f tsv,gff3 -goterms -iprlookup -verbose
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=40G -wd ${OUTPUT_DIR} -N interpro.${SPECIES_NAME}

