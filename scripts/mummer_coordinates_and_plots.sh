#!/bin/bash

#--------------------------------------------------------
# Set Usage Parametres
#--------------------------------------------------------
## This script will take two genome fasta files, align them
## with either mummer or promer, convert the delta file to coords file
## with the options -qlT, and then generate a png alignment plot
##
## How to run: ./mummer_coordinates_and_plots.sh -r /path/to/reference/fasta/file \
##                                               -q /path/to/query/fasta/file \
##                                               -a nucmer|promer \  
##                                               -p output_prefix \
##                                               -o /path/to/output/directory

programmename=$0
function usage {
    echo "$(basename $programmename) [-r <REFERENCE>] [-q <QUERY>] [-a <ALIGNMENT>] [-p <PREFIX>] [-o <OUTPUT_DIRECTORY>] [-h <help>]

    This script will generate nucmer or promer alignments between two genome fasta files. The two output files will be a mummerplot and a tab-deliminated coordinates file

    -r  Reference FASTA file, Required

    -q  Query FASTA file, Required

    -a  Alignment programe, must be either nucmer or promer, Required

    -p  File name prefix, Required

    -o  Output directory, Required

    ** FULL PATHS MUST BE SPECIFIED
    "
}

if [ "$#" -eq 0  ]; then
    usage
    exit
fi

if [ "$#" -eq 1 ] && [ "$1" == "-h" ]; then
    usage
    exit
fi


while getopts r:q:a:p:o:h flag
do 
    case ${flag} in
        q) REF=${OPTARG};;
        r) QUERY=${OPTARG};;
        a) ALIGNMENT=${OPTARG};;
        p) PREFIX=${OPTARG};;
        o) OUTPUT_DIR=${OPTARG};;
    esac
done


PACKAGE_DIR=${PKG}
MUMMER_DIR=${PACKAGE_DIR}/mummer-3.23

cd ${OUTPUT_DIR}

## Nucmer or promer alignment
${MUMMER_DIR}/${ALIGNMENT} -p ${PREFIX} ${REF} ${QUERY}

## Convert delta file to parsable coords file: -q = sort by query, -l = include sequence length -T = tab deliminated
${PACKAGE_DIR}/mummer-3.23/show-coords -qlT ${PREFIX}.delta > ${PREFIX}.coords

## Generate stock mummerplot
${MUMMER_DIR}/mummerplot --filter --color --layout --png -p ${PREFIX} ${PREFIX}.delta


