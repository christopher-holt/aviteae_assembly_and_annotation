# Genome Assembly and Annotation of *Acanthocheilonema viteae* FR3 Strain

## Table Of Contents
* [Genome Assembly](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#genome-assembly)
    * [Canu Assembly](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#canu-assembly)
    * [Hifi Polishing with Racon](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#hifi-polishing-with-racon)
    * [Mitochondria Assembly](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#mitochondria-assembly)
* [Genome Annotation](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#genome-annotation)
    * [Annotate and Mask Repeats](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#annotate-and-mask-repeats)
    * [Align Reads](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#align-reads)
    * [Generate Structural Annotation](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#generate-structural-annotation)
        * [Generate List of Bam Files](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#generate-list-of-bam-files)
        * [Run Braker](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#run-braker)
        * [Combine Annotations](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#combine-annotations)
    * [Generate Updated GFF File](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#generate-updated-gff-file)
    * [Generate Functional Annotation](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#generate-functional-annotation)
    * [Add Locus Tags to GFF File](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#add-locus-tags-to-gff-file)
    * [Generate New Polypeptide File](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#generate-new-polypeptide-file)
    * [Check for Genbank Compatability](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#check-for-genbank-compatability)
* [Generate Geneinfo](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#generate-geneinfo)
    * [Download GO and IPR Terms](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#download-go-and-ipr-terms)
    * [Run InterproScan](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#run-interproscan)
    * [Generate Counts](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#generate-counts)
    * [Convert InterproScan to GeneInfo](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#convert-interproscan-to-geneinfo)
* [Mummerplot](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#mummerplot)
* [Comparison Tables](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#comparison-tables)
    * [Assembly](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#assembly)
    * [Annotation](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#annotation)
        * [InterproScan and Geneinfo](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#interproscan-and-geneinfo)
* [Sequencing Depth](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#sequencing-depth)
    * [Generate Depth File](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#generate-depth-file)
        * [Inputs](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#inputs)
        * [Commands](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#commands)
            * [Illumina](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#illumina)
            * [CLR, HiFi, and ONT Long Reads](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#clr-hifi-and-ont-long-reads)
    * [Plot Sequencing Depth](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#plot-sequencing-depth)

## Genome Assembly
### Canu Assembly
```bash
#!/bin/bash

## Set Prefix and Output dirs
PROJECT_DIR=
PACKAGE_DIR=
FASTQ_DIR=${PROJECT_DIR}/fastq

CANU_BIN_DIR=${PACKAGE_DIR}/canu-2.1.1/bin

## Create assembly using ONT data

OUTPUT_PREFIX=aviteae.ont
OUTPUT_DIR=${PROJECT_DIR}/aviteae_assembly/canu.ont.min106
ONT_FASTQ=${FASTQ_DIR}/DNA_female_ONT_min106.fastq.gz

mkdir -p ${OUTPUT_DIR}

## Original run had genome size of 90m
${CANU_BIN_DIR}/canu -p ${OUTPUT_PREFIX} -d ${OUTPUT_DIR} \
    genomeSize=90m \
    gridEngineThreadsOption="-pe thread THREADS" \
    gridEngineMemoryOption="-l mem_free=MEMORY" \
    gridOptions="-P jdhotopp-lab -q threaded.q" \
    -nanopore ${ONT_FASTQ}

```
### Hifi Polishing with Racon
```bash
PROJECT_DIR=
PACKAGE_DIR=

RACON_BIN_DIR=${PACKAGE_DIR}/racon-1.3.1/bin
SAMTOOLS_BIN_DIR=${PACKAGE_DIR}/samtools-1.11/bin
MINIMAP2_BIN_DIR=${PACKAGE_DIR}/minimap2-2.24/bin/

FASTQ_DIR=${PROJECT_DIR}/fastq
ASSEMBLY_DIR=${PROJECT_DIR}/aviteae_assembly
OUTPUT_DIR=${ASSEMBLY_DIR}/polishing/racon_polishing

HIFI_READS=${FASTQ_DIR}/DNA_female_1_HiFi.fastq.gz

mkdir -p ${OUTPUT_DIR}

THREADS=10


if ! [[ -L ${OUTPUT_DIR}/aviteae.ont.contigs.racon.000.fasta" && -e ${OUTPUT_DIR}/aviteae.ont.contigs.racon.000.fasta" ]]; then
    ln -s ${ASSEMBLY_DIR}/FINAL_CONTIGS/Aviteae_Potential_Final_Nuclear_Genome.fasta ${OUTPUT_DIR}/aviteae.ont.contigs.racon.000.fasta
fi

## Starts indexing at 0, so do 1 round of polishing
for round in {0..0}; do

    polishing_round=$(printf "%03d" $round)
    new_fasta=$(printf "%03d" $((round+1)))
    old_fasta=$(printf "%03d" $((round-1)))

    CANU_ASSEMBLY=${OUTPUT_DIR}/aviteae.ont.contigs.racon.${polishing_round}.fasta
    SAM=${OUTPUT_DIR}/aviteae.ont.contigs.racon.${polishing_round}.sam

    echo -e "
    ${MINIMAP2_BIN_DIR}/minimap2 -a ${CANU_ASSEMBLY} ${HIFI_READS} | gzip -1 > ${SAM}
    " | qsub -V -q threaded.q -pe thread ${THREADS} -P jhotopp-gcid-proj4b-filariasis -N minimap_overlaps.${polishing_round} -wd ${OUTPUT_DIR} -l mem_free=20G -hold_jid racon_correction.${old_fasta}

    echo -e "
    ${RACON_BIN_DIR}/racon ${HIFI_READS} ${SAM} ${CANU_ASSEMBLY} > ${OUTPUT_DIR}/aviteae.ont.contigs.racon.${new_fasta}.fasta &&
        rm ${SAM}
    " | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=100G -wd ${OUTPUT_DIR} -N racon_correction.${polishing_round} -hold_jid minimap_overlaps.${polishing_round}

done

```
### Mitochondria Assembly
```bash
PROJECT_DIR=
PACKAGE_DIR=
PYTHON3_BIN_DIR=${PACKAGE_DIR}/python-3.8.2/bin/
FASTQ_DIR=${PROJECT_DIR}/fastq

MITO_ASSEMBLY_DIR=${PROJECT_DIR}/aviteae_assembly/mitochondria_assembly/

ONT_READS=${FASTQ}/DNA_female_ONT_min106.fastq.gz
UNICYCLER_ASSEMBLY_DIR=${MITO_ASSEMBLY_DIR}/unicycler_assembly

mkdir -p ${MITO_ASSEMBLY_DIR} ${ASSEMBLY_DIR} ${UNICYCLER_ASSEMBLY_DIR}

THREADS=10

## Long Read Only
echo -e "
${PYTHON3_BIN_DIR}/unicycler -t ${THREADS} -l ${ONT_READS} -o ${UNICYCLER_ASSEMBLY_DIR} --keep 0
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=70G -wd ${MITO_ASSEMBLY_DIR} -N unicycler_assembly 
```

## Genome Annotation 
### Annotate and Mask Repeats
```bash
PACKAGE_DIR=
PROJECT_DIR=

ANNOTATION_DIR=${PROJECT_DIR}/annotation_final_rc/annotate_and_mask_repeats
REFERENCE_DIR=${PROJECT_DIR}/aviteae_assembly_final/polishing_final/racon_polishing/

REPEATMODELER_OUTPUT_DIR=${ANNOTATION_DIR}/repeatmodeler
REPEATMASKER_SOFTMASKED_OUTPUT_DIR=${ANNOTATION_DIR}/repeatmasker_softmasked
REPEATMASKER_MASKED_OUTPUT_DIR=${ANNOTATION_DIR}/repeatmasker_masked

REPEATMASKER_DIR=${PACKAGE_DIR}/repeatmasker-4.0.7
REPEATMODELER_BIN_DIR=${PACKAGE_DIR}/repeatmodeler-1.0.11/

AVITEAE_REFERENCE=${REFERENCE_DIR}/aviteae_racon_polished_final_rc.fasta

mkdir -p ${REPEATMODELER_OUTPUT_DIR} ${REPEATMASKER_SOFTMASKED_OUTPUT_DIR} ${REPEATMASKER_MASKED_OUTPUT_DIR}

THREADS=8

## Build Database
echo -e "
${REPEATMODELER_BIN_DIR}/BuildDatabase -name ${REPEATMODELER_OUTPUT_DIR}/aviteae.ncbi.db ${AVITEAE_REFERENCE} -engine ncbi -dir ${REPEATMODELER_OUTPUT_DIR}
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -N build_database -l mem_free=30G -wd ${REPEATMODELER_OUTPUT_DIR} -pe thread ${THREADS} -q threaded.q

## Model Repeats
echo -e "
${REPEATMODELER_BIN_DIR}/RepeatModeler -pa ${THREADS} -database ${REPEATMODELER_OUTPUT_DIR}/aviteae.ncbi.db -engine ncbi -dir ${REPEATMODELER_OUTPUT_DIR}
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -N repeatmodeler -l mem_free=50G -wd ${REPEATMODELER_OUTPUT_DIR} -pe thread ${THREADS} -q threaded.q -hold_jid build_database

## Generate SoftMasked Fasta
echo -e "
${REPEATMASKER_DIR}/RepeatMasker -pa ${THREADS} -xsmall -gff -e ncbi -lib ${REPEATMODELER_OUTPUT_DIR}/consensi.fa.classified -s -dir ${REPEATMASKER_SOFTMASKED_OUTPUT_DIR} ${AVITEAE_REFERENCE}
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -N soft_repeatmasker -l mem_free=50G -wd ${REPEATMASKER_SOFTMASKED_OUTPUT_DIR} -pe thread ${THREADS} -q threaded.q -hold_jid repeatmodeler

## Generate HardMasked Fasta
echo -e "
${REPEATMASKER_DIR}/RepeatMasker -pa ${THREADS} -gff -e ncbi -lib ${REPEATMODELER_OUTPUT_DIR}/consensi.fa.classified -s -dir ${REPEATMASKER_MASKED_OUTPUT_DIR} ${AVITEAE_REFERENCE}
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -N masked_repeatmasker -l mem_free=50G -wd ${REPEATMASKER_MASKED_OUTPUT_DIR} -pe thread ${THREADS} -q threaded.q -hold_jid repeatmodeler

```

### Align Reads
```bash
PROJECT_DIR=
PACKAGE_DIR=
ILLUMINA_LIST=${PROJECT_DIR}/annotation_final_rc/srr.id.list
THREADS=10

READS_DIR=${PROJECT_DIR}/fastq/
TRIMMED_READS_DIR=${READS_DIR}/trimmed_reads

JDK_DIR=${PACKAGE_DIR}/jdk/bin/
TRIMMOMATIC_DIR=${PACKAGE_DIR}/trimmomatic-0.38/
HISAT2_DIR=${PACKAGE_DIR}/hisat2-2.1.0
SAMTOOLS_BIN_DIR=${PACKAGE_DIR}/samtools-1.9/bin

ANNOTATION_DIR=${PROJECT_DIR}/annotation_final_rc/annotate_and_mask_repeats/repeatmasker_softmasked
REFERENCE=${ANNOTATION_DIR}/aviteae_racon_polished_final_rc.fasta.masked
INDEXED_REFERENCE=${ANNOTATION_DIR}/aviteae_racon_polished_final_rc.fasta.masked.indexed

BAM_DIR=${PROJECT_DIR}/annotation_final_rc/bam


mkdir -p ${TRIMMED_READS_DIR} ${BAM_DIR}


echo -e " 
${HISAT2_DIR}/hisat2-build ${REFERENCE} ${INDEXED_REFERENCE} 
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=20G -wd ${REFERENCE_DIR} -N hisat.build


for SRR in $(cat ${ILLUMINA_LIST}); do
    SAMPLE=$(echo ${SRR} | cut -d',' -f3)
    FASTQ1=$(echo ${READS_DIR}/${SAMPLE}_R1.fastq.gz)
    FASTQ2=$(echo ${READS_DIR}/${SAMPLE}_R2.fastq.gz)
    SORTED_BAM_FILE=${BAM_DIR}/${SAMPLE}.sorted.bam 

    ## Trim Reads
    echo -e "
    ${JDK_DIR}/java -Xmx20g -jar ${TRIMMOMATIC_DIR}/trimmomatic-0.38.jar PE -threads ${THREADS} -phred33 -trimlog ${TRIMMED_READS_DIR}/trim_RNA.log \
        ${FASTQ1} ${FASTQ2} ${TRIMMED_READS_DIR}/${SAMPLE}_paired_R1.fastq.gz ${TRIMMED_READS_DIR}/${SAMPLE}_unpaired_R1.fastq.gz ${TRIMMED_READS_DIR}/${SAMPLE}_paired_R2.fastq.gz  \
        ${TRIMMED_READS_DIR}/${SAMPLE}_unpaired_R2.fastq.gz \
        ILLUMINACLIP:${TRIMMOMATIC_DIR}/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
    " | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=30G -wd ${TRIMMED_READS_DIR} -N trimmomatic.${SAMPLE}

    ## Align Reads
    echo -e "
    ${HISAT2_DIR}/hisat2 -p ${THREADS} --rna-strandness RF --max-intronlen 50000 -x ${INDEXED_REFERENCE} -1 ${TRIMMED_READS_DIR}/${SAMPLE}_paired_R1.fastq.gz -2 ${TRIMMED_READS_DIR}/${SAMPLE}_paired_R2.fastq.gz | \
        ${SAMTOOLS_BIN_DIR}/samtools sort -@ ${THREADS} -o ${SORTED_BAM_FILE} && \
       ${SAMTOOLS_BIN_DIR}/samtools index ${SORTED_BAM_FILE}
    " | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=20G -wd ${BAM_DIR} -N hisat2.align.all.${SAMPLE} -hold_jid trimmomatic.${SAMPLE},hisat.build

done

```

### Generate Structural Annotation
#### Generate List of Bam Files
```bash
PROJECT_DIR=
BAM_DIR=${PROJECT_DIR}/annotation_final_rc/bam

RNA_FEMALE=${BAM_DIR}/RNA_female.sorted.bam
RNA_L3=${BAM_DIR}/RNA_L3.sorted.bam
RNA_MALE=${BAM_DIR}/RNA_male.sorted.bam
RNA_MF=${BAM_DIR}/RNA_mf.sorted.bam


echo -e "${RNA_FEMALE},${RNA_L3},${RNA_MALE},${RNA_MF}" > ${PROJECT_DIR}/list_of_bam_files
```


#### Run Braker
```bash
### Set Directories and File Paths
PROJECT_DIR=
PACKAGE_DIR=
USER_PACKAGE_DIR=

REFERENCE_DIR=${PROJECT_DIR}/annotation_final_rc/annotate_and_mask_repeats/repeatmasker_softmasked
REFERENCE_FASTA=${REFERENCE_DIR}/aviteae_racon_polished_final_rc.fasta.masked

## Custom Installs
### Genemark
GENEMARK_DIR=${USER_PACKAGE_DIR}/genemark-es-et-4.62 ## Needed to edit perl shebang line for some files


## Set Paths for programs and files
BRAKER_SCRIPT_DIR=${PACKAGE_DIR}/braker-2.1.6/scripts

PROTEIN_OUTPUT_DIR=${PROJECT_DIR}/annotation_final_rc/structural_annotation/braker_protein.final
RNASEQ_OUTPUT_DIR=${PROJECT_DIR}/annotation_final_rc/structural_annotation/braker_rnaseq.final

AUGUSTUS_PKG_DIR=${PACKAGE_DIR}/augustus-3.4.0/
AUGUSTUS_BIN_DIR=${AUGUSTUS_PKG_DIR}/bin/
AUGUSTUS_SCRIPTS_DIR=${AUGUSTUS_PKG_DIR}/scripts/
BAMTOOLS_BIN_DIR=${PACKAGE_DIR}/bamtools-2.5.1/bin/
PYTHON3_BIN_DIR=${PACKAGE_DIR}/python-3.8.2/bin/
DIAMOND_PATH=${USER_PACKAGE_DIR}/diamond-0.9.24/
CDBFASTA_BIN_DIR=${PACKAGE_DIR}/cdbfasta/
SAMTOOLS_BIN_DIR=${PACKAGE_DIR}/samtools-1.9/bin/
BLAST_BIN_DIR=${PACKAGE_DIR}/ncbi-blast+-2.8.1/bin/
PROTHINT_BIN_DIR=${GENEMARK_DIR}/ProtHint-2.6.0/bin


SPECIES=Acanthocheilonemaviteae

export LD_LIBRARY_PATH="${PACKAGE_DIR}/gsl-2.4/lib:$LD_LIBRARY_PATH"

## List of RNA-Seq Bam files for braker
BAM_FILES=${PROJECT_DIR}/annotation_final_rc/list_of_bam_files


THREADS=25

## Make sure gm_key is present
if [ ! -f ${HOME}/.gm_key ]; then
 	ln -s ${GENEMARK_DIR}/gm_key ${HOME}/.gm_key
fi

### Set Augustus Dirs with Permissions in Working Dir
mkdir -p ${PROTEIN_OUTPUT_DIR} ${RNASEQ_OUTPUT_DIR}


## RNA-Seq Structrual Annotation
AUGUSTUS_WD_DIR=${RNASEQ_OUTPUT_DIR}/augustus-3.4.0/
AUGUSTUS_CONFIG_PATH=${AUGUSTUS_WD_DIR}/config/
mkdir -p ${AUGUSTUS_WD_DIR}
if [ ! -d ${AUGUSTUS_WD_DIR}/config ]; then

    cp -r ${AUGUSTUS_PKG_DIR}/config ${AUGUSTUS_WD_DIR}
    chmod -R a+rwx ${RNASEQ_OUTPUT_DIR}/augustus-3.4.0

fi

echo -e " 
${BRAKER_SCRIPT_DIR}/braker.pl --species=${SPECIES} \
    --genome=${REFERENCE_FASTA} \
    --bam=$(cat ${BAM_FILES}) \
    --cores=${THREADS} \
    --AUGUSTUS_CONFIG_PATH=${AUGUSTUS_CONFIG_PATH} \
    --AUGUSTUS_BIN_PATH=${AUGUSTUS_BIN_DIR} \
    --AUGUSTUS_SCRIPTS_PATH=${AUGUSTUS_SCRIPTS_DIR} \
    --workingdir=${RNASEQ_OUTPUT_DIR} \
    --BAMTOOLS_PATH=${BAMTOOLS_BIN_DIR} \
    --PYTHON3_PATH=${PYTHON3_BIN_DIR} \
    --GENEMARK_PATH=${GENEMARK_DIR} \
    --DIAMOND_PATH=${DIAMOND_PATH} \
    --CDBTOOLS_PATH=${CDBFASTA_BIN_DIR} \
    --SAMTOOLS_PATH=${SAMTOOLS_BIN_DIR} \
    --useexisting \
    --gff3 \
    --softmasking 
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -N ${SPECIES}.rnaseq_braker -l mem_free=100G -wd ${RNASEQ_OUTPUT_DIR} -pe thread ${THREADS} -q threaded.q




## Protein Structural Annotation
AUGUSTUS_WD_DIR=${PROTEIN_OUTPUT_DIR}/augustus-3.4.0/
if [ ! -d ${AUGUSTUS_WD_DIR}/config ]; then

    cp -r ${AUGUSTUS_PKG_DIR}/config ${AUGUSTUS_WD_DIR}
    chmod -R a+rwx ${PROTEIN_OUTPUT_DIR}/augustus-3.4.0

fi

echo -e " 
${BRAKER_SCRIPT_DIR}/braker.pl --species=${SPECIES} \
    --genome=${REFERENCE_FASTA} \
    --cores=${THREADS} \
    --AUGUSTUS_CONFIG_PATH=${AUGUSTUS_CONFIG_PATH} \
    --AUGUSTUS_BIN_PATH=${AUGUSTUS_BIN_DIR} \
    --AUGUSTUS_SCRIPTS_PATH=${AUGUSTUS_SCRIPTS_DIR} \
    --workingdir=${PROTEIN_OUTPUT_DIR} \
    --BAMTOOLS_PATH=${BAMTOOLS_BIN_DIR} \
    --PYTHON3_PATH=${PYTHON3_BIN_DIR} \
    --GENEMARK_PATH=${GENEMARK_DIR} \
    --DIAMOND_PATH=${DIAMOND_PATH} \
    --CDBTOOLS_PATH=${CDBFASTA_BIN_DIR}  \
    --BLAST_PATH=${BLAST_BIN_DIR} \
    --SAMTOOLS_PATH=${SAMTOOLS_BIN_DIR} \
    --PROTHINT_PATH=${PROTHINT_BIN_DIR} \
    --prot_seq=${PROJECT_DIR}/databases/proteins.fasta \
    --epmode \
    --useexisting \
    --gff3 \
    --softmasking 
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -N ${SPECIES}.protein_braker -l mem_free=50G -wd ${PROTEIN_OUTPUT_DIR} -pe thread ${THREADS} -q threaded.q



```
#### Combine Annotations
```bash
PROJECT_DIR=
USER_PACKAGE_DIR=

TSEBRA_DIR=${USER_PACKAGE_DIR}/TSEBRA-1.0.3/
ANNOTATION_DIR=${PROJECT_DIR}/annotation_final_rc/structural_annotation
PROT_ONLY=${ANNOTATION_DIR}/braker_protein.final/
RNA_SEQ_ONLY=${ANNOTATION_DIR}/braker_rnaseq.final

${TSEBRA_DIR}/bin/tsebra.py -g ${RNA_SEQ_ONLY}/augustus.hints.gtf,${PROT_ONLY}/augustus.hints.gtf -c ${TSEBRA_DIR}/config/default.cfg -e ${RNA_SEQ_ONLY}/hintsfile.gff,${PROT_ONLY}/hintsfile.gff -o ${ANNOTATION_DIR}/aviteae_structural_annotation.final.gtf

```

### Generate Updated GFF File
```bash
SCRIPT_DIR=
PROJECT_DIR=
PACKAGE_DIR=
USER_PACKAGE_DIR=

ANNOTATION_DIR=${PROJECT_DIR}/annotation_final_rc/structural_annotation/
AVITEAE_AA=${ANNOTATION_DIR}/aviteae_structural_annotation.final.aa
GTF=${ANNOTATION_DIR}/aviteae_structural_annotation.final.gtf
GFF=${ANNOTATION_DIR}/aviteae_structural_annotation.final.gff3
REFERENCE_DIR=${PROJECT_DIR}/aviteae_assembly_final/polishing_final/racon_polishing/
AVITEAE_FASTA=${REFERENCE_DIR}/aviteae_racon_polished_final_rc.fasta

THREADS=1

echo -e "
${PACKAGE_DIR}/r-4.1.2/bin/Rscript ${SCRIPT_DIR}/convert_tsebra_gtf_to_gff.R --gtf ${GTF} --gff ${GFF}
" | qsub -V -q threaded.q -pe thread ${THREADS} -P jhotopp-gcid-proj4b-filariasis -N convert_file_format -wd ${ANNOTATION_DIR} -l mem_free=50G



echo -e "
${USER_PACKAGE_DIR}/gffread-0.12.7/gffread -y ${AVITEAE_AA} -g ${AVITEAE_FASTA} ${GFF}
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=50G -wd ${ANNOTATION_DIR} -N gff_read -hold_jid convert_file_format

```

### Generate Functional Annotation
```bash
SCRIPT_DIR=
PROJECT_DIR=

ANNOTATION_DIR=${PROJECT_DIR}/annotation_final_rc
AA_HINTS=${ANNOTATION_DIR}/structural_annotation/aviteae_structural_annotation.final.aa
AVITEAE_ANNOTATION=${ANNOTATION_DIR}/structural_annotation/aviteae_structural_annotation.final.gff3
FUNCTIONAL_ANNOTATION_DIR=${ANNOTATION_DIR}/functional_annotation/tsebra.final
REFERENCE_DIR=${PROJECT_DIR}/aviteae_assembly_final/polishing_final/racon_polishing/
AVITEAE_FASTA=${REFERENCE_DIR}/aviteae_racon_polished_final_rc.fasta


sh ${SCRIPT_DIR}/tsebra_functional_annotation.sh \
    -a ${AA_HINTS} \
    -s ${AVITEAE_ANNOTATION} \
    -o ${FUNCTIONAL_ANNOTATION_DIR} \
    -f ${AVITEAE_FASTA} \
    -g ${FUNCTIONAL_ANNOTATION_DIR}/aviteae_functional_annotation.final_rc.gff3
```

### Add Locus Tags to GFF File
The ACH3XW tag was provided by NCBI
```bash
PROJECT_DIR=
ANNOTATION_DIR=${PROJECT_DIR}/annotation_final_rc/functional_annotation/tsebra.final
## BIOCODE_DIR obtained from https://github.com/jorvis/biocode
BIOCODE_DIR=


${BIOCODE_DIR}/gff/add_gff3_locus_tags.py -i ${ANNOTATION_DIR}/aviteae_functional_annotation.final.gff3 -o ${ANNOTATION_DIR}/aviteae_functional_annotation.final.with_locus_tags.gff3 -p ACH3XW -a 4 -n 5 -s 5
```

### Generate New Polypeptide File
```bash
PROJECT_DIR=
USER_PACKAGE_DIR=
PACKAGE_DIR=

ANNOTATION_DIR=${PROJECT_DIR}/annotation_final_rc/functional_annotation/tsebra.final

AVITEAE_AA=${ANNOTATION_DIR}/aviteae_functional_annotation.final_rc.with_locus_tags.aa
AVITEAE_FNA=${ANNOTATION_DIR}/aviteae_functional_annotation.final_rc.with_locus_tags.fna
AVITEAE_RESIDUES=${ANNOTATION_DIR}/aviteae_functional_annotation.final_rc.with_locus_tags.residues
AVITEAE_GFF=${ANNOTATION_DIR}/aviteae_functional_annotation.final_rc.with_locus_tags.gff3
AVITEAE_FASTA=${PROJECT_DIR}/aviteae_assembly_final/polishing_final/racon_polishing/aviteae_racon_polished_final_rc.fasta

THREADS=1

echo -e "
${USER_PACKAGE_DIR}/gffread-0.12.7/gffread -x ${AVITEAE_FNA} -g ${AVITEAE_FASTA} ${AVITEAE_GFF}
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=50G -wd ${ANNOTATION_DIR} -N gff_read -hold_jid convert_file_format

echo -e "
${PACKAGE_DIR}/emboss-6.6.0/bin/transeq -sequence ${AVITEAE_FNA} -outseq ${AVITEAE_AA} -table 1 -clean
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=50G -wd ${ANNOTATION_DIR} -N emboss -hold_jid gff_read
```

### Check for Genbank Compatability
```bash
THREADS=1

PROJECT_DIR=
USER_PACKAGE_DIR=
GB_DIR=${PROJECT_DIR}/annotation_final_rc/fix_for_ncbi

SBT=${GB_DIR}/Aviteae.sbt
CMT=${GB_DIR}/Aviteae.cmt
FSA=${GB_DIR}/Aviteae_nucgenome_polished_rc.fsa ## Symlink of ${PROJECT_DIR}/aviteae_assembly_final/polishing_final/racon_polishing/aviteae_racon_polished_final_rc.fasta
GFF3=${PROJECT_DIR}/functional_annotation/tsebra.final/aviteae_functional_annotation.final_rc.with_locus_tags.gff3

mkdir -p ${GB_DIR}/table2asn_output
cd ${GB_DIR}/table2asn_output
${USER_PACKAGE_DIR}/tbl2asn/linux64.table2asn_GFF -J -c w -t ${SBT} -w ${CMT} -i ${FSA} -f ${GFF3} -j "[gcode=1]" -M n -V vb -Z -euk -outdir ${GB_DIR}/table2asn_output
```

## Generate GeneInfo
The code below uses the reference and annotation file from [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_046563165.1/)
### Download GO and IPR Terms
```bash
## Download GO Terms
PACKAGE_DIR=
PROJECT_DIR=
${PACKAGE_DIR}/r-4.1.2/bin/Rscript ${SCRIPTS}/download_go_terms.R --output ${PROJECT_DIR}/goid_and_description.map

## Download IPR Terms
PYTHON_BIN_DIR=${PACKAGE_DIR}/python-3.8.2/bin

## Code is available at https://www.ebi.ac.uk/interpro/result/download/#/entry/InterPro/|tsv
${PYTHON_BIN_DIR}/python download_ipr_terms.py > ${PROJECT_DIR}/interproid_and_description.map

```
### Run InterproScan
```bash
## Run Interproscan
PROJECT_DIR=
PACKAGE_DIR=
USER_PACKAGE_DIR=
INTERPROSCAN_DIR=${USER_PACKAGE_DIR}/interproscan-5.56-89.0
GFFREAD_DIR=${USER_PACKAGE_DIR}/gffread-0.12.7
OUTPUT_DIR=${PROJECT_DIR}/interproscan
AVITEAE_AA=${PROJECT_DIR}/reference/GCA_046563165.1_ASM4656316v1_genomic.aa
GFF=${PROJECT_DIR}/reference/GCA_046563165.1_ASM4656316v1_genomic.gff
REFERENCE_GENOME=${PROJECT_DIR}/reference/GCA_046563165.1_ASM4656316v1_genomic.fna



mkdir -p ${OUTPUT_DIR}


if [ ! -e ${HOME}/lib/libpcre.so ]; then
    ln -s /usr/lib64/libpcre.so.1.2.10 ${HOME}/lib/libpcre.so
fi
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}::${HOME}/lib
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
export PATH=${PACKAGE_DIR}/jdk-19/bin:$PATH

THREADS=15

${GFFREAD_DIR}/gffread -y ${AVITEAE_AA} -g ${REFERENCE_GENOME} ${GFF}

echo -e "
${INTERPROSCAN_DIR}/interproscan.sh -cpu ${THREADS} -i ${AVITEAE_AA}  -d ${OUTPUT_DIR} -f tsv,gff3 -goterms -iprlookup -verbose
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=40G -wd ${OUTPUT_DIR} -N interpro

```
### Generate Counts
Using bam files generated in [Align Reads](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#align-reads)
```bash
PACKAGE_DIR=
PROJECT_DIR=

PYTHON_BIN_DIR=${PACKAGE_DIR}/python-3.8.2/bin

THREADS=16

## Set Directories and Output Files
BAM_DIR=${PROJECT_DIR}/bam_final
COUNTS_DIR=${PROJECT_DIR}/counts
COMBINED_COUNTS=${COUNTS_DIR}/combined.final.counts
SRR_LIST=${PROJECT_DIR}/annotation_final_rc/srr.id.list
REFERENCE_DIR=${PROJECT_DIR}/reference
GFF3=${REFERENCE_DIR}/GCA_046563165.1_ASM4656316v1_genomic.gff

mkdir -p ${COUNTS_DIR}

for SAMPLE in $(cat ${SRR_LIST}); do
    SAMPLE_NAME=$(echo ${SAMPLE} | cut -d',' -f3)
    SORTED_BAM_FILE=${BAM_DIR}/${SAMPLE_NAME}.sorted.bam 

    echo -e "
    ${PYTHON_BIN_DIR}/htseq-count -n ${THREADS} -s reverse --max-reads-in-buffer 3000000000 -r pos --nonunique none -f bam -m union -t gene --idattr ID ${SORTED_BAM_FILE} ${GFF3} \
    | awk -v a=${SAMPLE} '{print \$0, a}' | sed -e 's/ /\t/g' >> ${COMBINED_COUNTS}
    " | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=40G -wd ${OUTPUT_DIR} -N generate_counts.${SAMPLE_NAME}
done
```


### Convert InterproScan to GeneInfo
```bash
PROJECT_DIR=
COUNTS_DIR=${PROJECT_DIR}/counts
OUTPUT_DIR=${PROJECT_DIR}/interproscan
INTERPROSCAN_TSV=${OUTPUT_DIR}/GCA_046563165.1_ASM4656316v1_genomic.aa.tsv
INTERPROSCAN_FORMATTED=${OUTPUT_DIR}/GCA_046563165.1_ASM4656316v1_genomic.aa.tsv.fixed
GOMAP=${PROJECT_DIR}/goid_and_description.map
IPRMAP=${PROJECT_DIR}/interproid_and_description.map
COMBINED_COUNTS=${COUNTS_DIR}/combined.final.counts

## Create rna to gene conversion table
cat ${OUTPUT_DIR}/GCA_046563165.1_ASM4656316v1_genomic.gff | awk -F'\t' '$3=="mRNA"' | cut -f9 | cut -d';' -f1,2 | sed -e 's/ID=//g' | sed -e 's/Parent=//g' > ${OUTPUT_DIR}/rna_to_gene.txt
## Swap transcript for gene names
${SCRIPTS}/format_interproscan.R --interproscan ${INTERPROSCAN_TSV} --map ${OUTPUT_DIR}/rna_to_gene.txt --output ${INTERPROSCAN_FORMATTED}


THREADS=1
## Convert interproscan to gene.info file
echo -e "
${SCRIPTS}/interproscan_to_geneinfo.R --interproscan ${INTERPROSCAN_FORMATTED} --counts ${COMBINED_COUNTS} --gomap ${GOMAP} --iprmap ${IPRMAP} --out ${OUTPUT_DIR}/aviteae_gene.info
" | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=50G -wd ${OUTPUT_DIR} -N interproscan_to_geneinfo

```
## Mummerplot
### Inputs
```bash
PROJECT_DIR=
OUTPUT_DIR=${PROJECT_DIR}/mummerplot

mkdir -p ${OUTPUT_DIR}

## Nucmer matches of A. viteae and B. malayi Nuclear Genome
QUERY=${PROJECT_DIR}/reference/GCF_000002995.4_B_malayi-4.0_genomic.fna
REF=${PROJECT_DIR}/reference/GCA_046563165.1_ASM4656316v1_genomic.fna
PREFIX=bmalayi_vs_aviteae
ALIGNMENT=nucmer

## Promer matches of A. viteae and its mitochondria
QUERY=${PROJECT_DIR}/reference/GCA_046563165.1_ASM4656316v1_genomic_mitochondria.fna
REF=${PROJECT_DIR}/reference/GCA_046563165.1_ASM4656316v1_genomic.fna
PREFIX=mito_vs_aviteae_promer
ALIGNMENT=promer

## Promer matches of B. malayi and its mitochondria
QUERY=${PROJECT_DIR}/reference/GCF_NC004298.1_B_malayi_mitochondria_genomic.fna
REF=${PROJECT_DIR}/reference/GCF_000002995.4_B_malayi-4.0_genomic.fna
PREFIX=mito_vs_bmalayi_promer
ALIGNMENT=promer

## Promer matches of A. viteae and wBm (B. malayi mitochondria)
QUERY=${PROJECT_DIR}/reference/b_malayi.AE017321.genomic.fa
REF=${PROJECT_DIR}/reference/GCA_046563165.1_ASM4656316v1_genomic.fna
PREFIX=wbm_vs_aviteae_promer
ALIGNMENT=promer

## Promer matches of B. malayi and wBm
QUERY=${PROJECT_DIR}/reference/b_malayi.AE017321.genomic.fa
REF=${PROJECT_DIR}/reference/GCF_000002995.4_B_malayi-4.0_genomic.fna
PREFIX=wbm_vs_bmalayi_promer
ALIGNMENT=promer
```

### Command
These mummerplots are shown in Fig 4.1 and Figs 4.4-4.7 respectively
```bash
sh ${SCRIPTS}/mummer_coordinates_and_plots.sh -r ${REF} -q ${QUERY} -a ${ALIGNMENT}$ -p ${PREFIX}$ -o ${OUTPUT_DIR}

## For Aviteae and Bmalayi, ggplot of the nucmer alignments
${RSCRIPT_BIN_DIR}/Rscript ${SCRIPTS}/plot_nucmer_coords.R --nucmer_coords ${OUTPUT_DIR}/bmalayi_vs_aviteae.coords --outdir ${OUTPUT_DIR}
## Plotting promer coordinates
### Aviteae vs mitochondria
${RSCRIPT_BIN_DIR}/Rscript ${SCRIPTS}/plot_promer_coords.R --promer_coords ${OUTPUT_DIR}/mito_vs_aviteae_promer.coords --query Aviteae --ref mitochondria --outdir ${OUTPUT_DIR}
### Bmalayi vs mitochondria
${RSCRIPT_BIN_DIR}/Rscript ${SCRIPTS}/plot_promer_coords.R --promer_coords ${OUTPUT_DIR}/mito_vs_bmalayi_promer.coords --query Bmalayi --ref mitochondria --outdir ${OUTPUT_DIR}
### Aviteae vs wBm
${RSCRIPT_BIN_DIR}/Rscript ${SCRIPTS}/plot_promer_coords.R --promer_coords ${OUTPUT_DIR}/wbm_vs_aviteae_promer.coords --query Aviteae --ref Wolbachia --outdir ${OUTPUT_DIR}
### Bmalayi vs wBm
${RSCRIPT_BIN_DIR}/Rscript ${SCRIPTS}/plot_promer_coords.R --promer_coords ${OUTPUT_DIR}/wbm_vs_bmalayi_promer.coords --query Bmalayi --ref Wolbachia --outdir ${OUTPUT_DIR}
```

## Comparison Tables
### Assembly
#### Inputs
```bash
PROJECT_DIR=
PACKAGE_DIR=
USER_PACKAGE_DIR=

BUSCO_BIN_DIR=${USER_PACKAGE_DIR}/busco/bin/
QUAST_DIR=${USER_PACKAGE_DUR}/quast
BUSCO_OUTPUT_DIR=${PROJECT_DIR}/busco_final
QUAST_OUTPUT_DIR=${PROJECT_DIR}/quast
THREADS=15

mkdir -p ${BUSCO_OUTPUT_DIR} ${QUAST_OUTPUT_DIR}

source ${PACKAGE_DIR}/miniconda3/etc/profile.d/conda.sh
conda activate busco

AUGUSTUS_WD_DIR=${PROJECT_DIR}/busco_final/augustus-3.4.0/
AUGUSTUS_CONFIG_PATH=${AUGUSTUS_WD_DIR}/config/
AUGUSTUS_SPECIES_PATH=${AUGUSTUS_WD_DIR}/config/species

export LD_LIBRARY_PATH="${PACKAGE_DIR}/gsl-2.4/lib:$LD_LIBRARY_PATH"
export AUGUSTUS_CONFIG_PATH=${AUGUSTUS_CONFIG_PATH}
export AUGUSTUS_SPECIES_PATH=${AUGUSTUS_SPECIES_PATH}

## A. viteae
ASSEMBLY=${PROJECT_DIR}/reference/GCA_046563165.1_ASM4656316v1_genomic.fna
CONFIG=${OUTPUT_DIR}/config_aviteae.ini

## A. viteae (Blaxter)
ASSEMBLY=${PROJECT_DIR}/reference/Aviteae_GCA_900537255.1_ASM90053725v1_genomic.fna
CONFIG=${OUTPUT_DIR}/config_aviteae_blaxter.ini

## B. malayi
ASSEMBLY=${PROJECT_DIR}/reference/GCF_000002995.4_B_malayi-4.0_genomic.gff
CONFIG=${OUTPUT_DIR}/config_bmalayi.ini

## O. volvulus
ASSEMBLY=${PROJECT_DIR}/reference/ovolvulus.PRJEB513.WBPS15.genomic.fa
CONFIG=${OUTPUT_DIR}/config_ovolvulus.ini

## D. immitis
ASSEMBLY=${PROJECT_DIR}/reference/Dimmitis_GCA_024305405.1_ICBAS_JMDir_1.0_genomic.fna
CONFIG=${OUTPUT_DIR}/config_dimmitis.ini

## Li. sigmodontis
ASSEMBLY=${PROJECT_DIR}/reference/Lsigmodontis_GCA_963070105.1_nxLitSigm11.1_genomic.fna
CONFIG=${OUTPUT_DIR}/config_lsigmodontis.ini

## Loa loa
ASSEMBLY=${PROJECT_DIR}/reference/Lloa_GCA_000183805.3_Loa_loa_V3.1_genomic.fna
CONFIG=${OUTPUT_DIR}/config_loa.ini

```
#### Commands
The BUSCO scores and quast outputs are in Table 4.1
```bash
## BUSCO, need a new output dir for each assembly
cd ${BUSCO_OUTPUT_DIR}
${BUSCO_BIN_DIR}/busco -c ${THREADS} -f --config ${CONFIG} -m genome -i ${ASSEMBLY} --download_path ${BUSCO_OUTPUT_DIR}
## Quast, need a new output dir for each assembly
cd ${QUAST_OUTPUT_DIR}
${QUAST_DIR}/quast.py -t ${THREADS} -o ${QUAST_OUTPUT_DIR} ${REFERENCE_FILE}
```


### Annotation
#### Inputs
```bash
PROJECT_DIR=
PACKAGE_DIR=

OUTPUT_DIR=${PROJECT_DIR}/annotation_statistics
R_SCRIPT_BIN_DIR=${PACKAGE_DIR}/r-4.1.2/bin

## A. viteae 
GFF3=${PROJECT_DIR}/reference/GCA_046563165.1_ASM4656316v1_genomic.gff
ANNOTATION_STATS_OUTPUT=${OUTPUT_DIR}/aviteae_statistics.txt

## A. viteae (Blaxter)
GFF3=${PROJECT_DIR}/reference/Aviteae_GCA_900537255.1_ASM90053725v1_genomic.gff
ANNOTATION_STATS_OUTPUT=${OUTPUT_DIR}/aviteae_blaxter_statistics.txt

## B. malayi
GFF3=${PROJECT_DIR}/reference/GCF_000002995.4_B_malayi-4.0_genomic.gff
ANNOTATION_STATS_OUTPUT=${OUTPUT_DIR}/bmalayi_statistics.txt

## O. volvulus 
GFF3=${PROJECT_DIR}/reference/ovolvulus.PRJEB513.WBPS15.annotations.no.comment.no_mitochondria.wormbase.gff3
ANNOTATION_STATS_OUTPUT=${OUTPUT_DIR}/ovolvulus_statistics.txt

## D. immitis
GFF3=${PROJECT_DIR}/reference/Dimmitis_GCA_024305405.1_ICBAS_JMDir_1.0_genomic.gff
ANNOTATION_STATS_OUTPUT=${OUTPUT_DIR}/dimmitis_statistics.txt

## Li. sigmodontis
GFF3=${PROJECT_DIR}/reference/litomosoides_sigmodontis.PRJEB3075.WBPS19.annotations.gff3
ANNOTATION_STATS_OUTPUT=${OUTPUT_DIR}/lsigmodontis_statistics.txt

## Loa loa
GFF3=${PROJECT_DIR}/reference/Lloa_GCA_000183805.3_Loa_loa_V3.1_genomic.gff
ANNOTATION_STATS_OUTPUT=${OUTPUT_DIR}/lloa_statistics.txt
```

#### Commands
The annnotation statistics are shown in Table 4.3
```bash
${RSCRIPT_BIN_DIR}/Rscript ${SCRIPTS}/annotation_stats.R --annotation ${GFF3} > ${ANNOTATION_STATS_OUTPUT}
```
### InterproScan and Geneinfo
The *A. viteae* InterproScan and Geneinfo results were generated in [Run InterproScan](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#run-interproscan) and [Convert InterproScan to GeneInfo](https://github.com/christopher-holt/aviteae_assembly_and_annotation?tab=readme-ov-file#convert-interproscan-to-geneinfo). Geneinfo files are available in [input_files/annotation_table](https://github.com/christopher-holt/aviteae_assembly_and_annotation/tree/main/input_files/annotation_table). The *B. malayi* geneinfo is available [here](https://raw.githubusercontent.com/matt-chung/lf_transcriptome/refs/heads/master/input_data_files/bmalayi_gene.info).

#### InterproScan
```bash
## D. immitis
${SCRIPTS}/generate_interproscan.sh Dimmitis_GCA_024305405.1_ICBAS_JMDir_1.0_genomic.fna Dimmitis_GCA_024305405.1_ICBAS_JMDir_1.0_genomic.gff Dimmitis
## O. volvulus
${SCRIPTS}/generate_interproscan.sh ovolvulus.PRJEB513.WBPS15.genomic.fa ovolvulus.PRJEB513.WBPS15.annotations.no.comment.no_mitochondria.wormbase.gff3 Ovolvulus
## Loa loa
${SCRIPTS}/generate_interproscan.sh Lloa_GCA_000183805.3_Loa_loa_V3.1_genomic.fna Lloa_GCA_000183805.3_Loa_loa_V3.1_genomic.gff Lloa
## A. viteae (Blaxter)
${SCRIPTS}/generate_interproscan.sh Aviteae_GCA_900537255.1_ASM90053725v1_genomic.fna Aviteae_GCA_900537255.1_ASM90053725v1_genomic.gff Aviteae
## L. sigmodontis
${SCRIPTS}/generate_interproscan.sh litomosoides_sigmodontis.PRJEB3075.WBPS19.genomic.fa litomosoides_sigmodontis.PRJEB3075.WBPS19.annotations.gff3 Lsigmodontis
```

#### Geneinfo
The *_gene_names.txt files are available in the input_files folder. They were created by parsing the gene id from the 9th column of the GFF file for the gene features
```bash
${SCRIPTS}/generate_geneinfo.sh Dimmitis.gene_sequences.aa.tsv Dimmitis_gene_names.txt
${SCRIPTS}/generate_geneinfo.sh Lloa.gene_sequences.aa.tsv Lloa_gene_names.txt
${SCRIPTS}/generate_geneinfo.sh Lsigmodontis.gene_sequences.aa.tsv Lsigmodontis_gene_names.txt
${SCRIPTS}/generate_geneinfo.sh Aviteae.gene_sequences.aa.tsv Aviteae_blaxter_gene_names.txt
${SCRIPTS}/generate_geneinfo.sh Ovolvulus.gene_sequences.aa.tsv Ovolvulus_PRJEB513_gene_names_OVOC.txt
```
## Sequencing Depth
The sequencing depth plot is shown in Fig 4.2 and the median sequencing depths are shown in Table 4.2
### Generate Depth File
#### Inputs
```bash
## Set up directories
PACKAGE_DIR=
PROJECT_DIR=

FASTQ_DIR=${PROJECT_DIR}/fastq
BAM_DIR=${PROJECT_DIR}/final_assembly/bam
COVERAGE_DIR=${PROJECT_DIR}/final_assembly/coverage_depth/
AVITEAE_REF=${PROJECT_DIR}/reference/GCA_046563165.1_ASM4656316v1_genomic.fna

mkdir -p ${COVERAGE_DIR} ${BAM_DIR}

## Package Paths
SAMTOOLS_BIN_DIR=${PACKAGE_DIR}/samtools-1.9/bin
HISAT_DIR=${PACKAGE_DIR}/hisat2-2.1.0/
MINIMAP2_BIN_DIR=${PACKAGE_DIR}/minimap2-2.17/bin
```
#### Commands
##### Illumina
```bash
## Map DNA Illumina 
FASTQ1=$(echo ${FASTQ_DIR}/DNA_male_1_R1.fastq.gz)
FASTQ2=$(echo ${FASTQ_DIR}/DNA_male_1_R2.fastq.gz)
DEPTH_FILE=${COVERAGE_DIR}/DNA_male_1.depth.coverage 

INDEXED_REFERENCE=${PROJECT_DIR}/reference/GCA_046563165.1_ASM4656316v1_genomic

THREADS=10



## Generate Hisat Index
${HISAT_DIR}/hisat2-build ${AVITEAE_REF} ${INDEXED_REFERENCE}

## Align reads to genome, index bam file, calculate depth at each position
echo -e " 
${HISAT_DIR}/hisat2 -p ${THREADS} --no-spliced-alignment -x ${INDEXED_REFERENCE} -1 ${FASTQ1} -2 ${FASTQ2} | ${SAMTOOLS_BIN_DIR}/samtools sort -@ ${THREADS} -o ${BAM_DIR}/DNA_male_1.sorted.bam && \
    ${SAMTOOLS_BIN_DIR}/samtools index ${BAM_DIR}/DNA_male_1.sorted.bam &&\
    ${SAMTOOLS_BIN_DIR}/samtools depth -aa ${BAM_DIR}/DNA_male_1.sorted.bam > ${DEPTH_FILE}

" | qsub -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=30G -wd ${COVERAGE_DIR} -N depth_file.${SAMPLE}
```
##### CLR, HiFi, and ONT Long Reads
```bash
## Map CLR, HiFi, ONT Long Reads with minimap2
SAMPLE=DNA_male_2_CLR
MAPPING_TYPE=map-pb
DEPTH_FILE=${COVERAGE_DIR}/${SAMPLE}.depth.coverage 

SAMPLE=DNA_female_1_HiFi
MAPPING_TYPE=asm20
DEPTH_FILE=${COVERAGE_DIR}/${SAMPLE}.depth.coverage 

SAMPLE=DNA_female_ONT_min106
MAPPING_TYPE=map-ont
DEPTH_FILE=${COVERAGE_DIR}/${SAMPLE}.depth.coverage 

echo -e " 
${MINIMAP2_BIN_DIR}/minimap2 -t ${THREADS} -ax ${MAPPING_TYPE} ${AVITEAE_REF} ${FASTQ_DIR}/${SAMPLE}.fastq.gz | ${SAMTOOLS_BIN_DIR}/samtools sort -@ ${THREADS} -o ${BAM_DIR}/${SAMPLE}.sorted.bam && \
${SAMTOOLS_BIN_DIR}/samtools index ${BAM_DIR}/${SAMPLE}.sorted.bam && \
${SAMTOOLS_BIN_DIR}/samtools depth -aa ${BAM_DIR}/${SAMPLE}.sorted.bam > ${DEPTH_FILE}
" | qsub -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=30G -wd ${COVERAGE_DIR} -N minimap2.${SAMPLE}
```

### Plot Sequencing Depth
```bash
## Plot Depth File (Illumina Only)
${PACKAGE_DIR}/r-4.1.2/bin/Rscript ${SCRIPTS}/plot_sequencing_depth.R --depthfile ${DEPTH_FILE} --outdir ${BAM_DIR}
## Plot Median Sequencing Depth
${PACKAGE_DIR}/r-4.1.2/bin/Rscript ${SCRIPTS}/plot_sequencing_depth.R --depthfile ${DEPTH_FILE} --outdir ${BAM_DIR} --median_sequencing_depth
```
