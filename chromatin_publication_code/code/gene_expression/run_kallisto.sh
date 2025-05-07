#!/bin/bash -l

# index the transcriptome
KALINDEX=data/functional_annotation/gene_expression/kallisto/index

mkdir -p $KALINDEX

/home/tom/Documents/tools/kallisto_0.48.0/kallisto index \
	-i "$KALINDEX"/index \
	data/reference_genome/annotation/Ahypochondriacus_v3.spliced_exons.fasta


##### Perform quantification
# create array of read fastq files (R1 only):
SOURCE_DIR=../Ahyp_v2_2/raw_data/Clouse_short_reads/
FILES=("$SOURCE_DIR"SRR*_1.fastq.gz)
OUTDIR=data/functional_annotation/gene_expression/kallisto/
TISSUE_NAMES=("young_seed" "stem" "flower" "leaf" "root" "water_stressed" "cotyledones" "mature_seed")

# kallisto after indexing
for (( i=0; i<=7; i++))
do
	/home/tom/Documents/tools/kallisto_0.48.0/kallisto quant -i "$KALINDEX"/index -o "$OUTDIR""${TISSUE_NAMES[$i]}" --bias --plaintext -t 6 --verbose "${FILES[$i]}" "${FILES[$i]/_1.fastq.gz/_2.fastq.gz}"
done
