#!/bin/bash -l


# Beforehand:
#Short read data downloaded from SRA using the following accession/run numbers:

#SRA Accession numbers:
#Floral tissue: SRX722058 SRR1598911
#Leaf tissue: SRX722059 SRR1598912
#Root tissue: SRX722060 SRR1598913
#Stem tissue: SRX722057 SRR1598910
#Water stressed tissue sample: SRX722061 SRR1598914
#Immature seeds: SRX722056 SRR1598909
#Mature seeds: SRX722063 SRR1598916
#Green Cotyledone: SRX722062 SRR1598915

#using the following commands:
#Download (show progress):
#/home/twinkle1/tools/sratoolkit.2.11.2-centos_linux64/bin/prefetch -p -O Clouse_short_reads/ SRR1598916
#Converted to fastq (reads separated into two files, with an additional file for unpaired reads):
#/home/twinkle1/tools/sratoolkit.2.11.2-centos_linux64/bin/fastq-dump --split-3 --outdir /scratch/twinkle1/Clouse_short_reads/ Clouse_short_reads/SRR1598909>


module load star/2.7.8a

# create array of read fastq files (R1 only):
SOURCE_DIR=../Ahyp_v2_2/raw_data/Clouse_short_reads
FILES=("$SOURCE_DIR"/SRR*_1.fastq.gz)

# run STAR after genome index creation
mkdir -p data/STAR/STAR_mappings/

STAR --runThreadN 8 \
	--runMode alignReads \
	--outSAMtype BAM SortedByCoordinate \
	--genomeDir data/STAR/STAR_index \
	--readFilesCommand zcat \
	--outSAMstrandField intronMotif \
	--outFileNamePrefix data/STAR/STAR_mappings/SRR_"${SLURM_ARRAY_TASK_ID}"_ \
	--readFilesIn "${FILES["${SLURM_ARRAY_TASK_ID}"]}" "${FILES["${SLURM_ARRAY_TASK_ID}"]/_1.fastq.gz/_2.fastq.gz}"
