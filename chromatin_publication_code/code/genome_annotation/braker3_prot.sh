#!/bin/bash -l

#Braker2 protein input using the plant sequences from orthoDB (downloaded from: https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz) as well as the protein sequences from amaranthus cruentus (removed asterisks and space in fasta header)
#The downloaded dataset contains sequences from 117 embryophyte species.

#wget https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz
#tar -xvf odb10_plants_fasta.tar.gz
# write all into the same file:
#cat plants/Rawdata/* > protein_db.fasta

# also add the Cruentus sequences:
#cat /projects/ag-stetter/twinkle/Amaranthus_cruentus/Amacr_pep_nospace.fa >> protein_db.fasta
# removed asterisk
#sed 's/\*//' protein_db.fasta > protein_db.fa

#-> in total 3536219 plant protein sequences

REFERENCE="data/reference_genome/assembly/Ahypochondriacus_v3.softmasked.fasta"

# this script is used to run braker3 in protein mode

source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate /opt/rrzk/software/conda-envs/braker3


module load bedtools/2.31.0
export PATH=$PATH:/home/twinkle1/tools/gffread/
export PATH=$PATH:/home/twinkle1/tools/stringtie-2.2.1.Linux_x86_64/
module load samtools/1.18

mkdir -p data/braker/prot_braker3

# run braker in protein mode:

braker.pl --AUGUSTUS_CONFIG_PATH=/home/twinkle1/tools/config/ \
	--prot_seq=data/braker/input/protein/protein_db.fasta \
	--genome="$REFERENCE" \
	--softmasking \
	--species=assembly_v3_prot \
	--threads=8 \
	--workingdir=data/braker/prot_braker3
