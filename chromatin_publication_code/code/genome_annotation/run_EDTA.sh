#!/bin/bash -l

# this script runs EDTA as de nove TE annotation pipeline on the genome assembly v3

# load requirements
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate /scratch/twinkle1/EDTA2_env

module load samtools/1.18

# variables:
OUTDIR=/scratch/twinkle1/EDTA2_out/

# outdir
mkdir -p "$OUTDIR"

# capitalize input fasta
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' data/reference_genome/assembly/Ahypochondriacus_v3.softmasked.fasta > "$OUTDIR"Ahypochondriacus_v3.fasta

# EDTA outdir is wd
cd "$OUTDIR"

# create and set temporary directory
mkdir -p /scratch/twinkle1/temp_dir/
export TMPDIR=/scratch/twinkle1/temp_dir

# tmpdir dry run:
mktemp -u

# run EDTA
/home/twinkle1/tools/EDTA/EDTA.pl --genome Ahypochondriacus_v3.fasta \
	--step all \
	--cds /projects/ag-stetter/twinkle/projects/annotation_Ahyp_v3/data/reference_genome/annotation/Ahypochondriacus_v3.cds.fasta \
	--sensitive 1 \
	--anno 1 \
	--threads 12
