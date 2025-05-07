#!/bin/bash -l


# this script is used to run braker3 in both RNAseq and protein mode

# conda env
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate /opt/rrzk/software/conda-envs/braker3

# load modules
module load bedtools/2.31.0
export PATH=$PATH:/home/twinkle1/tools/gffread/
export PATH=$PATH:/home/twinkle1/tools/stringtie-2.2.1.Linux_x86_64/

module load samtools

# it is recommended to merge the separate bam files beforehand, as specifying many files can cause issues with braker
mkdir -p data/braker/input/rnaseq

samtools merge --threads 7 data/braker/input/rnaseq/clouse_reads_merged.bam \
	data/STAR/STAR_mappings/SRR_0_Aligned.sortedByCoord.out.bam \
	data/STAR/STAR_mappings/SRR_1_Aligned.sortedByCoord.out.bam \
	data/STAR/STAR_mappings/SRR_2_Aligned.sortedByCoord.out.bam \
	data/STAR/STAR_mappings/SRR_3_Aligned.sortedByCoord.out.bam \
	data/STAR/STAR_mappings/SRR_4_Aligned.sortedByCoord.out.bam \
	data/STAR/STAR_mappings/SRR_5_Aligned.sortedByCoord.out.bam \
	data/STAR/STAR_mappings/SRR_6_Aligned.sortedByCoord.out.bam \
	data/STAR/STAR_mappings/SRR_7_Aligned.sortedByCoord.out.bam

# it is unnecessary, however, to filter the bam files beforehand
# see (https://github.com/Gaius-Augustus/BRAKER/issues/241)

REFERENCE="data/reference_genome/assembly/Ahypochondriacus_v3.softmasked.fasta"

mkdir -p data/braker/prot_rna_braker3

# run braker in RNAseq and protein mode:
braker.pl --AUGUSTUS_CONFIG_PATH=/home/twinkle1/tools/config/ \
	--prot_seq=data/braker/input/protein/protein_db.fasta \
	--genome="$REFERENCE" \
	--softmasking \
	--species=assembly_v3_prot_rna \
	--bam=data/braker/input/rnaseq/clouse_reads_merged.bam \
	--threads=8 \
	--workingdir=data/braker/prot_rna_braker3
