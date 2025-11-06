#!/bin/bash -l
#SBATCH -D /scratch/cgraf2/methylation/code/bismarck/working_dir
#SBATCH -o /scratch/cgraf2/atac/logs/alig_methylation-%j.txt
#SBATCH -e /scratch/cgraf2/atac/logs/alig_methylation-%j.err
#SBATCH -t 120:00:00
#SBATCH -J map_methyl
#SBATCH --array=0-0
#SBATCH --nodes=1-1
#SBATCH --ntasks 8
#SBATCH --mem 48g
#SBATCH --mail-type=ALL # if you want emails, otherwise remove
#SBATCH --account=UniKoeln
#SBATCH --mail-user=cgraf2@uni-koeln.de # receive an email with updates

module load samtools
module load bowtie2

REFERENCE=/scratch/cgraf2/ref_genome/methylation/

INPUTPATH=/scratch/cgraf2/methylation/data
OUTPUTPATH=/scratch/cgraf2/methylation/data/run2
mkdir -p ${OUTPUTPATH}

#perl bismark --bowtie2 -n 1 -l 20 methylation/ -1 Plainsman_FKDN230460809-1A_ALL_L2_1_cut.fastq -2 Plainsman_FKDN230460809-1A_ALL_L2_2_cut.fastq -o ${OUTPUTPATH}


INPUTPATH=/scratch/cgraf2/methylation/data/run1
OUTPUTPATH=/scratch/cgraf2/methylation/data/run1/extraction3
mkdir -p ${OUTPUTPATH}

perl bismark_methylation_extractor -p --comprehensive ${INPUTPATH}/Plainsman_nickname.bam --output_dir ${OUTPUTPATH} --bedGraph

#perl bismark2bedGraph --dir ${OUTPUTPATH} -o Plainsman_nickname.bed ${OUTPUTPATH}/CHG_context_Plainsman_nickname.txt ${OUTPUTPATH}/CpG_context_Plainsman_nickname.txt ${OUTPUTPATH}/CHH_context_Plainsman_nickname.txt Plainsman_nickname_splitting_report.txt
