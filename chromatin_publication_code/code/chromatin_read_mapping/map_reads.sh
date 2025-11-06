#!/bin/bash -l
#SBATCH -D /scratch/cgraf2/ATAC_seq
#SBATCH -o /scratch/cgraf2/ATAC_seq/logs/mapping_atac-%j.txt
#SBATCH -e /scratch/cgraf2/ATAC_seq/logs/mapping_atac-%j.err
#SBATCH -t 24:00:00
#SBATCH -J map_atac
#SBATCH --array=0-0
#SBATCH --partition=smp
#SBATCH --mem 48g
#SBATCH --mail-type=ALL # if you want emails, otherwise remove
#SBATCH --account=ag-stetter
#SBATCH --mail-user=cgraf2@uni-koeln.de # receive an email with updates

####SLURM_ARRAY_TASK_ID=0



#module load bwamem2/2.2.1
#module load samtools/1.13
#module load openjdk/1.8.0_60

module load bio/bwa-mem2/2.2.1-intel-compilers-2023.1.0
module load bio/SAMtools/1.18-GCC-12.3.0

REFERENCE=/projects/ag-stetter/reference_genomes/Acruentus/V1_0/Amacr_genome.fasta
#REFERENCE=/projects/ag-stetter/cgraf2/ATAC_seq/cheops_data/data/reference_genomes/v3/Ahypochondriacus_v3.softmasked.fasta
#bwa-mem2 index $REFERENCE						# only need to do this once per system/person

PROVIDER=CCG

#INPUTPATH=/projects/ag-stetter/cgraf/atac_seq/data/raw_data/
#OUTPUTPATH=/scratch/cgraf2/atac/data/bam_files
#INPUTPATH=/projects/ag-stetter/amaranth_mapping_sequencing/mapping_F2_AM_01200_hyp/X205SC23081062-Z01-F002/01.RawData/AM_00208						# sequences to map; needs application to be strated from dir containing data
#OUTPUTPATH=/scratch/cgraf2/recomb/data/hyp/mapped/AM_00208 						# where to write bam files
#INPUTPATH=/projects/ag-stetter/cgraf2/ATAC_seq/cheops_data/data/raw_data
#OUTPUTPATH=/projects/ag-stetter/cgraf2/ATAC_seq/cheops_data/data/bam_files/v3
INPUTPATH=/projects/ag-stetter/cgraf2/ATAC_seq/cheops_data/data/wgs_data/reference_panel/sorted/fastq
OUTPUTPATH=/projects/ag-stetter/cgraf2/ATAC_seq/cheops_data/data/wgs_data/reference_panel/Acruv1_mapped

mkdir -p ${OUTPUTPATH}/metrics/

#FASTQFILESR1=($(ls -d $INPUTPATH/*1.fq.gz))				# path to R1 sequence files
#FASTQFILESR2=($(ls -d $INPUTPATH/*2.fq.gz))				# path to R2 sequence files
#FASTQFILESR1=($(ls -d $INPUTPATH/*1.fastq))                             # path to R1 sequence files
#FASTQFILESR2=($(ls -d $INPUTPATH/*2.fastq))
FASTQFILESR1=/projects/ag-stetter/cgraf2/ATAC_seq/cheops_data/data/wgs_data/reference_panel/sorted/fastq/PI642741.final_sorted.bam_1.fq
FASTQFILESR1=/projects/ag-stetter/cgraf2/ATAC_seq/cheops_data/data/wgs_data/reference_panel/sorted/fastq/PI642741.final_sorted.bam_2.fq

INFILE_R1="${FASTQFILESR1[$SLURM_ARRAY_TASK_ID]}"			# path to input file for R1
INFILE_R2="${FASTQFILESR2[$SLURM_ARRAY_TASK_ID]}"			# path to input file for R2
#INFILE_R1=/scratch/cgraf2/atac/data/wgs/fastq/SRR2106212_sorted_1.fq
#INFILE_R2=/scratch/cgraf2/atac/data/wgs/fastq/SRR2106212_sorted_2.fq

echo "there should be something below"
echo $INFILE_R1
echo $INFILE_R2
echo "there should be something above"

#INDNAME=$(basename $INFILE_R1 _R1.fq.gz)
INDNAME=$(basename $INFILE_R1 .final_sorted.bam_1.fq) 				# smiplify sample name for later use in naming files; removes _R1.fq.gz
echo $INDNAME " this should start with INDNAME"

echo maping reads of $INDNAME
TMPFOLDER=/scratch/cgraf2/ATAC_seq/tmp					# write sorted_bam files to scratch and other large files that will only be needed once
mkdir -p $TMPFOLDER							# create tmp folder on scratch, in case it has been cleared since last run

SORTED_NAME=${TMPFOLDER}/sortedV22_Acruv1_${INDNAME}.bam			# output path and name for sorted bam files

echo $SORTED_NAME


bwa-mem2 mem -t 8 -R '@RG\tID:'${INDNAME}'\tSM:'${INDNAME}'\tCN:'${PROVIDER}'\tPL:illumina' $REFERENCE $INFILE_R1 $INFILE_R2 | samtools sort -O bam -o ${SORTED_NAME}

echo mark duplicates
DEDUP_NAME=${OUTPUTPATH}/${INDNAME}_Acruv1.bam					# where to write and store depulicated bam files
METRICS_FILE=${OUTPUTPATH}/metrics/${INDNAME}_picard_Acru.txt		# where to write and store metric files of picard

mkdir -p $TMPFOLDER/picard   #might be necessary for picard to not run out of space

java -Xmx40g -jar /home/cgraf2/tools/picard.jar MarkDuplicates INPUT=${SORTED_NAME} OUTPUT=${DEDUP_NAME} METRICS_FILE=${METRICS_FILE} TMP_DIR=${TMPFOLDER}/picard 

echo index deduplicated
samtools index $DEDUP_NAME

echo calculate samtools flagstat
samtools flagstat $DEDUP_NAME > ${OUTPUTPATH}/metrics/${INDNAME}.flagstat

echo removing sorted bam
#rm $SORTED_NAME


#mkdir -p ${OUTPUTPATH}/gvcf
#GVCFFILE=${OUTPUTPATH}/gvcf/${INDNAME}_V22.g.vcf

#/projects/mstette2/tools/gatk-4.1.7.0/gatk --java-options "-Xmx48G" CreateSequenceDictionary -R $REFERENCE

#/projects/mstette2/tools/gatk-4.1.7.0/gatk --java-options "-Xmx48G" HaplotypeCaller \
#-R $REFERENCE \
#-I $DEDUP_NAME \
#-ERC GVCF \
#--output ${GVCFFILE}
