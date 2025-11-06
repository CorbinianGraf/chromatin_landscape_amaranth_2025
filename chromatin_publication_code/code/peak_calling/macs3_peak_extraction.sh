#!/bin/bash -l
#SBATCH -D /scratch/cgraf2/atac
#SBATCH -o /scratch/cgraf2/atac/logs/macs3_extLog-%j.txt
#SBATCH -e /scratch/cgraf2/atac/logs/macs3_extLog-%j.err
#SBATCH -t 00:10:00
#SBATCH -J macs3_extraction
#SBATCH --array=0-45
#SBATCH --nodes=1-1
#SBATCH --ntasks 8
#SBATCH --mem 8g

####SLURM_ARRAY_TASK_ID=0



INPUTPATH=data/macs3/Acru/depth_filtered                                                # sequences to map; needs application to be strated from dir containing d$
OUTPUTPATH=data/macs3/Acru/depth_filtered/extraction                                               # where to write bam files


mkdir -p ${OUTPUTPATH}/

BAMFILES=($(ls -d $INPUTPATH/*.xls))                           # path to all bam files


INFILE="${BAMFILES[$SLURM_ARRAY_TASK_ID]}"                       # path to current iterration bam file


INDNAME=$(basename $INFILE .xls)


OUTFILE=${OUTPUTPATH}/${INDNAME}


OUTFILE2=data/macs3/Acru/depth_filtered/extraction/summary
mkdir -p ${OUTFILE2}

# get only important info (chr, start, stop)
cut -d$'\t' -f1,4,5 ${INFILE} > ${OUTFILE}.csv

# get only peaks on scaffolds
grep "Scaffold" ${OUTFILE}.csv > ${OUTFILE}_nocont.csv
#grep "Contig" ${OUTFILE}.csv > ${OUTFILE}_cont.csv

# get total amount of open chrom per sample
#cat << "${INDNAME}" >> ${OUTFILE2}/lengthsum.csv
#cut -d$'\t' -f4 ${OUTFILE}_nocont.csv
awk '{s+=$2}END{print s}' ${OUTFILE}_nocont.csv >> ${OUTFILE2}/lengthsum.csv
echo "${INDNAME}" >> ${OUTFILE2}/lengthsum.csv

# get total amount of peaks per sample
grep -c "Scaffold" ${OUTFILE}_nocont.csv >> ${OUTFILE2}/peaktotal.csv
echo "${INDNAME}" >> ${OUTFILE2}/peaktotal.csv
