#!/bin/bash -l
#SBATCH -D /scratch/cgraf2/ATAC_seq
#SBATCH -o /scratch/cgraf2/ATAC_seq/logs/macs3-%j.txt
#SBATCH -e /scratch/cgraf2/ATAC_seq/logs/macs3-%j.err
#SBATCH -t 24:00:00
#SBATCH -J macs3
#SBATCH --array=0-45
#SBATCH --partition=smp
#SBATCH --mem 8g
#SBATCH --mail-type=ALL # if you want emails, otherwise remove
#SBATCH --account=ag-stetter
#SBATCH --mail-user=cgraf2@uni-koeln.de # receive an email with updates
###SBATCH --test-only

####SLURM_ARRAY_TASK_ID=0
source /scratch/cgraf2/ATAC_seq/miniconda3/bin/activate
conda activate /scratch/cgraf2/miniconda/envs/macs3_env #atac


INPUTPATH=/projects/ag-stetter/cgraf2/ATAC_seq/cheops_data/data/Acru_bam_files/depth_filtered
OUTPUTPATH=/scratch/cgraf2/ATAC_seq/data/peak_calling/macs3/Acru/df250


mkdir -p ${OUTPUTPATH}/

#BAMFILES=($(cat general_projects/data/reference_panle.txt))                           # path to all bam files
BAMFILES=($(ls -d ${INPUTPATH}/*.bam)) 

INFILE="${BAMFILES[$SLURM_ARRAY_TASK_ID]}"                       # path to current iterration bam file 
#INFILE=/scratch/cgraf2/atac/data/wgs/mapped/SRR2106212_Ahyp.bam

echo "current bam files" >&2
echo $BAMFILES >&2

#INDNAME=$(basename $INFILE _Ahypv3.bam)				#create base name (e.g. AM-00354_S1) for other filenames           
#INDNAME=$(basename $INFILE _Acru.bam)
INDNAME=$(basename $INFILE _depth_filtered250.bam)


echo " this is INFILE" >&2
echo $INFILE >&2

echo "this is INDNAME" >&2
echo $INDNAME >&2

#macs3 callpeak -f BAMPE -t ${INFILE} -g 430000000 --outdir ${OUTPUTPATH} -n ${INDNAME} -B -q 0.01
macs3 callpeak -f BAMPE -t ${INFILE} -g 370000000 --outdir ${OUTPUTPATH} -n ${INDNAME} -B -q 0.01



#chip_Seq course script

#macs2 callpeak --keep-dup all -t ${INFILE} -c control_file -n ${OUTPUTPATH}/${INDNAME} --format=BAM --gsize 500000000 --bw=300 --qvalue 0.05
# -t -> data from ChIP-Seq (bam files of replicates)
# -c -> control/mock data (inputs files)
# --format=BAM -> format of input files
# --gsize hs -> mappable genome size
# --bw=300 -> bandwith, used to scan genome for model building
# --qvalue 0.05 -> cutoff: q-value/minimum FDR
# --keep-dup all , means keep duplicats, because we removed duplicats before, here we should have only biologial duplicats

