#!/bin/bash -l
#SBATCH -D /scratch/cgraf2/atac
#SBATCH -o /scratch/cgraf2/atac/logs/index_referenceg_methyl-%j.txt
#SBATCH -e /scratch/cgraf2/atac/logs/index_referenceg_methyl-%j.err
#SBATCH -t 12:00:00
#SBATCH -J index_methylation
#SBATCH --array=0-0
#SBATCH --nodes=1-1
#SBATCH --ntasks 8
#SBATCH --mem 8g
