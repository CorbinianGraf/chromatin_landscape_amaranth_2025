#!/bin/bash -l

# analyse the computational genome annotation
# calculate busco score for braker run using only protein evidence
busco -m protein -i data/braker/prot_braker3/braker.aa -o busco_braker3_prot -l embryophyta_odb10 --out_path data/braker/prot_braker3 --download_path ../Ahyp_v2_2/data/annotation_analysis/busco/datasets/ -c 6

# fix braker gtf file:
/home/tom/Documents/tools/gffread-0.12.7/gffread -T data/braker/prot_braker3/braker.gtf > data/braker/prot_braker3/fixed_braker.gtf


# calculate busco score for braker run using protein and rna evidence
busco -m protein -i data/braker/prot_rna_braker3/braker.aa -o busco_braker3_prot_rna -l embryophyta_odb10 --out_path data/braker/prot_rna_braker3 --download_path ../Ahyp_v2_2/data/annotation_analysis/busco/datasets/ -c 6

# fix braker gtf file:
/home/tom/Documents/tools/gffread-0.12.7/gffread -T data/braker/prot_rna_braker3/braker.gtf > data/braker/prot_rna_braker3/fixed_braker.gtf
