# code/genome_annotation/

Generate computational genome annotation:
- braker3_prot.sh = initial braker3 run using only protein evidence
- index_STAR.sh = index reference assembly for STAR
- run_STAR.sh = map RNA-seq reads to the reference using STAR, guided by initial braker3 annotation
- braker3_prot_rna.sh = braker3 run using protein and RNA-seq evidence
- analyse_annotation.sh = assess annotation completeness using BUSCO
- finalise_and_analyse_annotation.Rmd = analyse annotation after combining computational annotation with Iso-seq transcripts. Rename transcripts to follow a consistent naming convention. Assess final annotation completeness.
