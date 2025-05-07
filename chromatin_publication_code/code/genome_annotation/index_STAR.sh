#!/bin/bash -l

module load star/2.7.8a

# Index the reference genome
# only run once per reference genome

# 8 threads, genome generation mode
# sjdbOverhang and genomeSAindexNbases settings specific for the genome
# more specific settings: use the polished, softmasked reference assembly
# as SJDB file, use the newly generated braker3 protein gtf file

REFERENCE="data/reference_genome/assembly/Ahypochondriacus_v3.softmasked.fasta"

mkdir -p data/STAR/STAR_index/

STAR --runThreadN 8 \
	--runMode genomeGenerate \
	--genomeDir data/STAR/STAR_index/ \
	--sjdbOverhang 89 \
	--genomeSAindexNbases 13 \
	--genomeFastaFiles "$REFERENCE" \
	--sjdbGTFfile data/braker/prot_braker3/braker.gtf
