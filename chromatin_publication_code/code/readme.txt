# code/

genome_analysis: Ensure assembly co-linearity to the previous reference genomes through reverse complementation, if necessary and plot the genomic structure as circos plot.
assemble_isoseq: Combine Iso-Seq reads from different sources, assemble into unique full-length transcripts and analyse.
genome_annotation: Annotate protein-coding genes in the assembly using BRAKER3 and Iso-Seq transcript models. Assess annotation quality using BUSCO. Annotate transposable elements using the EDTA pipeline.
gene_expression: Quantify gene expression in different tissues using kallisto.
