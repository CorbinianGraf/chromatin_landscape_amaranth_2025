
# required anaCogent conda environment
#conda activate anaCogent

# create output directory
MAPPINGOUT=data/isoseq/v2/
mkdir -p "$MAPPINGOUT"


### MAPPING

REFERENCE=data/reference_genome/assembly/Ahypochondriacus_v3.softmasked.fasta

# Align sequences to reference genome, minimap v2.26
INPUT=raw_data/Ahyp_v2_2/combined.merged_clustered.hq.fasta
OUTPUTMM2="$MAPPINGOUT"combined_aln.sam

/home/tom/Documents/tools/minimap2/minimap2 -t 6 -ax splice:hq -uf --secondary=no -a $REFERENCE $INPUT -o $OUTPUTMM2


# Before collapsing isoforms, sequences have to be sorted
OUTPUTSORT="$MAPPINGOUT"combined_aln_sorted.sam

# remove unmapped sequences from file and sort
sort -k 3,3 -k 4,4n $OUTPUTMM2 > $OUTPUTSORT



### COLLAPSE


# Collapsing isoforms

INPUTSORTED="$MAPPINGOUT"combined_aln_sorted.sam
OUTPUTCOLLAPSE="$MAPPINGOUT"combined

collapse_isoforms_by_sam.py --input $INPUT -s $INPUTSORTED -o $OUTPUTCOLLAPSE -c 0.95 -i 0.9 --max_3_diff 1000


# leave out count based filtering for now
# Cupcake support scripts after collapse
# First obtain associated count information

INPUTABUNDANCE="$MAPPINGOUT"combined.collapsed
CLUSTERREPORT=raw_data/Ahyp_v2_2/combined.merged_clustered.cluster_report.csv

get_abundance_post_collapse.py $INPUTABUNDANCE $CLUSTERREPORT


# add minimum read count, filtering only after the merge

filter_by_count.py --min_count 2 --dun_use_group_count $INPUTABUNDANCE


# filter away 5' degraded isoforms, not working without abundance information

OUTPUTFILTERED="$MAPPINGOUT"combined.collapsed.min_fl_2

filter_away_subset.py $OUTPUTFILTERED

# statistics after collapse
simple_stats_post_collapse.py "$OUTPUTFILTERED".filtered
#simple_stats_post_collapse.py "$INPUTABUNDANCE"

#mkdir -p data/isoseq/gmst

# clean transcript names
#sed 's/|.*//' data/isoseq/mapping_and_collapse/combined.collapsed.rep.fa > data/isoseq/mapping_and_collapse/combined.collapsed.renamed.fasta

# predict ORF using genemark-ST
#/home/tom/Documents/tools/gmst/gmst.pl --strand direct data/isoseq/mapping_and_collapse/combined.collapsed.renamed.fasta --output data/isoseq/gmst/gmst_collapsed.gff --format GFF --fnn --faa
# move log files
#mv GeneMark_hmm.mod data/isoseq/gmst/
#mv gms.log data/isoseq/gmst/

# calculate busco score for collapsed transcript set, protein mode
# conda activate busco
# busco -m protein -i data/isoseq/gmst/gmst_collapsed.gff.faa -o busco_gmst_collapsed -l embryophyta_odb10 --out_path data/isoseq/ --download_path ../Ahyp_v2_2/data/annotation_analysis/busco/datasets/ -c 6
