#!/usr/bin/perl

use strict;
use warnings;

# Aim: Using dRep to cluster phage genomes in this study and vOTU representatives

# Copy and rename all phages genomes in this study

`mkdir dRep_working_dir`;
`mkdir dRep_working_dir/phage_genomes`;

open OUT, ">tmp.copy_and_rename_all_phage_genomes.sh";
open IN, "find */vRhyme_best_bins_fasta_parsed/ -maxdepth 1  -name '*vRhyme*.fasta' | ";
while (<IN>){
	chomp;
	my $file = $_;
	my ($img_id, $file_name) = $file =~ /^(\d+?)\/vRhyme_best_bins.+?\/(33.+?\_\_vRhyme\_.+?)\.fasta/;
	print OUT "cp $file dRep_working_dir/phage_genomes/$file_name.fasta\n";
}
close IN;
close OUT;

`cat tmp.copy_and_rename_all_phage_genomes.sh | parallel -j 20`;

`rm tmp.copy_and_rename_all_phage_genomes.sh`;


# Copy vOTU representatives phage genomes to "dRep_working_dir/phage_genomes"

open OUT, ">tmp.copy_and_rename_vOTU_representatives_phage_genomes.sh";
open IN, "find /slowdata/databases/IMGVR-NCBI_phages/vOTU_representatives_phage_genomes -name '*.fasta' | ";
while (<IN>){
	chomp;
	my $file = $_;
	print OUT "cp $file dRep_working_dir/phage_genomes\n";
}
close IN;
close OUT;

`cat tmp.copy_and_rename_vOTU_representatives_phage_genomes.sh | parallel -j 20`;

`rm tmp.copy_and_rename_vOTU_representatives_phage_genomes.sh`;


# Run dRep 
`dRep dereplicate dRep_working_dir/Output_directory -p 50 -g /storage1/data11/TYMEFLIES_phage/dRep_working_dir/phage_genomes_list.txt -l 2000 --ignoreGenomeQuality -pa 0.8 -sa 0.95 -nc 0.85 --multiround_primary_clustering --run_tertiary_clustering -comW 0 -conW 0 -strW 0 -N50W 0 -sizeW 1 -centW 0`;

##########################################################################
# Details of options (from website https://drep.readthedocs.io/en/latest/module_descriptions.html#dereplicate):
# -l
#                       Minimum genome length (default: 50000)   
#
#                       >>> here set to 2000
# --ignoreGenomeQuality
#                       Don't run checkM or do any quality filtering. NOT
#                       RECOMMENDED! This is useful for use with
#                       bacteriophages or eukaryotes or things where checkM
#                       scoring does not work. Will only choose genomes based
#                       on length and N50 (default: False)
#
#                       >>> here set to True
# -pa
#                       ANI threshold to form primary (MASH) clusters
#                       (default: 0.9)
#
#                       >>> here set to 0.8
# -sa
#                       ANI threshold to form secondary clusters (default:
#                       (default: 0.99)
#
#                       >>> here set to 0.95
# -nc
#                       Minmum level of overlap between genomes when doing
#                       secondary comparisons (default: 0.1)
#
#                       >>> here set to 0.85
# GREEDY CLUSTERING OPTIONS
# These decrease RAM use and runtime at the expense of a minor loss in accuracy.
# Recommended when clustering 5000+ genomes:
# --multiround_primary_clustering
#                       Cluster each primary clunk separately and merge at the
#                       end with single linkage. Decreases RAM usage and
#                       increases speed, and the cost of a minor loss in
#                       precision and the inability to plot
#                       primary_clustering_dendrograms. Especially helpful
#                       when clustering 5000+ genomes. Will be done with
#                       single linkage clustering (default: False)
#
#                       >>> here set to True
# --primary_chunksize PRIMARY_CHUNKSIZE
#                       Impacts multiround_primary_clustering. If you have
#                       more than this many genomes, process them in chunks of
#                       this size. (default: 5000)
# --greedy_secondary_clustering
#                       Use a heuristic to avoid pair-wise comparisons when
#                       doing secondary clustering. Will be done with single
#                       linkage clustering. Only works for fastANI S_algorithm
#                       option at the moment (default: False)
# --run_tertiary_clustering
#                       Run an additional round of clustering on the final
#                       genome set. This is especially useful when greedy
#                       clustering is performed and/or to handle cases where
#                       similar genomes end up in different primary clusters.
#                       Only works with dereplicate, not compare. (default:
#                       False)
#
#                       >>> here set to True
# SCORING CRITERIA
# Based off of the formula:
# A*Completeness - B*Contamination + C*(Contamination * (strain_heterogeneity/100)) + D*log(N50) + E*log(size) + F*(centrality - S_ani)
#
# A = completeness_weight; B = contamination_weight; C = strain_heterogeneity_weight; D = N50_weight; E = size_weight; F = cent_weight:
# -comW COMPLETENESS_WEIGHT, --completeness_weight COMPLETENESS_WEIGHT
#                       completeness weight (default: 1)
#
#                       >>> here set to 0
# -conW CONTAMINATION_WEIGHT, --contamination_weight CONTAMINATION_WEIGHT
#                       contamination weight (default: 5)
#
#                       >>> here set to 0
# -strW STRAIN_HETEROGENEITY_WEIGHT, --strain_heterogeneity_weight STRAIN_HETEROGENEITY_WEIGHT
#                       strain heterogeneity weight (default: 1)
#
#                       >>> here set to 0
# -N50W N50_WEIGHT, --N50_weight N50_WEIGHT
#                       weight of log(genome N50) (default: 0.5)
#
#                       >>> here set to 0
# -sizeW SIZE_WEIGHT, --size_weight SIZE_WEIGHT
#                       weight of log(genome size) (default: 0)
#
#                       >>> here set to 1
# -centW CENTRALITY_WEIGHT, --centrality_weight CENTRALITY_WEIGHT
#                       Weight of (centrality - S_ani) (default: 1)
#
#                       >>> here set to 0



