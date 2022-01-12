#!/usr/bin/perl

use strict;
use warnings;

# Aim: Using dRep to cluster phage genomes in this study and vOTU representatives

# Copy and rename all phages genomes in this study
=pod
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


# Make phage_genomes_list.txt (contains both vOTU representatives and phage genomes from this study)
`find /storage1/data11/TYMEFLIES_phage/dRep_working_dir/phage_genomes -name "*.fasta" > /storage1/data11/TYMEFLIES_phage/dRep_working_dir/phage_genomes_list.txt`;
=cut
# Make individual phage genome lists for individual years
# Store all metagenomes info
my %Meta_info = (); # $img_id => $year (for example, "2000")
my %Years = (); # $year => 1
open IN, "TYMEFLIES_metagenome_info.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		my @tmp = split (/\t/);
		my $img_id = $tmp[0];
		my $date = $tmp[8];
		my ($year) = $date =~ /^(\d\d\d\d)\-/;
		$Years{$year} = 1;
		$Meta_info{$img_id} = $year;
	}
}
close IN;

# Make each year's genome list
my %IMGID2Gn_adr = (); # Store all genome adr; $img_id => $gn_adr collection (separated by "\t")
open IN, "/storage1/data11/TYMEFLIES_phage/dRep_working_dir/phage_genomes_list.txt";
while (<IN>){
	chomp;
	my $gn_adr = $_;
	if ($gn_adr =~ /vRhyme/){
		my ($img_id) = $gn_adr =~ /phage_genomes\/(\d+?)__vRhyme/;
		if (!exists $IMGID2Gn_adr{$img_id}){
			$IMGID2Gn_adr{$img_id} = $gn_adr;
		}else{
			$IMGID2Gn_adr{$img_id} .= "\t".$gn_adr;
		}
	}
}
close IN;

foreach my $year (sort keys %Years){
	my %Gn_adrs_for_this_year = (); # $gn_list => 1
	foreach my $img_id (sort keys %IMGID2Gn_adr){
		if ($Meta_info{$img_id} eq $year){
			my $gn_adrs = $IMGID2Gn_adr{$img_id};
			my @Gn_adrs = split (/\t/,$gn_adrs);
			foreach my $gn_adr (@Gn_adrs){
				$Gn_adrs_for_this_year{$gn_adr} = 1;
			}
		}
	}
	
	open OUT, ">/storage1/data11/TYMEFLIES_phage/dRep_working_dir/phage_genomes_list.$year.txt";
	foreach my $gn_adr (sort keys %Gn_adrs_for_this_year){
		print OUT "$gn_adr\n";
	}
	close OUT;
}

# Run dRep 
open OUT, ">tmp.run_dRep.sh";
foreach my $year (sort keys %Years){
	print OUT "dRep dereplicate dRep_working_dir/Output_directory.$year -p 15 -g /storage1/data11/TYMEFLIES_phage/dRep_working_dir/phage_genomes_list.$year.txt -l 2000 --ignoreGenomeQuality -pa 0.8 -sa 0.95 -nc 0.85 --multiround_primary_clustering --run_tertiary_clustering -comW 0 -conW 0 -strW 0 -N50W 0 -sizeW 1 -centW 0\n";
}
close OUT;

`cat tmp.run_dRep.sh | parallel -j 2`;

`rm tmp.run_dRep.sh`;


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



