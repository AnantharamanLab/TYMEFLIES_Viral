#!/usr/bin/perl

use strict;
use warnings;

# Aim: Conduct Fst analysis for summer and winter sample groups
# Usage of inStrain_lite: python3 /slowdata/data1/hydrothermal_plume_omics/Genome/inStrain_lite/bin/inStrain_lite bam fasta
=pod
# Step 1 Calculate the median insert size
# Note: This command should be run under BBTools_v37.62 env (conda activate BBTools_v37.62)
`reformat.sh in=Summer_vs_Winter_Fst_analysis/All_winter_IMG.viral_species_rep.id90.filtered.bam ihist=Summer_vs_Winter_Fst_analysis/All_winter_IMG.viral_species_rep.id90.filtered.ihist.txt reads=5000000`;
# The result is: Insert Size: Mean 292.178
`reformat.sh in=Summer_vs_Winter_Fst_analysis/All_summer_IMG.viral_species_rep.id90.filtered.bam ihist=Summer_vs_Winter_Fst_analysis/All_summer_IMG.viral_species_rep.id90.filtered.ihist.txt reads=5000000`;
# The result is: Insert Size: Mean 315.669
=cut
# Step 2 Filter reads by inStrain_lite
# Note: (1) Use 3 * median insert size as the cutoff for -l (Maximum insert size between two reads)
# Note: (1) This command should be run under inStrain_lite (conda activate inStrain_lite)
my $refer_fasta = "reference_fasta_for_metapop/All_phage_species_rep_gn.fasta"; # Use all species rep gn as the refer fasta
## Firstly make bam index
#`samtools index Summer_vs_Winter_Fst_analysis/All_winter_IMG.viral_species_rep.id90.filtered.bam`;
`python3 /slowdata/data1/hydrothermal_plume_omics/Genome/inStrain_lite/inStrain_lite/filter_reads.py -m 0.93 -q -1 -l 877 Summer_vs_Winter_Fst_analysis/All_winter_IMG.viral_species_rep.id90.filtered.bam $refer_fasta -g`;
#`samtools index Summer_vs_Winter_Fst_analysis/All_summer_IMG.viral_species_rep.id90.filtered.bam`;
`python3 /slowdata/data1/hydrothermal_plume_omics/Genome/inStrain_lite/inStrain_lite/filter_reads.py -m 0.93 -q -1 -l 947 Summer_vs_Winter_Fst_analysis/All_summer_IMG.viral_species_rep.id90.filtered.bam $refer_fasta -g`;

`mv All_winter_IMG_filtered_sort.bam Summer_vs_Winter_Fst_analysis`;
`samtools index Summer_vs_Winter_Fst_analysis/All_winter_IMG_filtered_sort.bam`;
`mv All_summer_IMG_filtered_sort.bam Summer_vs_Winter_Fst_analysis`;
`samtools index Summer_vs_Winter_Fst_analysis/All_summer_IMG_filtered_sort.bam`;

# Step 3 Run inStrain_lite
`inStrain_lite Summer_vs_Winter_Fst_analysis/All_winter_IMG_filtered_sort.bam $refer_fasta -p 10 -o Summer_vs_Winter_Fst_analysis/All_winter_IMG_inStrain_lite_out`;
`inStrain_lite Summer_vs_Winter_Fst_analysis/All_summer_IMG_filtered_sort.bam $refer_fasta -p 10 -o Summer_vs_Winter_Fst_analysis/All_summer_IMG_inStrain_lite_out`;

# Step 4 Run Fst
`mkdir Summer_vs_Winter_Fst_analysis/Fst_for_each_viral_gn`;
open IN, "ls Summer_vs_Winter_Fst_analysis/Genes/*.genes |";
while (<IN>){
	chomp;
	my $gene_file = $_;
	my ($viral_gn) = $gene_file =~ /Genes\/(.+?)\.genes/;
	`python3 /slowdata/data1/hydrothermal_plume_omics/Genome/inStrain_lite/inStrain_lite/fst.py Summer_vs_Winter_Fst_analysis/All_winter_IMG_inStrain_lite_out Summer_vs_Winter_Fst_analysis/All_summer_IMG_inStrain_lite_out -g $gene_file -o Summer_vs_Winter_Fst_analysis/Fst_for_each_viral_gn/summer_vs_winter.fst_result.$viral_gn`;
	
}
close IN;
=pod
# Step 5 Run inStrain
# Note: (1) This command should be run under inStrain (conda activate inStrain)
# Note: (2) "Scf2genome.stb" is the scaffold to genome mapping file manually created
`inStrain profile Summer_vs_Winter_Fst_analysis/All_winter_IMG_filtered_sort.bam Summer_vs_Winter_Fst_analysis/Viral_species_containing_four_AMGs.fasta -o Summer_vs_Winter_Fst_analysis/All_winter_IMG.IS -p 6 -g Summer_vs_Winter_Fst_analysis/Genes/All_viral_species_containing_four_AMGs.genes.fna -s Summer_vs_Winter_Fst_analysis/Scf2genome.stb`;
`inStrain profile Summer_vs_Winter_Fst_analysis/All_summer_IMG_filtered_sort.bam Summer_vs_Winter_Fst_analysis/Viral_species_containing_four_AMGs.fasta -o Summer_vs_Winter_Fst_analysis/All_summer_IMG.IS -p 6 -g Summer_vs_Winter_Fst_analysis/Genes/All_viral_species_containing_four_AMGs.genes.fna -s Summer_vs_Winter_Fst_analysis/Scf2genome.stb`;
=cut


