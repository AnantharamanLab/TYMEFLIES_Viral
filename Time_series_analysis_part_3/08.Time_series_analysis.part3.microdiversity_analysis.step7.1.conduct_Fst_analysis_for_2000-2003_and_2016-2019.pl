#!/usr/bin/perl

use strict;
use warnings;

# Aim: Conduct Fst analysis for year 2000-2003 and year 2016-2019
# Usage of inStrain_lite: python3 /slowdata/data1/hydrothermal_plume_omics/Genome/inStrain_lite/bin/inStrain_lite bam fasta

# Step 1 Merge bam files for year 2000-2003 and 2016-2019
# Note: This command should be run under inStrain_lite (conda activate inStrain_lite)
`samtools merge viral_species_rep_bams.for_each_year/2000-2003.viral_species_rep.id90.bam viral_species_rep_bams.for_each_year/2000.viral_species_rep.id90.bam viral_species_rep_bams.for_each_year/2001.viral_species_rep.id90.bam viral_species_rep_bams.for_each_year/2002.viral_species_rep.id90.bam viral_species_rep_bams.for_each_year/2003.viral_species_rep.id90.bam`;
`samtools merge viral_species_rep_bams.for_each_year/2016-2019.viral_species_rep.id90.bam viral_species_rep_bams.for_each_year/2016.viral_species_rep.id90.bam viral_species_rep_bams.for_each_year/2017.viral_species_rep.id90.bam viral_species_rep_bams.for_each_year/2018.viral_species_rep.id90.bam viral_species_rep_bams.for_each_year/2019.viral_species_rep.id90.bam`;

# Step 2 Calculate the median insert size
# Note: This command should be run under BBTools_v37.62 env (conda activate BBTools_v37.62)
`reformat.sh in=viral_species_rep_bams.for_each_year/2000-2003.viral_species_rep.id90.bam ihist=viral_species_rep_bams.for_each_year/2000-2003.viral_species_rep.id90.ihist.txt reads=5000000`;
# The result is: Insert Size: Mean 314.066
`reformat.sh in=viral_species_rep_bams.for_each_year/2016-2019.viral_species_rep.id90.bam ihist=viral_species_rep_bams.for_each_year/2016-2019.viral_species_rep.id90.ihist.txt reads=5000000`;
# The result is: Insert Size: Mean 318.678

# Step 3 Filter reads by inStrain_list
# Note: (1) Use 3 * median insert size as the cutoff for -l (Maximum insert size between two reads)
# Note: (2) This command should be run under inStrain_lite (conda activate inStrain_lite)
my $refer_fasta = "viral_species_rep_bams.for_each_year/Viral_species_containing_four_AMGs.fasta";
## Firstly make bam index
`samtools index viral_species_rep_bams.for_each_year/2000-2003.viral_species_rep.id90.bam`;
`samtools index viral_species_rep_bams.for_each_year/2016-2019.viral_species_rep.id90.bam`;
`python3 /slowdata/data1/hydrothermal_plume_omics/Genome/inStrain_lite/inStrain_lite/filter_reads.py -m 0.95 -q -1 -l 942 viral_species_rep_bams.for_each_year/2000-2003.viral_species_rep.id90.bam $refer_fasta -g`;
`python3 /slowdata/data1/hydrothermal_plume_omics/Genome/inStrain_lite/inStrain_lite/filter_reads.py -m 0.95 -q -1 -l 957 viral_species_rep_bams.for_each_year/2016-2019.viral_species_rep.id90.bam $refer_fasta -g`;

`mv 2000-2003_filtered_sort.bam viral_species_rep_bams.for_each_year`;
`samtools index viral_species_rep_bams.for_each_year/2000-2003_filtered_sort.bam`;
`mv 2016-2019_filtered_sort.bam viral_species_rep_bams.for_each_year`;
`samtools index viral_species_rep_bams.for_each_year/2016-2019_filtered_sort.bam`;

# Step 4 Run inStrain_lite
`inStrain_lite viral_species_rep_bams.for_each_year/2000-2003_filtered_sort.bam $refer_fasta -p 10 -o viral_species_rep_bams.for_each_year/2000-2003_inStrain_lite_out`;
`inStrain_lite viral_species_rep_bams.for_each_year/2016-2019_filtered_sort.bam $refer_fasta -p 10 -o viral_species_rep_bams.for_each_year/2016-2019_inStrain_lite_out`;

# Step 5 Run Fst
`mkdir viral_species_rep_bams.for_each_year/Fst_for_each_viral_gn`;
open IN, "ls viral_species_rep_bams.for_each_year/Genes/*.genes |";
while (<IN>){
	chomp;
	my $gene_file = $_;
	my ($viral_gn) = $gene_file =~ /Genes\/(.+?)\.genes/;
	`python3 /slowdata/data1/hydrothermal_plume_omics/Genome/inStrain_lite/inStrain_lite/fst.py viral_species_rep_bams.for_each_year/2000-2003_inStrain_lite_out viral_species_rep_bams.for_each_year/2016-2019_inStrain_lite_out -g $gene_file -o viral_species_rep_bams.for_each_year/Fst_for_each_viral_gn/2000-2003_vs_2016-2019.fst_result.$viral_gn`;	
}
close IN;

# Step 6 Run inStrain
# Note: (1) This command should be run under inStrain (conda activate inStrain)
# Note: (2) "Scf2genome.stb" is the scaffold to genome mapping file manually created
`inStrain profile viral_species_rep_bams.for_each_year/2000-2003_filtered_sort.bam viral_species_rep_bams.for_each_year/Viral_species_containing_four_AMGs.fasta -o viral_species_rep_bams.for_each_year/2000-2003.IS -p 6 -g viral_species_rep_bams.for_each_year/Genes/All_viral_species_containing_four_AMGs.genes.fna -s viral_species_rep_bams.for_each_year/Scf2genome.stb`;
`inStrain profile viral_species_rep_bams.for_each_year/2016-2019_filtered_sort.bam viral_species_rep_bams.for_each_year/Viral_species_containing_four_AMGs.fasta -o viral_species_rep_bams.for_each_year/2016-2019.IS -p 6 -g viral_species_rep_bams.for_each_year/Genes/All_viral_species_containing_four_AMGs.genes.fna -s viral_species_rep_bams.for_each_year/Scf2genome.stb`;