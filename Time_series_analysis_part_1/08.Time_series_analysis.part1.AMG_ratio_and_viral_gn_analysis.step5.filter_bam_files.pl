#!/usr/bin/perl

use strict;
use warnings;

# Aim: Move all *.viral_species_rep.id90.bam files into a new folder and 
# filter bams by scaffolds of viral species representatives

# Note: Need to first enter vRhyme conda env by "conda activate vRhyme"

# Step 1 Store the address of all bam files
my %Bam_files = (); # $bam_file => 1
open IN, "find /storage1/data11/TYMEFLIES_phage/33* -maxdepth 2 -name '*.viral_species_rep.id90.bam'| ";
while (<IN>){
	chomp;
	my $bam_file = $_;
	$Bam_files{$bam_file} = 1;
}
close IN;

# Step 2 Move bam files
`mkdir viral_species_rep_bams`;
`mkdir viral_species_rep_bams/original`;

foreach my $bam_file (sort keys %Bam_files){
	`mv $bam_file viral_species_rep_bams/original`;
}

# Step 3 File bams
## The scaffold name of viral species representatives, which will be used as the reference for bam filtering
my $scaf_name_file = "/storage1/data11/TYMEFLIES_phage/viral_species_representative_scaf_name.txt";

`/slowdata/scripts/python_scripts/filter_bam_by_reference.py -b viral_species_rep_bams/original/*.bam -r $scaf_name_file -j 20`;



