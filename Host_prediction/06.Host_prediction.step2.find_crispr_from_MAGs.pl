#!/usr/bin/perl

use strict;
use warnings;

# AIM: Find crispr spacers from MAGs (both TYMEFLIES and GEM MAGs)
# NOTE: To run MinCED, it should conda activate it first before run this script: conda activate MinCED
# NOTE: To run PILER-CR, it should conda activate it first before run this script: conda activate PILER-CR_v1.06

# Step 1 Write down and run crispr searching command for TYMEFLIES MAGs and GEM MAGs
`mkdir Host_prediction/find_crispr_out`;

## Step 1.1 Write down crispr searching command by MinCED for TYMEFLIES MAGs and GEM MAGs
open OUT, ">tmp.find_crispr_by_minced_1.sh";
open IN, "find /storage1/data11/TYMEFLIES_phage/Binning_Data/33* -name '*.fna' |";
while (<IN>){
	chomp;
	my $fna_file = $_;
	my ($out_name) = $fna_file =~ /^.+\/(.+?)\.fna/;
	$out_name = "TYMEFLIES_MAGs_".$out_name;
	print OUT "minced $fna_file Host_prediction/find_crispr_out/$out_name.minced_out.crisprs\n";
}
close IN;
close OUT;

open OUT, ">tmp.find_crispr_by_minced_2.sh";
open IN, "find /slowdata/databases/GEM/fna -name '*.fna'|";
while (<IN>){
	chomp;
	my $fna_file = $_;
	my ($out_name) = $fna_file =~ /^.+\/(.+?)\.fna/;
	$out_name = "GEM_MAGs_".$out_name;
	print OUT "minced $fna_file Host_prediction/find_crispr_out/$out_name.minced_out.crisprs\n";
}
close IN;
close OUT;

# Step 1.2 Run the crispr searching command by MinCED with 10 cpus
`cat tmp.find_crispr_by_minced_1.sh | parallel -j 10`;
`cat tmp.find_crispr_by_minced_2.sh | parallel -j 10`;

`rm tmp.find_crispr_by_minced_1.sh`;
`rm tmp.find_crispr_by_minced_2.sh`;


## Step 1.3 Write down crispr searching command by PILER-CR for TYMEFLIES MAGs and GEM MAGs
open OUT, ">tmp.find_crispr_by_pilercr_1.sh";
open IN, "find /storage1/data11/TYMEFLIES_phage/Binning_Data/33* -name '*.fna' |";
while (<IN>){
	chomp;
	my $fna_file = $_;
	my ($out_name) = $fna_file =~ /^.+\/(.+?)\.fna/;
	$out_name = "TYMEFLIES_MAGs_".$out_name;
	print OUT "pilercr -in $fna_file -out Host_prediction/find_crispr_out/$out_name.pilercr_out.crisprs\n";
}
close IN;
close OUT;


open OUT, ">tmp.find_crispr_by_pilercr_2.sh";
open IN, "find /slowdata/databases/GEM/fna -name '*.fna'|";
while (<IN>){
	chomp;
	my $fna_file = $_;
	my ($out_name) = $fna_file =~ /^.+\/(.+?)\.fna/;
	$out_name = "GEM_MAGs_".$out_name;
	print OUT "pilercr -in $fna_file -out Host_prediction/find_crispr_out/$out_name.pilercr_out.crisprs\n";
}
close IN;
close OUT;

# Step 1.4 Run the crispr searching command by PILER-CR with 10 cpus
`cat tmp.find_crispr_by_pilercr_1.sh | parallel -j 10`;
`cat tmp.find_crispr_by_pilercr_2.sh | parallel -j 10`;

`rm tmp.find_crispr_by_pilercr_1.sh`;
`rm tmp.find_crispr_by_pilercr_2.sh`;







