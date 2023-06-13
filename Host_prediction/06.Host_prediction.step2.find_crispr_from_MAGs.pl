#!/usr/bin/perl

use strict;
use warnings;

# AIM: Find crispr spacers from MAGs (both TYMEFLIES and GEM MAGs)
# NOTE: To run MinCED, it should conda activate it first before run this script: conda activate MinCED
# NOTE: To run PILER-CR, it should conda activate it first before run this script: conda activate PILER-CR_v1.06

# Step 1 Write down and run crispr searching command for TYMEFLIES MAGs
`mkdir Host_prediction/find_crispr_out`;

## Step 1.1 Write down crispr searching command by MinCED for TYMEFLIES MAGs and GEM MAGs
open OUT, ">tmp.find_crispr_by_minced_1.sh";
open IN, "find /storage1/data11/TYMEFLIES_phage/Robin_MAGs/All_passed_MAGs -name '*.fasta' |";
while (<IN>){
	chomp;
	my $fasta_file = $_;
	my ($out_name) = $fasta_file =~ /^.+\/(.+?)\.fasta/;
	$out_name = "TYMEFLIES_MAGs_".$out_name;
	print OUT "minced $fasta_file Host_prediction/find_crispr_out/$out_name.minced_out.crisprs\n";
}
close IN;
close OUT;

# Step 1.2 Run the crispr searching command by MinCED with 20 cpus
`cat tmp.find_crispr_by_minced_1.sh | parallel -j 20`;
`rm tmp.find_crispr_by_minced_1.sh`;



## Step 1.3 Write down crispr searching command by PILER-CR for TYMEFLIES MAGs
open OUT, ">tmp.find_crispr_by_pilercr_1.sh";
open IN, "find /storage1/data11/TYMEFLIES_phage/Robin_MAGs/All_passed_MAGs -name '*.fasta' |";
while (<IN>){
	chomp;
	my $fasta_file = $_;
	my ($out_name) = $fasta_file =~ /^.+\/(.+?)\.fasta/;
	$out_name = "TYMEFLIES_MAGs_".$out_name;
	print OUT "pilercr -in $fasta_file -out Host_prediction/find_crispr_out/$out_name.pilercr_out.crisprs\n";
}
close IN;
close OUT;

# Step 1.4 Run the crispr searching command by PILER-CR with 20 cpus
`cat tmp.find_crispr_by_pilercr_1.sh | parallel -j 20`;
`rm tmp.find_crispr_by_pilercr_1.sh`;

