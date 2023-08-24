#!/usr/bin/perl

use strict;
use warnings;

# Aim: Find non-cyanobacteria MAGs that connect with cyanophage

# Step 1 Store the list of MAGs from "by_crispr_match.txt"
my %MAGs = (); # $mag => 1
open IN, "Check_non-cyanobacteria/Viral_gn_with_psbAD2host_mag.by_crispr_match.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my @tmp2 = split (/\|/, $tmp[1]);
	my @tmp3 = split (/\,/, $tmp2[0]);
	
	foreach my $mag (@tmp3){
		$MAGs{$mag} = 1;
	}
}
close IN;

# Step 2 Write down the list of MAGs
open OUT, ">Check_non-cyanobacteria/Non-cyanobacteria_MAGs_that_connect_with_cyanophage.by_crispr_match.txt";
foreach my $mag (sort keys %MAGs){
	print OUT "$mag\n";
}
close OUT;

# Step 3 Store the list of MAGs from "by_sequence_identity"
my %MAGs2 = (); # $mag => 1
open IN, "Check_non-cyanobacteria/Viral_gn_with_psbAD2host_mag.by_sequence_identity.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my @tmp2 = split (/\|/, $tmp[1]);
	my @tmp3 = split (/\,/, $tmp2[0]);
	
	foreach my $mag (@tmp3){
		$MAGs2{$mag} = 1;
	}
}
close IN;

# Step 4 Write down the list of MAGs
open OUT, ">Check_non-cyanobacteria/Non-cyanobacteria_MAGs_that_connect_with_cyanophage.by_sequence_identity.txt";
foreach my $mag (sort keys %MAGs2){
	print OUT "$mag\n";
}
close OUT;


