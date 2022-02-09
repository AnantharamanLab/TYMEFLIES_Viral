#!/usr/bin/perl

use strict;
use warnings;

# AIM: Check if all the 587 VOG marker HMMs exisit and make a subset of HMM

# Step 1. To store all the VOG marker HMMs
my %VOG_marker = (); # $vog => 1;
open IN, "VOG_marker_table.txt";
while (<IN>){
	chomp;
	if (!/VOG\tFunction/){
		my @tmp = split (/\t/);
		$VOG_marker{$tmp[0]} = 1;
	}
}
close IN;

# Step 2. To check whether all the VOG markers are present
my %VOGDB_HMM = (); # $vog => 1
open IN, "VOGDB97.HMM";
while (<IN>){
	chomp;
	if (/^NAME/){
		my ($vog) = $_ =~ /(VOG\d+?)$/;
		$VOGDB_HMM{$vog} = 1;
	}
}
close IN;

foreach my $vog (sort keys %VOG_marker){
	if (! exists $VOGDB_HMM{$vog}){
		print "$vog\.hmm is not present in VOGDB97.HMM\n";
	}
}

