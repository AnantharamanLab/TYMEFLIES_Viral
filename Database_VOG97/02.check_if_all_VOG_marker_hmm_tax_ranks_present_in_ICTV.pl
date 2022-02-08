#!/usr/bin/perl

use strict;
use warnings;

# AIM: Check if all the 587 VOG marker HMMs exisit

# Step 1. To store all the VOG marker HMMs and tax
my %VOG_marker = (); # $vog => $tax;
my %Tax = (); # $tax => 1;
my %Ranks = (); # $rank => 1; # Here only ranks over family level are included
open IN, "VOG_marker_table.txt";
while (<IN>){
	chomp;
	if (!/VOG\tFunction/){
		my @tmp = split (/\t/);
		$VOG_marker{$tmp[0]} = $tmp[2];
		$Tax{$tmp[2]} = 1;
		
		my @Ranks = split (/\;/, $tmp[2]);
		foreach my $rank (@Ranks){
			if ($rank ne "NA"){
				$Ranks{$rank} = 1;
			}
		}
	}
}
close IN;

# Step 2. To store all the ranks in ICTV
my %ICTV_ranks = (); # $rank => 1
                     # Here only ranks over family level are included
open IN, "/slowdata/databases/NCBI_RefSeq_viral/ICTV_Master_Species_List_2020.v1.txt";
while (<IN>){
	chomp;
	if (!/^Sort/){
		my @tmp = split (/\t/);
		$ICTV_ranks{$tmp[1]} = 1; # Realm
		$ICTV_ranks{$tmp[3]} = 1; # Kingdom
		$ICTV_ranks{$tmp[5]} = 1; # Phylum
		$ICTV_ranks{$tmp[7]} = 1; # Class
		$ICTV_ranks{$tmp[9]} = 1; # Order
		$ICTV_ranks{$tmp[11]} = 1; # Family
	}
}
close IN;

# Step 3. Compare two hash to find VOG marker tax ranks are not present in ICTV (MSL #36) ranks
foreach my $rank (sort keys %Ranks){
	if (! exists $ICTV_ranks{$rank}){
		print "Tax rank $rank from \"VOG_marker_table.txt\" is not present in ICTV (MSL #36)\n";
	}
}
