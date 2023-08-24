#!/usr/bin/perl
  
use strict;
use warnings;

# Aim: Make the MAG family to full family tax map

# Step 1 Store the new GTDB-Tk result
my %GTDB_result_new = (); # $mag => $tax
my %Tax_full = (); # $tax => 1; Store the dereplicated full tax information
open IN, "/storage1/data11/TYMEFLIES_phage/Binning_Data/All_MAGs_gtdbtk_result/gtdbtk.bac120.summary.tsv";
while (<IN>){
	chomp;
	if (!/^user/){
		my @tmp = split (/\t/);
		my $mag = $tmp[0];
		my $tax = $tmp[1];
		$GTDB_result_new{$mag} = $tax;
		$Tax_full{$tax} = 1;
	}
}
close IN;

# Step 2 Make the MAG family to full family tax map hash
my %MAG_family2full_family = (); # $family => $full_family
foreach my $tax_full (sort keys %Tax_full){
	my @Tax_full = split (/\;/, $tax_full);
	my $family = $Tax_full[3]."\;".$Tax_full[4]; # for instance, "o__Burkholderiales;f__Burkholderiaceae"
	my $full_family = $Tax_full[0]."\;".$Tax_full[1]."\;".$Tax_full[2]."\;".$Tax_full[3]."\;".$Tax_full[4];
	if ($family ne "o__;f__"){
		$MAG_family2full_family{$family} = $full_family;
	}
}

# Step 3 Write down MAG_family2full_family_tax_map.txt
open OUT, ">Binning_Data/MAG_family2full_family_tax_map.txt";
foreach my $family (sort keys %MAG_family2full_family){
	print OUT "$family\t$MAG_family2full_family{$family}\n";
}
close OUT;
