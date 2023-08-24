#!/usr/bin/perl
  
use strict;
use warnings;

# Aim: To see how many (or percentage) more MAGs get species level assignment
# after using the lastest GTDB-Tk (v2.1.1)

# Step 1 Store the old GTDB-Tk result
my %GTDB_result_old = (); # $mag => [0] $tax [1] $completeness [2] $contamination [3] $scfs 
open IN, "/storage1/data11/TYMEFLIES_phage/Binning_Data/TYMEFLIES_all_MAGs_stat.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		chomp;
		my @tmp = split (/\t/);
		my $mag = $tmp[0];
		my $tax = $tmp[1];
		my $completeness = $tmp[2];
		my $contamination = $tmp[3];
		my $scfs = $tmp[4];
		$GTDB_result_old{$mag}[0] = $tax;
		$GTDB_result_old{$mag}[1] = $completeness;
		$GTDB_result_old{$mag}[2] = $contamination;
		$GTDB_result_old{$mag}[3] = $scfs;
	}
}
close IN;

open IN, "/storage1/data11/TYMEFLIES_phage/Binning_Data/TYMEFLIES_all_MAGs_stat.additional_49_metagenomes.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		chomp;
		my @tmp = split (/\t/);
		my $mag = $tmp[0];
		my $tax = $tmp[1];
		my $completeness = $tmp[2];
		my $contamination = $tmp[3];
		my $scfs = $tmp[4];
		$GTDB_result_old{$mag}[0] = $tax;
		$GTDB_result_old{$mag}[1] = $completeness;
		$GTDB_result_old{$mag}[2] = $contamination;
		$GTDB_result_old{$mag}[3] = $scfs;
	}
}
close IN;

# Step 2 Store the new GTDB-Tk result
my %GTDB_result_new = (); # $mag => $tax
open IN, "/storage1/data11/TYMEFLIES_phage/Binning_Data/All_MAGs_gtdbtk_result/gtdbtk.bac120.summary.tsv";
while (<IN>){
	chomp;
	if (!/^user/){
		my @tmp = split (/\t/);
		my $mag = $tmp[0];
		my $tax = $tmp[1];
		$GTDB_result_new{$mag} = $tax;
	}
}
close IN;

# Step 3 Compare the two results and write down the comparison result
my $total_mags_num = 0;
my $num_old_gtdbtk_tax_down_to_species_level = 0; # The number of MAGs that have species level GTDB-Tk tax (for old GTDB-Tk assignment)
my $num_new_gtdbtk_tax_down_to_species_level = 0; # The number of MAGs that have species level GTDB-Tk tax (for new GTDB-Tk assignment)

foreach my $mag (sort keys %GTDB_result_old){
	$total_mags_num++;
	# See the old GTDB tax species 
	my $tax_old = $GTDB_result_old{$mag}[0];
	my @Tax_old = split (/\;/, $tax_old);
	my $tax_old_species = $Tax_old[-1];
	if ($tax_old_species ne "s\_\_"){
		$num_old_gtdbtk_tax_down_to_species_level++;
	}
	# See the new GTDB tax species 
	my $tax_new = $GTDB_result_new{$mag};
	my @Tax_new = split (/\;/, $tax_new);
	my $tax_new_species = $Tax_new[-1];
	if ($tax_new_species ne "s\_\_"){
		$num_new_gtdbtk_tax_down_to_species_level++;
	}
}

open OUT, ">/storage1/data11/TYMEFLIES_phage/Binning_Data/New_old_GTDB_tax_assignment_comparison.txt";
print OUT "Total MAG\tNumber of MAGs down to species level for old\tNumber of MAGs down to species level for new\n";
print OUT "$total_mags_num\t$num_old_gtdbtk_tax_down_to_species_level\t$num_new_gtdbtk_tax_down_to_species_level\n";
close OUT;

