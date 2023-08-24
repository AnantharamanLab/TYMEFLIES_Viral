#!/usr/bin/perl

use strict;
use warnings;

# AIM: Get the statistics for all MAGs from TYMEFLIES project

# Note: This script has been updated by using tax from GTDB-Tk v2.1.1 to replace the old tax

# Store the mbin_datafile_${img_id}.txt files
my %MAG2stat = (); # $mag => [0] GTDB tax [1] Genome completeness [2] Genome contamination [3] scaffolds (separated by ",")

`cat 33*/mbin*.txt | grep "^33" > tmp.all_mbin_datafiles.txt`;

open IN, "tmp.all_mbin_datafiles.txt";
while (<IN>){
	chomp;
	if (/^33/){
		my @tmp = split (/\t/);
		my $mag = $tmp[0];
		my $gtdb_tax = $tmp[3];
			
		# Change GTDB tax from Bacteria;Proteobacteria;Gammaproteobacteria;Betaproteobacteriales;Burkholderiaceae;Polynucleobacter;None
		# to d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Betaproteobacteriales;f__Burkholderiaceae;g__Polynucleobacter;s__
		my @GTDB_tax = split (/\;/,$gtdb_tax);
		my $domain = $GTDB_tax[0]; 
		if ($domain ne "None"){
			$domain = "d\_\_".$domain;
		}else{
			$domain = "d\_\_";
		}
		my $phylum = $GTDB_tax[1]; 
		if ($phylum ne "None"){
			$phylum = "p\_\_".$phylum;
		}else{
			$phylum = "p\_\_";
		}			
		my $class = $GTDB_tax[2]; 
		if ($class ne "None"){
			$class = "c\_\_".$class;
		}else{
			$class = "c\_\_";
		}				
		my $order = $GTDB_tax[3]; 
		if ($order ne "None"){
			$order = "o\_\_".$order;
		}else{
			$order = "o\_\_";
		}				
		my $family = $GTDB_tax[4]; 
		if ($family ne "None"){
			$family = "f\_\_".$family;
		}else{
			$family = "f\_\_";
		}				
		my $genus = $GTDB_tax[5]; 
		if ($genus ne "None"){
			$genus = "g\_\_".$genus;
		}else{
			$genus = "g\_\_";
		}				
		my $species = $GTDB_tax[6]; 
		if ($species ne "None"){
			$species = "s\_\_".$species;
		}else{
			$species = "s\_\_";
		}
		$gtdb_tax = "$domain;$phylum;$class;$order;$family;$genus;$species";
			
		my $completeness = $tmp[4];
		my $contamination = $tmp[5];
		my $scaffolds = $tmp[-1];
			
		$MAG2stat{$mag}[0] = $gtdb_tax;
		$MAG2stat{$mag}[1] = $completeness;
		$MAG2stat{$mag}[2] = $contamination;
		$MAG2stat{$mag}[3] = $scaffolds;
	}
}
close IN;

`rm tmp.all_mbin_datafiles.txt`;

# Store new GTDB tax
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

# Write down all stat
open OUT, ">TYMEFLIES_all_MAGs_stat.txt";
print OUT "IMG Bin ID\tGTDB-TK lineage\tBin Completeness\tBin Contamintation\tScaffolds\n";
foreach my $key (sort keys %MAG2stat){
	# Note that here we used the new GTDB tax
	print OUT "$key\t$GTDB_result_new{$key}\t$MAG2stat{$key}[1]\t$MAG2stat{$key}[2]\t$MAG2stat{$key}[3]\n";
}
close OUT;



