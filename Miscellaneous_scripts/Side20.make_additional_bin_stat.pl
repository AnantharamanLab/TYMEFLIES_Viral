#!/usr/bin/perl

use strict;
use warnings;

# Aim: Make stat file for additional MAGs from 49 metagenomes

# Note: This script has been updated by using tax from GTDB-Tk v2.1.1 to replace the old tax

# Step 1 Find the $img_id of 49 left metagenomes
my %IMGID = (); # $img_id => 1
open IN, "ls -l | sed -r s'/ +/\t/g' | grep ^d | cut -f 9 | grep '^33' |";
while (<IN>){
	chomp;
	my $img_id = $_;
	$IMGID{$img_id} = 1;
}
close IN;

my %IMGID_left = %IMGID; # $img_id => 1
open IN, "/storage1/data11/TYMEFLIES_phage/Binning_Data/TYMEFLIES_all_MAGs_stat.txt";
while (<IN>){
	chomp;
	if (!/IMG/){
		my @tmp = split (/\t/);
		my ($img_id) = $tmp[0] =~ /^(.+?)\_/;
		delete $IMGID_left{$img_id};
	}
}
close IN;

# Step 2 Store CheckM result
my %MAG2checkM_result = (); # $mag => $completeness and $contamination, separated by "\t"
foreach my $img_id (sort keys %IMGID_left){
	open IN, "/storage1/data11/TYMEFLIES_phage/Binning_Data/$img_id/checkm_result.txt";
	while (<IN>){
		chomp;
		my $line = $_;
		if ($line =~ /root/ or $line =~ /__/){
			$line =~ s/(^\s+|\s+$)//g;
			$line =~ s/ +/ /g;
			my @tmp = split (/\s/,$line);
			my $mag = $tmp[0];
			my $completeness = $tmp[12];
			my $contamination = $tmp[13];
			$MAG2checkM_result{$mag} = "$completeness\t$contamination";
		}
	}
	close IN;
}

# Step 3 Store new GTDB tax
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

# Step 4 Make additional MAG stat file
my %MAG2stat = (); # $mag => [0] $tax [1] $completeness [2] $contamination [3] $scfs
my %MAGs_to_delete = (); # $mag => 1
foreach my $mag (sort keys %MAG2checkM_result){
	my ($completeness, $contamination) = $MAG2checkM_result{$mag} =~ /^(.+?)\t(.+?)$/;
	if ($completeness >= 50 and $contamination <= 10){
		my $tax = $GTDB_result_new{$mag};
		my @Scfs = (); # Store all the scfs within this MAG
		my ($img_id) = $mag =~ /^(.+?)\_/;
		open IN, "/storage1/data11/TYMEFLIES_phage/Binning_Data/$img_id/$mag.fna";
		while (<IN>){
			chomp;
			if (/^>/){
				my ($scf) = $_ =~ /^>(.+?)$/;
				push @Scfs, $scf;
			}
		}
		close IN;
		my $scfs = join("\,",@Scfs);
		$MAG2stat{$mag}[0] = $tax;
		$MAG2stat{$mag}[1] = $completeness;
		$MAG2stat{$mag}[2] = $contamination;
		$MAG2stat{$mag}[3] = $scfs;
	}else{
		$MAGs_to_delete{$mag} = 1;
	}
}

# Step 5 Write down additional MAG stat file
open OUT, ">/storage1/data11/TYMEFLIES_phage/Binning_Data/TYMEFLIES_all_MAGs_stat.additional_49_metagenomes.txt";
print OUT "IMG Bin ID\tGTDB-TK lineage\tBin Completeness\tBin Contamintation\tScaffolds\n";
foreach my $mag (sort keys %MAG2stat){
	print OUT "$mag\t$MAG2stat{$mag}[0]\t$MAG2stat{$mag}[1]\t$MAG2stat{$mag}[2]\t$MAG2stat{$mag}[3]\n";
}
close OUT;
=pod
# Step 6 Write down which MAGs to delete
open OUT, ">MAGs_to_delete.txt";
foreach my $mag (sort keys %MAGs_to_delete){
	my ($img_id) = $mag =~ /^(.+?)\_/;
	print OUT "rm /storage1/data11/TYMEFLIES_phage/Binning_Data/$img_id/$mag.fna\n";
}
close OUT;

