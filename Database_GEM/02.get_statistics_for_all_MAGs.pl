#!/usr/bin/perl

use strict;
use warnings;

# AIM: Get the statistics for all MAGs from GEM project

my %MAG2stat = (); # $mag => [0] GTDB tax [1] Genome completeness [2] Genome contamination [3] scaffolds (separated by ",")
open IN, "genome_metadata.tsv";
while (<IN>){
	chomp;
	if (/^33/){
		my @tmp = split(/\t/);
		my $mag = $tmp[0];
		my $gtdb_tax = $tmp[14];	
		my $completeness = $tmp[9];
		my $contamination = $tmp[10];
		my $scaffolds = ""; 
		
		# Store scaffolds from each MAG
		my %Seq = _store_seq("/slowdata/databases/GEM/fna/$mag.fna");
		my @Scaffolds = (); 
		foreach my $key (sort keys %Seq){
			my ($key_clean) = $key =~ /^>(.+?)$/;
			push @Scaffolds, $key_clean;
		}
		$scaffolds = join(',',@Scaffolds);
		
		$MAG2stat{$mag}[0] = $gtdb_tax;
		$MAG2stat{$mag}[1] = $completeness;
		$MAG2stat{$mag}[2] = $contamination;
		$MAG2stat{$mag}[3] = $scaffolds;		
	}
}
close IN;

# Write down all stat
open OUT, ">GEM_all_MAGs_stat.txt";
print OUT "IMG Bin ID\tGTDB-TK lineage\tBin Completeness\tBin Contamintation\tScaffolds\n";
foreach my $key (sort keys %MAG2stat){
	print OUT "$key\t$MAG2stat{$key}[0]\t$MAG2stat{$key}[1]\t$MAG2stat{$key}[2]\t$MAG2stat{$key}[3]\n";
}
close OUT;



# Subroutine

sub _store_seq{
	my $file = $_[0];
	my %Seq = (); my $head = "";
	open _IN, "$file";
	while (<_IN>){
		chomp;
		if (/>/){
			if (/\s/){
				($head) = $_ =~ /^(>.+?)\s/;
				$Seq{$head} = "";
			}else{
				($head) = $_ =~ /^(>.+?)$/;
				$Seq{$head} = "";
			}
		}else{
			$Seq{$head} .= $_;
		}
	}
	close _IN;
	return %Seq;
}