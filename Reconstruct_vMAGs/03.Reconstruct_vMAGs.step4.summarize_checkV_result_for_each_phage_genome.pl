#!/usr/bin/perl

use strict;
use warnings;

# NOTE: This script should be run after 11.run_checkV_for_each_phage_genome.pl run has been completed

# Store all metagenomes
my %IMGID = (); # $img_id => 1
open IN, "ls -l | sed -r s'/ +/\t/g' | grep ^d | cut -f 9 | grep '^33' |";
while (<IN>){
	chomp;
	my $img_id = $_;
	$IMGID{$img_id} = 1;
}
close IN;

# TEST script
#%IMGID = ();
#$IMGID{"3300020486"} = 1;

# Store and rewrite CheckV result
my %Contamination = (); # $bin => $line; $bin contains $img_id in the front
my $contamination_table_header = ""; # The header line of contamination.tsv

my %Completeness = (); # $bin => $line; $bin contains $img_id in the front
my $completeness_table_header = ""; # The header line of completeness.tsv

my %Complete_genomes = (); # $bin => $line; $bin contains $img_id in the front
my $complete_genomes_table_header = ""; # The header line of complete_genomes.tsv
	
my %Quality_summary = (); # $bin => $line; $bin contains $img_id in the front
my $quality_summary_table_header = ""; # The header line of quality_summary.tsv

foreach my $img_id (sort keys %IMGID){	
	# Delete tmp
	`rm -r $img_id/CheckV_phage_bin/tmp`;
	
	# Store "contamination" info
	open INN, "$img_id/CheckV_phage_bin/contamination.tsv";
	while (<INN>){
		chomp;
		if (/^contig_id/){
			$contamination_table_header = $_;
		}else{
			my $line = $_;
			my @tmp = split (/\t/, $line);
			my $bin = $tmp[0];
			$Contamination{$bin} = $line;
		}
	}
	close INN;
	
	# Store "completeness" info
	open INN, "$img_id/CheckV_phage_bin/completeness.tsv";
	while (<INN>){
		chomp;
		if (/^contig_id/){
			$completeness_table_header = $_;
		}else{
			my $line = $_;
			my @tmp = split (/\t/, $line);
			my $bin = $tmp[0];
			$Completeness{$bin} = $line;
		}
	}
	close INN;

	# Store "complete_genomes" info
	open INN, "$img_id/CheckV_phage_bin/complete_genomes.tsv";
	while (<INN>){
		chomp;
		if (/^contig_id/){
			$complete_genomes_table_header = $_;
		}else{
			my $line = $_;
			my @tmp = split (/\t/, $line);
			my $bin = $tmp[0];
			$Complete_genomes{$bin} = $line;
		}
	}
	close INN;
	



	# Store "quality_summary" info
	open INN, "$img_id/CheckV_phage_bin/quality_summary.tsv";
	while (<INN>){
		chomp;
		if (/^contig_id/){
			$quality_summary_table_header = $_;
		}else{
			my $line = $_;
			my @tmp = split (/\t/, $line);
			my $bin = $tmp[0];
			$Quality_summary{$bin} = $line;
		}
	}
	close INN;
}

`mkdir CheckV_phage_bin_all`;
# Print "contamination" info
open OUT, ">CheckV_phage_bin_all/contamination.tsv";
print OUT $contamination_table_header."\n";
foreach my $key (sort keys %Contamination){
	print OUT "$Contamination{$key}\n";
}
close OUT;
	
# Print "completeness" info
open OUT, ">CheckV_phage_bin_all/completeness.tsv";
print OUT $completeness_table_header."\n";
foreach my $key (sort keys %Completeness){
	print OUT "$Completeness{$key}\n";
}
close OUT;	
	
# Print "complete_genomes" info
open OUT, ">CheckV_phage_bin_all/complete_genomes.tsv";
print OUT $complete_genomes_table_header."\n";
foreach my $key (sort keys %Complete_genomes){	
	print OUT "$Complete_genomes{$key}\n";
}
close OUT;
	
# Print "quality_summary" info
open OUT, ">CheckV_phage_bin_all/quality_summary.tsv";
print OUT $quality_summary_table_header."\n";
foreach my $key (sort keys %Quality_summary){
	print OUT "$Quality_summary{$key}\n";
}
close OUT;	