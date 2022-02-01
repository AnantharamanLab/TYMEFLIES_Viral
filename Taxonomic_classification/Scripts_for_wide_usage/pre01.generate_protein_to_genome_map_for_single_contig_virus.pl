#!/usr/bin/perl

use strict;
use warnings;

# Aim: Generate protein to genome map for single-contig virus

# Usage: The input should be the all virus protein file (ended with "faa")
# perl pre01.generate_protein_to_genome_map_for_single_contig_virus.pl [virus_protein_sequence.faa] [map_file.txt] 
# The first item is the input, the second item is the output (the map file)

# Note: The input faa file should be a prodigal-formate protein sequence file
# Note: Note that the virus genome is made up of just one contig/scaffold

my $input = $ARGV[0];
my $output = $ARGV[1];

# Step 1 Store all protein to scaffold map
my %Seq= _store_seq("$input");

my %Pro2scf = (); # $pro => $scf
foreach my $key (sort keys %Seq){
	my ($key_clean) = $key =~ /^>(.+?)$/;
	my ($scf) = $key_clean =~ /^(.+)\_/;
	$Pro2scf{$key_clean} = $scf;
}

# Step 2 Store the map file
open OUT, ">$output";
foreach my $pro (sort keys %Pro2scf){
	print OUT "$pro\t$Pro2scf{$pro}\n";
}
close OUT;



## Subroutines
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