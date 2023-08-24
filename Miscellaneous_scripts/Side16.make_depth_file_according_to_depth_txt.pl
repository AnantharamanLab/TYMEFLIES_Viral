#!/usr/bin/perl

use strict;
use warnings;

# Aim: Make depth file according to "a.depth.txt" file in 433 metagenome folders

open IN, "find /storage1/data11/TYMEFLIES_phage/ -maxdepth 2 -type f -name '*.a.depth.txt' | ";
while (<IN>){
	chomp;
	my $file = $_;
	my ($img_id) = $file =~ /^.+\/(.+?)\.a\.depth\.txt/;
	
	# Step 1 Store the covstat file info
	my %Scf2depth = (); # $scf => $depth
	open INN, "$file";
	while (<INN>){
		chomp;
		if (!/^ID/){
			my @tmp = split (/\t/);
			$Scf2depth{$tmp[0]} = $tmp[1];
		}
	}
	close INN;	
	
	# Step 2 Store the scf length info
	my %Scf2length = (); # $scf => $length
	open INN, "/storage1/data11/TYMEFLIES_phage/$img_id/$img_id.a.fna.seq_length.txt";
	while (<INN>){
		chomp;
		if (!/^sequence/){
			my @tmp = split (/\t/);
			$Scf2length{$tmp[0]} = $tmp[1];
		}
	}
	close INN;	
	
	# Step 3 Make depth file 
	open OUT, ">/storage1/data11/TYMEFLIES_phage/$img_id/$img_id.for_metabat.depth";
	print OUT "contigName\tcontigLen\ttotalAvgDepth\n";
	foreach my $scf (sort keys %Scf2depth){
		print OUT "$scf\t$Scf2length{$scf}\t$Scf2depth{$scf}\n";
	}
	close OUT;
}
close IN;



