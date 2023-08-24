#!/usr/bin/perl

use strict;
use warnings;

my %Depth = (); # $scaf =>  [0] bowtie2_depth [1] img_depth
my $img_id = "3300042358";
open IN, "$img_id/$img_id\.covstat";
while (<IN>){
	chomp;
	if (!/^#/){
		my @tmp = split (/\t/);
		my $scaf = $tmp[0];
		my $bowtie2_depth = $tmp[1];
		$Depth{$scaf}[0] = $bowtie2_depth;
	}
}
close IN;

open IN, "$img_id/$img_id\.a\.depth\.txt";
while (<IN>){
	chomp;
	if (!/^ID/){
		my @tmp = split (/\t/);
		my $scaf = $tmp[0];
		my $img_depth = $tmp[1];
		$Depth{$scaf}[1] = $img_depth;
	}
}
close IN;

open OUT, ">$img_id\.depth\_compare\.txt";
foreach my $scaf (sort keys %Depth){
	print OUT "$scaf\t$Depth{$scaf}[0]\t$Depth{$scaf}[1]\n";
}
close OUT;

