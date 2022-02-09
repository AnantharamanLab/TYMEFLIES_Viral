#!/usr/bin/perl

use strict;
use warnings;

# Step 1 Get list of 587 VOG markers
my %VOG_marker = (); # $vog => 1;
open IN, "VOG_marker_table.txt";
while (<IN>){
	chomp;
	if (!/VOG\tFunction/){
		my @tmp = split (/\t/);
		$VOG_marker{$tmp[0]} = 1;
	}
}
close IN;

open OUT, ">VOG_marker_list.txt";
foreach my $vog (sort keys %VOG_marker){
	print OUT "$vog\n";
}
close OUT;

# Step 1 Get a subset of HMM of 587 VOG markers
`hmmfetch -o VOGDB97.587_marker_mdfed.HMM -f VOGDB97.HMM VOG_marker_list.txt`;

`hmmpress VOGDB97.587_marker_mdfed.HMM`;

`rm VOG_marker_list.txt`;
