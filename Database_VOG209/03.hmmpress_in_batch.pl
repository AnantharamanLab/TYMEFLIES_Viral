#!/usr/bin/perl

use strict;
use warnings;

open OUT, ">tmp.run_hmmpress_in_batch.sh";
open IN, "ls *.hmm |";
while (<IN>){
	chomp;
	my $hmm = $_;
	#my ($hmm_name) = $hmm =~ /^(.+?)\.hmm/;
	print OUT "hmmpress $hmm\n";
}
close IN;
close OUT;

`cat tmp.run_hmmpress_in_batch.sh | parallel -j 20`;

`rm tmp.run_hmmpress_in_batch.sh`;