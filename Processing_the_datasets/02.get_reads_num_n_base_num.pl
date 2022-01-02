#!/usr/bin/perl

use strict;
use warnings;

# need to activate conda env "BBTools_v37.62" first (will use reformat.sh)

`mkdir tmp.run_reformat.intermediate_folder`;

open OUT, ">tmp.run_reformat.sh";
open IN, "fastq_list.2.txt";
while (<IN>){
	chomp;
	my ($fq_name) = $_ =~ /^(.+?)\.fq\.gz/;
	print OUT "reformat.sh in=$fq_name.fq.gz 2> tmp.run_reformat.intermediate_folder/$fq_name.reformat.out\n";
}
close IN;
close OUT;

`cat tmp.run_reformat.sh | parallel -j 40`;

`rm tmp.run_reformat.sh`;

my %Result = (); # $fq_name => [0] read num [1] read base num
open IN, "ls tmp.run_reformat.intermediate_folder/*.out |";
while (<IN>){
	chomp;
	my $file = $_;
	my ($fq_name) = $file =~ /\/(.+?)\.reformat\.out/;
	open INN, $file;
	while (<INN>){
		chomp;
		if (/^Input/){
			my $line = $_;
			my ($read_num,$read_base_num) = $line =~ /\s(\d+?)\sreads\s+?(\d+?)\sbases/;
			$Result{$fq_name}[0] = $read_num;
			$Result{$fq_name}[1] = $read_base_num;
		}
	}
	close INN;
}
close IN;

# print output table
open OUT, ">Read_num_n_base_num.2.txt";
foreach my $fq_name (sort keys %Result){
	print OUT "$fq_name\t$Result{$fq_name}[0]\t$Result{$fq_name}[1]\n";
}
close OUT;

`rm -r tmp.run_reformat.intermediate_folder`;

