#!/usr/bin/perl

use strict;
use warnings;

# AIM: Get all bin member number distribution in percentage
my %Bin_member_num2freq = (); # $bin_member_num => $freq (in percentage)
my %Bin_member_num2num= (); # $bin_member_num => $num 
my $bin_num = 0; # Total bin number

open IN, "cat /storage1/data11/TYMEFLIES_phage/33*/vRhyme_best_bins_fasta_parsed/Each_bin_info.txt | sort -k 2 -n | grep -v 'un' | grep -v 'bin' | ";
while (<IN>){
	chomp;
	my $line = $_;
	my @tmp = split (/\t/,$line);
	my $bin_member_num = $tmp[1];
	$Bin_member_num2num{$bin_member_num}++;
	$bin_num++;
}
close IN;

foreach my $bin_member_num (sort keys %Bin_member_num2num){
	my $num = $Bin_member_num2num{$bin_member_num}; 
	my $freq = $num / $bin_num * 100; # Here the $freq is the percentage (%)
	$Bin_member_num2freq{$bin_member_num} = $freq;
}

# Write down result
open OUT, ">Bin_member_number_distribution.txt";
foreach my $bin_member_num (sort keys %Bin_member_num2freq){
	print OUT "Bin${bin_member_num}\t$Bin_member_num2freq{$bin_member_num}\n";
}
close OUT;

open OUT, ">Bin_member_number_distribution.num.txt";
foreach my $bin_member_num (sort keys %Bin_member_num2num){
	print OUT "Bin${bin_member_num}\t$Bin_member_num2num{$bin_member_num}\n";
}
close OUT;