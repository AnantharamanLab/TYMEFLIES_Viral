#!/usr/bin/perl

use strict;
use warnings;

my %File = (); # $file => $file_name; the absolute path to each tar.gz file
open IN, "/mnt/researchdrive/TYMEFLIES_faa_fna_mags_etc/Binning_Data.tar.gz.file_list.txt";
while (<IN>){
	chomp;
	my $file = $_;
	my ($file_name) = $file =~ /Binning\_Data\/(.+?)\.tar\.gz/;
	$File{$file} = $file_name;
}
close IN;

# batch copy and dzip
open OUT, ">tmp.batch_copy_n_dzip.sh";
foreach my $file (sort keys %File){
	my $file_name = $File{$file};
	print OUT "cp $file ./; tar xzf $file_name.tar.gz\n";
}
close OUT;

`cat tmp.batch_copy_n_dzip.sh | parallel -j 10`;

`rm tmp.batch_copy_n_dzip.sh`;


