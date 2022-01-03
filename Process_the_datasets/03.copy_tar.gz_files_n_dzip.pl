#!/usr/bin/perl

use strict;
use warnings;

my %File = (); # $file => $img_id; the absolute path to each tar.gz file
open IN, "/mnt/researchdrive/TYMEFLIES_faa_fna_mags_etc/tar.gz.file_list.txt";
while (<IN>){
	chomp;
	my $file = $_;
	my ($img_id) = $file =~ /IMG\_Data\/(.+?)\.tar\.gz/;
	$File{$file} = $img_id;
}
close IN;

# batch copy and dzip
open OUT, ">tmp.batch_copy_n_dzip.sh";
foreach my $file (sort keys %File){
	my $img_id = $File{$file};
	print OUT "cp $file ./; tar xzf $img_id.tar.gz\n";
}
close OUT;

`cat tmp.batch_copy_n_dzip.sh | parallel -j 30`;

`rm tmp.batch_copy_n_dzip.sh`;


