#!/usr/bin/perl

use strict;
use warnings;

my %Fastq = (); # full path to each fastq gz file
# $fastq => $folder_name
open IN, "fastq_file_list.txt"; # contains the first 180 fastq gz files
while (<IN>){
	chomp;
	my @tmp = split (/ /);
	my $fastq_tmp = $tmp[-1];
	my $fastq = "/mnt/researchdrive/TYMEFLIES_fastq/".$fastq_tmp;
	my ($folder_name) = $fastq =~ /\/mnt\/researchdrive\/TYMEFLIES\_fastq\/(.+?)\/Sequence/;
	$Fastq{$fastq} = $folder_name;
}
close IN;

# store the map
my %Reads2IMG_ID_map = (); # $folder_name => $img_id
open IN, "Reads2IMG_ID_map.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $img_id = $tmp[0];
	my $folder_name = $tmp[2];
	$Reads2IMG_ID_map{$folder_name} = $img_id;
}
close IN;

# print batch copy and rename script
open OUT, ">tmp.batch_copy_and_rename.sh";
foreach my $fastq (sort keys %Fastq){
	my $folder_name = $Fastq{$fastq};
	my $img_id = $Reads2IMG_ID_map{$folder_name};
	print OUT "cp $fastq $img_id.filtered.fq.gz\n";
}
close OUT;

`cat tmp.batch_copy_and_rename.sh | parallel -j 20`;

`rm tmp.batch_copy_and_rename.sh`;
