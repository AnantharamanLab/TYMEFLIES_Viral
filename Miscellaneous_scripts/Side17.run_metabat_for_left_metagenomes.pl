#!/usr/bin/perl

use strict;
use warnings;

# Aim: Run metabat for scaffolds in 49 left metagenome folders

# Step 1 Find the $img_id of 49 left metagenomes
my %IMGID = (); # $img_id => 1
open IN, "ls -l | sed -r s'/ +/\t/g' | grep ^d | cut -f 9 | grep '^33' |";
while (<IN>){
	chomp;
	my $img_id = $_;
	$IMGID{$img_id} = 1;
}
close IN;

my %IMGID_left = %IMGID; # $img_id => 1
open IN, "/storage1/data11/TYMEFLIES_phage/Binning_Data/TYMEFLIES_all_MAGs_stat.txt";
while (<IN>){
	chomp;
	if (!/IMG/){
		my @tmp = split (/\t/);
		my ($img_id) = $tmp[0] =~ /^(.+?)\_/;
		delete $IMGID_left{$img_id};
	}
}
close IN;

# Step 2 Run metabat for 49 metagenomes with depth file as "$img_id.for_metabat.depth"
open OUT, ">tmp.run_metabat_in_parallel.sh";
foreach my $img_id (sort keys %IMGID_left){
	`mkdir /storage1/data11/TYMEFLIES_phage/Binning_Data/$img_id`;
	print OUT "/storage1/data11/metabat/metabat -i $img_id/$img_id.a.fna -a $img_id/$img_id.for_metabat.depth -o /storage1/data11/TYMEFLIES_phage/Binning_Data/$img_id/$img_id --superspecific -v -m 3000 -t 1\n";	
}
close OUT;

`cat tmp.run_metabat_in_parallel.sh | parallel -j 10`;

`rm tmp.run_metabat_in_parallel.sh`;









