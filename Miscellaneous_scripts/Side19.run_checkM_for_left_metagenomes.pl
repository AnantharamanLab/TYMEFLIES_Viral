#!/usr/bin/perl

use strict;
use warnings;

# Aim: Run checkM for genomes in 49 left metagenome folders
# Note: This should be run under conda env by "conda activate CheckM_v1.1.3"

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

# Step 2 Run checkM for 49 metagenomes 
open OUT, ">tmp.run_checkM_in_parallel.sh";
foreach my $img_id (sort keys %IMGID_left){
	print OUT "mkdir /storage1/data11/TYMEFLIES_phage/Binning_Data/$img_id/tmpdir_checkm;";
	print OUT "checkm lineage_wf -x .fna /storage1/data11/TYMEFLIES_phage/Binning_Data/$img_id  /storage1/data11/TYMEFLIES_phage/Binning_Data/$img_id/checkm_result -f /storage1/data11/TYMEFLIES_phage/Binning_Data/$img_id/checkm_result.txt -t 1 --tmpdir /storage1/data11/TYMEFLIES_phage/Binning_Data/$img_id/tmpdir_checkm;";
	print OUT "rm -r /storage1/data11/TYMEFLIES_phage/Binning_Data/$img_id/tmpdir_checkm\n";
	
	
}
close OUT;

`cat tmp.run_checkM_in_parallel.sh | parallel -j 10`;

`rm tmp.run_checkM_in_parallel.sh`;







