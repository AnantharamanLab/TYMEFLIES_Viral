#!/usr/bin/perl

use strict;
use warnings;

# Aim: Run gtdbtk for scaffolds in 49 left metagenome folders
# Note: This should be run under conda env by "conda activate gtdbtk_v0.2.2"

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

# Step 2 Run gtdbtk for all MAGs in 49 metagenomes
`mkdir /storage1/data11/TYMEFLIES_phage/Binning_Data/IMGID_left_bins`;
foreach my $img_id (sort keys %IMGID_left){
	`cp /storage1/data11/TYMEFLIES_phage/Binning_Data/$img_id/*.fna /storage1/data11/TYMEFLIES_phage/Binning_Data/IMGID_left_bins`;
}

open OUT, ">tmp.run_gtdbtk.sh";
foreach my $img_id (sort keys %IMGID_left){
	print OUT "gtdbtk classify_wf --genome_dir /storage1/data11/TYMEFLIES_phage/Binning_Data/IMGID_left_bins --out_dir /storage1/data11/TYMEFLIES_phage/Binning_Data/IMGID_left_bins/gtdbtk_output_dir -x .fna --cpus 20\n";
}
close OUT;

`bash tmp.run_gtdbtk.sh`;

`rm tmp.run_gtdbtk.sh`;







