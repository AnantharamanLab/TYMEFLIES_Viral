#!/usr/bin/perl

use strict;
use warnings;

# NOTE: This script should be run within conda env "CheckV_v0.8.1"
# NOTE: Before using, do this command to indicate database: export CHECKVDB=/slowdata/databases/checkv-db-v0.6

# Step 1. Link bin scaffolds with Ns to make temporary bin fasta file
# Step 2. Run checkV for temporary bin fasta file

# Store all metagenomes
my %IMGID = (); # $img_id => 1
open IN, "ls -l | sed -r s'/ +/\t/g' | grep ^d | cut -f 9 | grep '33' |";
while (<IN>){
	chomp;
	my $img_id = $_;
	$IMGID{$img_id} = 1;
}
close IN;

# TEST script
#%IMGID = ();
#$IMGID{"3300020486"} = 1;

# Step 1. Link bin scaffolds with Ns to make temporary bin fasta file
foreach my $img_id (sort keys %IMGID){
	`python3 /storage1/data11/TYMEFLIES_phage/vRhyme_v1.0.0/vRhyme/aux/link_bin_sequences.py -i $img_id/vRhyme_best_bins_fasta_parsed -e fasta -o $img_id/vRhyme_best_bins_fasta_parsed/Nlinked`;	
	`cat $img_id/vRhyme_best_bins_fasta_parsed/Nlinked/*.fasta > $img_id/vRhyme_best_bins_fasta_parsed.Nlinked.fasta`;
	`rm -r $img_id/vRhyme_best_bins_fasta_parsed/Nlinked`;
}

# Step 2. Run checkV for vRhyme_best_bins_fasta_parsed.Nlinked.fasta
open OUT, ">tmp.run_checkV.sh";
foreach my $img_id (sort keys %IMGID){
	`mkdir $img_id/CheckV_phage_bin`;  # Make the output folder for CheckV for phage bins
	print OUT "checkv end_to_end $img_id/vRhyme_best_bins_fasta_parsed.Nlinked.fasta $img_id/CheckV_phage_bin -t 1\n";
}
close OUT;

`cat tmp.run_checkV.sh | parallel -j 20`;

`rm tmp.run_checkV.sh`;
