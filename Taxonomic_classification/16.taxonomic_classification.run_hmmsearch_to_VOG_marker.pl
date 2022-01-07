#!/usr/bin/perl

use strict;
use warnings;

# AIM: Use all phage faa files to compare against VOG 587 marker HMM with hmmsearch

# Step 1. Write down the tmp file for running hmmsearch in batch
# Cat all phage genome faa files into one
`find /storage1/data11/TYMEFLIES_phage/33*/vRhyme_best_bins_fasta_parsed -name '*.faa' | xargs cat > All_phage_genome.faa`;

`mkdir Taxonomic_classification/tmp`;

my %VOG_marker2tax = (); # $vog => $tax (only to the family level)
open IN, "/slowdata/databases/VOG209/VOG_marker_table.mdfed.txt";
while (<IN>){
	chomp;
	if (!/^VOG\tFunction/){
		my @tmp = split (/\t/);
		$VOG_marker2tax{$tmp[0]} = $tmp[2];
	}
}
close IN;

open OUT, ">tmp.run_hmmsearch_to_VOG_marker.sh";
foreach my $vog (sort keys %VOG_marker2tax){
	print OUT "hmmsearch -E 0.00001 --cpu 1 --tblout Taxonomic_classification/tmp/$vog.hmmsearch_result.txt /slowdata/databases/VOG209/$vog.hmm All_phage_genome.faa\n";
}
close OUT;

# Step 2. Run hmmsearch in batch
`cat tmp.run_hmmsearch_to_VOG_marker.sh | parallel -j 10`;

`rm tmp.run_hmmsearch_to_VOG_marker.sh`;

`rm All_phage_genome.faa`; # Delete the concatenated phage genome file

# Step 3. Filter hmmsearch result to get protein hits to VOG marker hash
my %Pro2vog = ();
open IN, "cat Taxonomic_classification/tmp/*hmmsearch_result.txt |";
while (<IN>){
	chomp;
	if (!/^#/){
		my $line = $_;
		$line =~ s/ +/ /g;
		my @tmp = split (/\s/,$line);
		my $pro = $tmp[0];
		my $bit_score = $tmp[5];
		my $vog = $tmp[2];
		if ($bit_score >= 40){
			$Pro2vog{$pro} = $vog;
		}
	}
}
close IN;

# Step 4. Find a consensus taxonomy for each bin (simple plurality rule)
my %Phage_gn = (); # $phage_gn => $pro_hits (separeted by "\t"); Store phage genome with pro hits
foreach my $pro (sort keys %Pro2vog){
	my ($phage_gn) = $pro =~ /^(.+?\_\_.+?)\_\_/;
	if (!exists $Phage_gn{$phage_gn}){
		$Phage_gn{$phage_gn} = $pro;
	}else{
		$Phage_gn{$phage_gn} .= "\t".$pro; 
	}
}

my %Phage_gn2consensus_tax = (); # $phage_gn => $consensus_tax (simple plurality rule)
foreach my $phage_gn (sort keys %Phage_gn){
	my @Pro_hits = split(/\t/, $Phage_gn{$phage_gn});
	
	my %Tax_freq = (); # The frequency of each tax
	foreach my $pro (@Pro_hits){
		my $vog = $Pro2vog{$pro}; 
		my $tax = $VOG_marker2tax{$vog}; 
		$Tax_freq{$tax} = 1;
	}
	
	my $consensus_tax = ""; my $consensus_tax_freq = 0;
	foreach my $tax (sort keys %Tax_freq){
		if ($Tax_freq{$tax} > $consensus_tax_freq){
			$consensus_tax_freq = $Tax_freq{$tax};
			$consensus_tax = $tax;
		}
	}
	
	$Phage_gn2consensus_tax{$phage_gn} = $consensus_tax;
}

## Step 5 Write down the result
open OUT, ">Taxonomic_classification/Each_bin_consensus_tax_by_VOG_marker_HMM_searching.txt";
foreach my $key (sort keys %Phage_gn2consensus_tax){
	print OUT "$key\t$Phage_gn2consensus_tax{$key}\n";
}
close OUT;

`rm -r Taxonomic_classification/tmp`;


