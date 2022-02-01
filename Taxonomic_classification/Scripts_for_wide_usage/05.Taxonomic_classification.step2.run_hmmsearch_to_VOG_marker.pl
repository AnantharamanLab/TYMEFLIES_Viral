#!/usr/bin/perl

use strict;
use warnings;

# AIM: Use all phage faa files to compare against VOG 587 marker HMM with hmmsearch, and get virus taxonomy
# Usage: The input should be two files:
#        1) all virus protein sequences (file ended with "faa")
#        2) a map file of protein to genome; this map file should contain two columns: 1. each protein 2. the corresponding genome of this protein
#        An example here:
#        protein_id	genome_id
#		 protein_A	genome_A
#		 protein_B	genome_A
#		 protein_C	genome_A
#		 protein_D	genome_B
#		 protein_E	genome_B
#		 protein_F	genome_B	

#        perl 05.Taxonomic_classification.step2.run_hmmsearch_to_VOG_marker.pl [virus_protein_sequence.faa] [map_file.txt] 
#        The first two items are the inputs

my $input_faa = $ARGV[0];
my $map = $ARGV[1];
my $input_faa_name = "";
if ($input_faa !~ /\//){
	($input_faa_name) = $input_faa =~ /^(.+?)\.faa/;
}else{
	($input_faa_name) = $input_faa =~ /^.+\/(.+?)\.faa/;
}

# Step 1. Write down the tmp file for running hmmsearch in batch
# Cat all phage genome faa files into one
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
	print OUT "hmmsearch -E 0.00001 --cpu 1 --tblout Taxonomic_classification/tmp/$vog.hmmsearch_result.txt /slowdata/databases/VOG209/$vog.hmm $input_faa\n";
}
close OUT;

# Step 2. Run hmmsearch in batch
`cat tmp.run_hmmsearch_to_VOG_marker.sh | parallel -j 10`;

`rm tmp.run_hmmsearch_to_VOG_marker.sh`;

`rm $input_faa`; # Delete the concatenated phage genome file

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
open IN, "$map";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $pro = $tmp[0];
	my $gn = $tmp[1];
	if (exists $Pro2vog{$pro}){ # If this protein has a VOG hit
		if (! exists $Phage_gn{$gn}){
			$Phage_gn{$gn} = $pro;
		}else{
			$Phage_gn{$gn} .= "\t".$pro;
		}
	}
}
close IN;

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


