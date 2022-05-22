#!/usr/bin/perl

use strict;
use warnings;

# Aim: Grep gene files (in prodigal format) for "All_phage_species_rep_gn_containing_AMG.fasta"

# Step 1 Store viruses species rep genomes that contain AMG
my %Viruses_species_rep_gn = (); # $gn => 1
open IN, "All_phage_species_rep_gn_containing_AMG.fasta";
while (<IN>){
	chomp;
	if (/^>/){
		my ($gn) = $_ =~ /^>(.+?\_\_.+?)\_\_/;
		$Viruses_species_rep_gn{$gn} = 1;
	}
}
close IN;

# Step 2 Store the ffn file (gene file) 
## Step 2.1 store all IMG ID
my %IMGID = (); # $img_id => 1
open IN, "ls -l | sed -r s'/ +/\t/g' | grep ^d | cut -f 9 | grep '^33' |";
while (<IN>){
	chomp;
	my $img_id = $_;
	$IMGID{$img_id} = 1;
}
close IN;

## Step 2.2 store the ffn file (only store genes from Viruses_species_rep_gn)
my %All_ffn_seq = (); # Store all the ffn sequences that belong to %Viruses_species_rep_gn
#`find /storage1/data11/TYMEFLIES_phage/*/vRhyme_best_bins_fasta_parsed -name '*.ffn' -exec cat {} + > TYMEFLIES_All_ffn_seq.ffn`;
%All_ffn_seq = _store_seq("TYMEFLIES_All_ffn_seq.ffn");
#`rm TYMEFLIES_All_ffn_seq.ffn`;
foreach my $key (sort keys %All_ffn_seq){
	my ($gn) = $key =~ /^>(.+?\_\_.+?)\_\_/;
	if (!exists $Viruses_species_rep_gn{$gn}){
		delete $All_ffn_seq{$key};
	}
}

# Step 3 Store all the prodigal ffn file headers of all metagenomes (only the gene on viral scaffolds)
my %All_metagenome_ffn_seq_headers = (); # $clean_gene_id => $header; 
foreach my $img_id (sort keys %IMGID){
	my %Viral_scf = (); # Store the viral scaffolds in this metagenome (excluding fragment_X)
	open IN, "$img_id/VIBRANT_$img_id.a/VIBRANT_phages_$img_id.a/$img_id.a.phages_combined.fna";
	while (<IN>){
		chomp;
		if (/^>/){
			my $header = $_;
			if ($header !~ /fragment/){
				my ($viral_scf) = $header =~ /^>(.+?)$/;
				$Viral_scf{$viral_scf} = 1;
			}else{
				my ($viral_scf) = $header =~ /^>(.+?)_fragment/;
				$Viral_scf{$viral_scf} = 1;
			}
		}
	}
	close IN;
	
	open IN, "$img_id/VIBRANT_$img_id.a/$img_id.a.prodigal.ffn";
	while (<IN>){
		chomp;
		if (/^>/){
			my $header = $_;
			my ($clean_gene_id) = $header =~ /^>(.+?)\s/;
			my ($scf) = $clean_gene_id =~ /^(.+)\_/;
			if (exists $Viral_scf{$scf}){
				$All_metagenome_ffn_seq_headers{$clean_gene_id} = $header;
			}
		}
	}
	close IN;
}

# Step 4 Make a new hash - %All_ffn_seq2
my %All_ffn_seq2 = ();
foreach my $header1 (sort keys %All_ffn_seq){
	my $header1_transferred = ""; # Transfer back the header (gene id) to the original one in the metagenome
	if ($header1 !~ /fragment/){
		($header1_transferred) = $header1 =~ /^>33.+?\_\_.+?\_\_(.+?)$/;
	}else{
		my $header1_transferred_1 = "";my $header1_transferred_2 = "";
		($header1_transferred_1, $header1_transferred_2) = $header1 =~ /^>33.+?\_\_.+?\_\_(.+?)\_fragment\_\d+?(\_\d+?)$/;
		$header1_transferred = $header1_transferred_1.$header1_transferred_2;
	}
	
	my $header2 = "";
	if (exists $All_metagenome_ffn_seq_headers{$header1_transferred}){
		$header2 = $All_metagenome_ffn_seq_headers{$header1_transferred};
	}
	
	my ($header1_clean) = $header1 =~ /^>(.+?)$/;

	$header2 =~ s/$header1_transferred/$header1_clean/g;
	
	$All_ffn_seq2{$header2} = $All_ffn_seq{$header1};
}


## Step 5 Write down the ffn sequences (gene sequences)
open OUT, ">All_phage_species_rep_gn_containing_AMG.genes";
foreach my $key (sort keys %All_ffn_seq2){
	print OUT "$key\n$All_ffn_seq2{$key}\n";
}
close OUT;



## Subroutine

sub _store_seq{ 
	my $file = $_[0];
	my %Seq = (); my $head = "";
	open _IN, "$file";
	while (<_IN>){
		chomp;
		if (/>/){
			if (/\s/){
				($head) = $_ =~ /^(>.+?)\s/;
				$Seq{$head} = "";
			}else{
				($head) = $_ =~ /^(>.+?)$/;
				$Seq{$head} = "";
			}
		}else{
			$Seq{$head} .= $_;
		}
	}
	close _IN;
	return %Seq;
}