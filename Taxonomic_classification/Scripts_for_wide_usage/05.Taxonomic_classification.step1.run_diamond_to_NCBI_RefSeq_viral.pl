#!/usr/bin/perl

use strict;
use warnings;

# AIM: Use all phage faa files to compare against NCBI Viral RefSeq proteins with diamond, and get virus taxonomy
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

#        perl 05.Taxonomic_classification.step1.run_diamond_to_NCBI_RefSeq_viral.pl [virus_protein_sequence.faa] [map_file.txt] 
#        The first two items are the inputs

my $input_faa = $ARGV[0];
my $map = $ARGV[1];
my $input_faa_name = "";
if ($input_faa !~ /\//){
	($input_faa_name) = $input_faa =~ /^(.+?)\.faa/;
}else{
	($input_faa_name) = $input_faa =~ /^.+\/(.+?)\.faa/;
}

# Step 1. Run diamond
my $diamond_db = "/slowdata/databases/NCBI_RefSeq_viral/viral.protein.w_tax.dmnd";
`mkdir Taxonomic_classification`;
`mkdir Taxonomic_classification/tmp`;

## Run diamond with 10 cpus

`diamond blastp -q $input_faa -p 10 --db $diamond_db --evalue 0.00001 --query-cover 50 --subject-cover 50 -k 10000 -o Taxonomic_classification/tmp/$input_faa_name.diamond_out.txt -f 6 -b 16`;

# Step 2. Summarize the result
## Step 4.1 Store all faa sequence to gn file map
my %Faa_seq2gn = (); # $faa_seq => $gn 
my %Gn2all_seq = (); # $gn  => all the collection of $faa_seq (separated by "\t")
my %Gn2seq_num = (); # Store the number of sequences in each genome; $gn => $faa_seq_num
open IN, "$map";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $faa_seq = $tmp[0];
	my $gn = $tmp[1];
	$Faa_seq2gn{$faa_seq} = $gn;
	if (! exists $Gn2all_seq{$gn}){
		$Gn2all_seq{$gn} = $faa_seq;
	}else{
		$Gn2all_seq{$gn} .= "\t".$faa_seq;
	}
	$Gn2seq_num{$gn}++;
}
close IN;

## Step 4.2 Store the best hit and best hit taxonomy of each faa seq
### Step 4.2.1 Store the NCBI_RefSeq_viral protein to tax hash (ICTV tax with 8 ranks)
my %NCBI_RefSeq_viral_protein2tax = (); # $pro => $tax;
open IN, "/slowdata/databases/NCBI_RefSeq_viral/viral.protein.ictv_8_rank_tax.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $pro = $tmp[0];
	my $tax = $tmp[1];
	$NCBI_RefSeq_viral_protein2tax{$pro} = $tax;
}
close IN;

### Step 4.2.2 Store the best hits and to see whether >= 30% of the proteins for a faa have a hit to Viral RefSeq
my %Gn2best_hits = (); # $gn => the collection of best hits (separated by "\t"); Only record this if >= 30% of the proteins for a genome have a hit to Viral RefSeq

my $diamond_out_file = "Taxonomic_classification/tmp/$input_faa_name.diamond_out.txt";
if (-s $diamond_out_file){ # If the $diamond_out_file is not empty
	my %All_Pro2best_hit = _find_best_hits("$diamond_out_file"); # This is the result for all input faa
	my %Gn_involved = (); # $gn => 1; Store the $gn that have sequences inside have the best hit
	foreach my $pro (sort keys %All_Pro2best_hit){
		my $gn = $Faa_seq2gn{$pro};
		$Gn_involved{$gn} = 1;
	}
	# Split the best hit result into each bin ($gn)
	foreach my $gn (sort keys %Gn_involved){
		my %Pro2best_hit = (); # Store all the proteins that have best hits in this $gn
			
		foreach my $pro (sort keys %All_Pro2best_hit){
			if ($Faa_seq2gn{$pro} eq $gn){
				$Pro2best_hit{$pro} = $All_Pro2best_hit{$pro};
			}
		}
			
		my $pro_num_w_best_hit = scalar (keys %Pro2best_hit);   # The number of proteins within the $gn have best hits
		my $faa_seq_num = $Gn2seq_num{$gn}; # The total protein number from this $gn
		if ($pro_num_w_best_hit / $faa_seq_num >= 0.3){ # To see if >= 30% of the proteins for a genome have a hit to Viral RefSeq
			my @Best_hits = ();
			foreach my $pro (sort keys %Pro2best_hit){
				my $best_hit = $Pro2best_hit{$pro};
				push @Best_hits, $best_hit;
			}
			my $best_hits = join("\t",@Best_hits);
			$Gn2best_hits{$gn} = $best_hits; 
		}
	}
}

### Step 4.2.3 Get the consensus affiliation based on the best hits of individual proteins (>= 50 majority rule)
my %Gn2consensus_tax = (); # $gn => $consensus_tax

foreach my $gn (sort keys %Gn2best_hits){
	my @Best_hits = split (/\t/,$Gn2best_hits{$gn});
	
	my %Tax2freq = (); # The frequency of tax; $tax => $frequency
	foreach my $best_hit (@Best_hits){
		my $tax = $NCBI_RefSeq_viral_protein2tax{$best_hit};
		$Tax2freq{$tax}++;
	}
	
	my $tax_w_most_frequent = ""; my $highest_freq = 0;
	foreach my $tax (sort keys %Tax2freq){
		my $freq = $Tax2freq{$tax};
		if ($freq > $highest_freq){
			$tax_w_most_frequent = $tax;
			$highest_freq = $freq;
		}
	}
	
	my $perc_of_tax_w_most_frequent = 0;
	$perc_of_tax_w_most_frequent = $highest_freq / (scalar @Best_hits);
	
	if ($perc_of_tax_w_most_frequent >= 0.5){ # Only store the consensus tax if there is
		$Gn2consensus_tax{$gn} = $tax_w_most_frequent; 
	}
}

## Step 4.3 Write down the result
open OUT, ">Taxonomic_classification/Each_viral_genome_consensus_tax_by_NCBI_RefSeq_viral_protein_searching.txt";
foreach my $key (sort keys %Gn2consensus_tax){
	print OUT "$key\t$Gn2consensus_tax{$key}\n";
}
close OUT;

#`rm -r Taxonomic_classification/tmp`;



# Subroutine

sub _find_best_hits{ # Feasible even if the diamond out file is not ordered by bit score 
	my $file = $_[0]; # The input diamond out file
	my %Pro2best_hit = (); # $pro => [0] $best_hit (for example; 3300020575__vRhyme_28__Ga0208053_1000517_4 => YP_009226243.1)
						   # [1] $bit_score 
	
	open IN_, "$file";
	while (<IN_>){
		chomp;
		my @tmp = split (/\t/);
		my $pro = $tmp[0];
		my $hit = $tmp[1];
		my $bit_score = $tmp[-1];
		
		if (! exists $Pro2best_hit{$pro}[0]){
			$Pro2best_hit{$pro}[0] = $hit; # Store the best hit (temporary)
			$Pro2best_hit{$pro}[1] = $bit_score; # Store the bit score (temporary)
		}else{
			if ($bit_score > $Pro2best_hit{$pro}[1]){ # If the bit score of a second hit > the temporary best hit's bit score
				$Pro2best_hit{$pro}[0] = $hit;
				$Pro2best_hit{$pro}[1] = $bit_score;
			}
		}
	}
	close IN_;
	
	my %Result = (); # $pro => $best_hit
	foreach my $pro (sort keys %Pro2best_hit){
		my $best_hit = $Pro2best_hit{$pro}[0];
		$Result{$pro} = $best_hit;
	}
	
	return %Result;
}











