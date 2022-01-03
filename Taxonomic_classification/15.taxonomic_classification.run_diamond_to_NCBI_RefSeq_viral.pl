#!/usr/bin/perl

use strict;
use warnings;

# AIM: Use all phage faa files to compare against NCBI Viral RefSeq proteins with diamond

# Step 1. Make tmp file to run diamond in batch
my %Faa = (); # $faa_full_path (full path to the faa file) => $faa
open IN, "find /storage1/data11/TYMEFLIES_phage/33*/vRhyme_best_bins_fasta_parsed -name '*.faa' |";
while (<IN>){
	chomp;
	my $faa_full_path = $_;
	my ($faa) = $faa_full_path =~ /vRhyme_best_bins_fasta_parsed\/(.+?)\.faa/; # Store the name of faa like: 3300020480__vRhyme_unbinned31
	$Faa{$faa_full_path} = $faa;
}
close IN;

# Step 2. Write down the tmp file for running diamond in batch
my $diamond_db = "/slowdata/databases/NCBI_RefSeq_viral/viral.protein.w_tax.dmnd";
`mkdir Taxonomic_classification`;
`mkdir Taxonomic_classification/tmp`;

open OUT, ">tmp.run_diamond_to_NCBI_RefSeq_viral.sh";
foreach my $faa (sort keys %Faa){
	my $faa_name = $Faa{$faa};
	print OUT "diamond blastp -q $faa --db $diamond_db --evalue 0.00001 --query-cover 50 --subject-cover 50 -k 10000 -o Taxonomic_classification/tmp/$faa_name.diamond_out.txt -f 6\n";
}
close OUT;

# Step 3. Run diamond in batch
`cat tmp.run_diamond_to_NCBI_RefSeq_viral.sh | parallel -j 20`;

`rm tmp.run_diamond_to_NCBI_RefSeq_viral.sh`;

# Step 4. Summarize the result
## Step 4.1 Store all faa sequence to faa file map
my %Faa_seq2faa = (); # $faa_seq => $faa
my %Faa_seq_collection = (); # $faa => all the collection of $faa_seq (separated by "\t")
my %Faa2seq_num = (); # Store the number of sequences in each faa; $faa => $faa_seq_num
open IN, "find /storage1/data11/TYMEFLIES_phage/33*/vRhyme_best_bins_fasta_parsed -name '*.faa' |";
while (<IN>){
	chomp;
	my $faa_full_path = $_;
	my ($faa) = $faa_full_path =~ /vRhyme_best_bins_fasta_parsed\/(.+?)\.faa/; # Store the name of $faa like: 3300020480__vRhyme_unbinned31
	
	my $faa_seq_num = 0;
	open INN, "$faa_full_path";
	while (<INN>){
		chomp;
		if (/^>/){
			my ($faa_seq) = $_ =~ /^>(.+?)$/;
			$Faa_seq2faa{$faa_seq} = $faa;
			$faa_seq_num++;
			
			if (!exists $Faa_seq_collection{$faa}){
				$Faa_seq_collection{$faa} = $faa_seq;
			}else{
				$Faa_seq_collection{$faa} .= "\t".$faa_seq;
			}
		}
	}
	close INN;
	
	$Faa2seq_num{$faa} = $faa_seq_num;
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
my %Faa2best_hits = (); # $faa => the collection of best hits (separated by "\t"); Only record this if >= 30% of the proteins for a faa have a hit to Viral RefSeq
open IN, "find Taxonomic_classification/tmp -name '*.diamond_out.txt' |";
while (<IN>){
	chomp;
	my $diamond_out_file = $_;
	my ($faa) = $diamond_out_file =~ /Taxonomic_classification\/tmp\/(.+?)\.diamond/; # Here $faa is actually the phage bin name
	if (-s $diamond_out_file){ # If the $diamond_out_file is not empty
		my %Pro2best_hit = _find_best_hits("$diamond_out_file");
		my $pro_num_w_best_hit = scalar (keys %Pro2best_hit);   # The number of proteins within the $faa have best hits
		my $faa_seq_num = $Faa2seq_num{$faa}; # The total protein number from this $faa
		if ($pro_num_w_best_hit / $faa_seq_num >= 0.3){ # To see if >= 30% of the proteins for a faa have a hit to Viral RefSeq
			my @Best_hits = ();
			foreach my $pro (sort keys %Pro2best_hit){
				my $best_hit = $Pro2best_hit{$pro};
				push @Best_hits, $best_hit;
			}
			my $best_hits = join("\t",@Best_hits);
			$Faa2best_hits{$faa} = $best_hits; 
		}
	}
}
close IN;

### Step 4.2.3 Get the consensus affiliation based on the best hits of individual proteins (>= 50 majority rule)
my %Faa2consensus_tax = (); # $faa => $consensus_tax

foreach my $faa (sort keys %Faa2best_hits){
	my @Best_hits = split (/\t/,$Faa2best_hits{$faa});
	
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
		$Faa2consensus_tax{$faa} = $tax_w_most_frequent; 
	}
}

## Step 4.3 Write down the result
open OUT, ">Taxonomic_classification/Each_bin_consensus_tax_by_NCBI_RefSeq_viral_protein_searching.txt";
foreach my $key (sort keys %Faa2consensus_tax){
	print OUT "$key\t$Faa2consensus_tax{$key}\n";
}
close OUT;

`rm -r Taxonomic_classification/tmp`;



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











