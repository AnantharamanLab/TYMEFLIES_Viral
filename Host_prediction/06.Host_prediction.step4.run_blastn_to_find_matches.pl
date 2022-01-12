#!/usr/bin/perl

use strict;
use warnings;

# AIM: Find MAG crispr spacer to phage matches by blastn

# Step 1. Get all the spacers from TYMEFLIES and GEM MAGs
=pod
## Step 1.1 Store all MAG ID in the folder "Host_prediction/find_crispr_out"
my %MAG_ID = (); # Store all the MAG ID; $mag_id => 1; for example: GEM_MAGs_3300013103_44
open IN, "find /storage1/data11/TYMEFLIES_phage/Host_prediction/find_crispr_out -name '*.minced_out.crisprs' |";
while (<IN>){
	chomp;
	my $file = $_;
	my ($mag_id) = $file =~ /find_crispr_out\/(.+?)\.minced_out\.crisprs/;
	$MAG_ID{$mag_id} = 1;
}
close IN;

## Step 1.2 Store %Seq2spacers for MAG in a loop
my %Seq2spacers = (); # Store all spacers to seqs; # $seq => $spacers (combining both MinCED and PILER-CR results)
foreach my $img_id (sort keys %MAG_ID){
	my %Seq2spacers_this_MAG = (); # $seq => $spacers
	my %Seq2spacers_MinCED = ();
	if (-s "/storage1/data11/TYMEFLIES_phage/Host_prediction/find_crispr_out/$img_id.minced_out.crisprs"){
		%Seq2spacers_MinCED = _store_minced_crispr_out("/storage1/data11/TYMEFLIES_phage/Host_prediction/find_crispr_out/$img_id.minced_out.crisprs");
	}
	
	my %Seq2spacers_PILER_CR = ();
	`sed -n '/DETAIL REPORT/,/SUMMARY BY SIMILARITY/p' /storage1/data11/TYMEFLIES_phage/Host_prediction/find_crispr_out/$img_id.pilercr_out.crisprs > /storage1/data11/TYMEFLIES_phage/Host_prediction/find_crispr_out/$img_id.pilercr_out.crisprs.parsed`;
	if (-s "/storage1/data11/TYMEFLIES_phage/Host_prediction/find_crispr_out/$img_id.pilercr_out.crisprs.parsed"){
		%Seq2spacers_PILER_CR = _store_pilercr_crispr_out("/storage1/data11/TYMEFLIES_phage/Host_prediction/find_crispr_out/$img_id.pilercr_out.crisprs.parsed");
	}
	`rm /storage1/data11/TYMEFLIES_phage/Host_prediction/find_crispr_out/$img_id.pilercr_out.crisprs.parsed`;
	
	if (!(%Seq2spacers_MinCED) and !(%Seq2spacers_PILER_CR)){
		%Seq2spacers_this_MAG = ();
	}elsif (!(%Seq2spacers_MinCED) and %Seq2spacers_PILER_CR){
		%Seq2spacers_this_MAG =  %Seq2spacers_PILER_CR;
	}elsif (%Seq2spacers_MinCED and !(%Seq2spacers_PILER_CR)){
		%Seq2spacers_this_MAG = %Seq2spacers_MinCED;
	}elsif (%Seq2spacers_MinCED and %Seq2spacers_PILER_CR){
		%Seq2spacers_this_MAG = %Seq2spacers_MinCED;
		foreach my $seq (sort keys %Seq2spacers_this_MAG){
			if (exists $Seq2spacers_PILER_CR{$seq}){
				my $spacers_MinCED = $Seq2spacers_MinCED{$seq};
				my $spacers_PILER_CR = $Seq2spacers_PILER_CR{$seq};
				my $spacers_combined = $spacers_MinCED."\,".$spacers_PILER_CR;
				my @Spacers_combined = split (/\,/,$spacers_combined);
				my %Spacers_combined_drep = ();
				foreach my $key (@Spacers_combined){
					$Spacers_combined_drep{$key} = 1;
				}				
				$spacers_combined = join("\,", (keys %Spacers_combined_drep)); # dRep-ed spacers
				$Seq2spacers_this_MAG{$seq} = $spacers_combined;
			}
		}
		
		foreach my $seq (sort keys %Seq2spacers_PILER_CR){
			if (!exists $Seq2spacers_this_MAG{$seq}){
				$Seq2spacers_this_MAG{$seq} = $Seq2spacers_PILER_CR{$seq};
			}
		}
	}
	
	%Seq2spacers = (%Seq2spacers, %Seq2spacers_this_MAG);
}

## Step 1.3 Write down all the spacers into a fasta files
my %All_spacers = (); # $head => $spacer_seq
foreach my $seq (sort keys %Seq2spacers){
	my $spacers = $Seq2spacers{$seq};
	my @Spacers = split (/\,/,$spacers);
	for(my $i=0; $i<=$#Spacers; $i++){
		my $j = $i+1;
		$j = (sprintf "%04d", $j);
		my $head = ">$seq\_\_Spacer_$j";
		$All_spacers{$head} = $Spacers[$i];
	}
}

open OUT, ">/storage1/data11/TYMEFLIES_phage/Host_prediction/All_spacers.fasta";
foreach my $key (sort keys %All_spacers){
	print OUT "$key\n$All_spacers{$key}\n";
}
close OUT;

# Step 2. Make blastn db 

`makeblastdb -in Host_prediction/All_spacers.fasta -title All_spacers -dbtype nucl -out Host_prediction/All_spacers_blastndb`;

# Step 3. Run blastn to find spacer and phage association
`find /storage1/data11/TYMEFLIES_phage/33*/vRhyme_best_bins_fasta_parsed -name '*.fasta' -exec cat {} + > /storage1/data11/TYMEFLIES_phage/Host_prediction/All_phage_genomes.fasta`;
=cut
`mkdir Host_prediction/All_phage_genomes_split_fsa`;
`perl /storage1/data11/TYMEFLIES_phage/split_multifasta.pl --in Host_prediction/All_phage_genomes.fasta --output_dir=Host_prediction/All_phage_genomes_split_fsa --seqs_per_file=10000`;

open OUT, ">tmp.run_blastn_to_find_matches_in_parallel.sh";
open IN, "ls Host_prediction/All_phage_genomes_split_fsa/*.fsa |";
while (<IN>){
	chomp;
	my $fsa = $_;
	my ($fsa_name) = $fsa =~ /All_phage_genomes_split_fsa\/(.+?)\.fsa/;
	print OUT "blastn -task blastn-short -evalue 1 -dust no -word_size 7 -num_threads 1 -outfmt 6 -query $fsa -db Host_prediction/All_spacers_blastndb -out Host_prediction/all_phage_genome.${fsa_name}.2all_spacer.blastn_out.txt\n";
}
close IN;
close OUT;

`cat tmp.run_blastn_to_find_matches_in_parallel.sh | parallel -j 10`;
`rm tmp.run_blastn_to_find_matches_in_parallel.sh`;

=pod
# Step 4. Parse blastn results to store matches
# There are two categories of matches: 
# (1) have 0 or 1 mismatch over the entire spacer length (‘CRISPR (near)identical’) 
# (2) have ≥80% identity over the entire spacer length (‘CRISPR multiple partial’)

## Step 4.1 Store the spacer length
my %Spacer2length = (); # $spacer => $length; for example: 2004178001.a:gws2_d1_0103_30__Spacer_0001 => 32 
my %Hash = _store_seq("Host_prediction/All_spacers.fasta");
foreach my $key (sort keys %Hash){
	my ($key_clean) = $key =~ /^>(.+?)$/;
	my $length = length($Hash{$key});
	$Spacer2length{$key_clean} = $length;
}

## Step 4.2 Parse blastn results
my %CRISPR_near_identical_matches = (); # $match => 1; an example for a match: 2004178001.a:gws2_d1_0103_30__Spacer_0006|3300044846__vRhyme_449__Ga0453141_0000520
my %CRISPR_partial_matches = (); # $match => 1;
open IN, "Host_prediction/All_spacer2all_phage_genome.blastn_out.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $spacer = $tmp[0];
	my $phage_scf = $tmp[1];
	my $spacer_length = $Spacer2length{$spacer};
	
	my $identity = $tmp[2];
	my $mismatches = $tmp[4];
	my $gap_openings = $tmp[5];
	
	my $q_start = $tmp[6];
	my $q_end = $tmp[7];
	my $spacer_matched_length = abs($q_end - $q_start);
	
	if ($spacer_matched_length eq $spacer_length){
		if ($mismatches <= 1 and $gap_openings == 0){
			$CRISPR_near_identical_matches{"$spacer\|$phage_scf"} = 1;
		}elsif ($identity >= 80){
			$CRISPR_partial_matches{"$spacer\|$phage_scf"} = 1;
		}
	}
}
close IN;

my $num_CRISPR_near_identical_matches = scalar (keys (%CRISPR_near_identical_matches));
print "There are $num_CRISPR_near_identical_matches CRISPR (near)identical matches\n";
my $num_CRISPR_partial_matches = scalar (keys (%CRISPR_partial_matches));
print "There are $num_CRISPR_partial_matches CRISPR partial matches\n";

# Subroutine

sub _store_minced_crispr_out{
	my $file = $_[0];
	my %Hash = (); # $spacer => $seq
	my %Seq = (); # $seq => 1
	my $seq = "";
	open IN_, $file;
	while (<IN_>){
		chomp;
		if (/^Sequence/){
			my $line = $_;
			($seq) = $line =~ /\'(.+?)\'/;
			$Seq{$seq} = 1;
		}elsif(/\[/){
			my $line = $_;
			my @tmp = split (/\t/,$line);
			my $spacer = $tmp[3];
			$Hash{$spacer} = $seq;
		}
	}
	close IN_;
	
	my %Hash2 = (); # $seq => $spacer collection separated by ","
	foreach my $spacer (sort keys %Hash){
		my $seq = $Hash{$spacer};
		if (!exists $Hash2{$seq}){
			$Hash2{$seq} = $spacer;
		}else{
			$Hash2{$seq} .= "\,".$spacer;
		}
	}
	
	return %Hash2; 
}

sub _store_pilercr_crispr_out{
	my $file = $_[0];
	my %Hash = (); # $spacer => $seq
	my %Seq = (); # $seq => 1
	my $seq = "";
	open IN_, $file;
	while (<IN_>){
		chomp;
		if (/^>/){
			my $line = $_;
			($seq) = $line =~ /^>(.+?)$/;
			$Seq{$seq} = 1;
		}elsif(/\.\./){
			my $line = $_;
			$line =~ s/ +/ /g; # Change multiple spaces into one space
			$line =~ s/^ //g; # Delete the first space in the front
			my @tmp = split (/\s/,$line);
			my $spacer_length = $tmp[3];
			if ($spacer_length =~ /\d\d/){ # Get to the spacer line
				my $spacer = $tmp[-1];
				$Hash{$spacer} = $seq;				
			}						
		}
	}
	close IN_;
	
	
	my %Hash2 = (); # $seq => $spacer collection separated by ","
	foreach my $spacer (sort keys %Hash){
		my $seq = $Hash{$spacer};
		if (!exists $Hash2{$seq}){
			$Hash2{$seq} = $spacer;
		}else{
			$Hash2{$seq} .= "\,".$spacer;
		}
	}
	
	return %Hash2; 
}

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