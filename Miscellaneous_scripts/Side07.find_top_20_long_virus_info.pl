#!/usr/bin/perl

use strict;
use warnings;

# AIM: Get the length, bin member number, tax for the top 20 longest viruses

my $n_top = 20; # Get the top 20 longest virues

# Step 1 Get the bin member number hash
my %Bin2bin_member_num = (); # $bin => $bin_member_num
open IN, "cat /storage1/data11/TYMEFLIES_phage/33*/vRhyme_best_bins_fasta_parsed/Each_bin_info.txt | sort -k 2 -n | grep -v 'un' | grep -v 'bin' | ";
while (<IN>){
	chomp;
	my $line = $_;
	my @tmp = split (/\t/,$line);
	my $bin = $tmp[0];
	my $bin_member_num = $tmp[1];
	$Bin2bin_member_num{$bin} = $bin_member_num;
}
close IN;

# Step 2 Get the bin length
## Step 2.1 Get all scaffolds length and completeness information
my %Scf2length_n_completeness = (); # $scf => [0] $length [1] $completeness
open IN, "ls /storage1/data11/TYMEFLIES_phage/33*/CheckV_phage_scaffold/quality_summary.tsv |";
while (<IN>){
	chomp;
	my $file = $_;
	open INN, "$file";
	while (<INN>){
		chomp;
		if (!/^contig_id/){
			my @tmp = split (/\t/);
			my $scf = $tmp[0];
			my $length = $tmp[1];
			my $completeness = $tmp[9];
			$Scf2length_n_completeness{$scf}[0] = $length;
			$Scf2length_n_completeness{$scf}[1] = $completeness;
		}
	}
	close INN;
}
close IN;

# Step 2.2 Get scf 2 bin information
my %Scf_full2bin = (); # $scf_full => $bin
my %Bin = (); # $bin => $scf_full collection, separated by "\t"
my %Scf_full2scf = (); # $scf_full => $scf
my %Scf2scf_full = (); # $scf => $scf_full

open IN, "find /storage1/data11/TYMEFLIES_phage/33* -maxdepth 1 -name 'vRhyme_best_bins_fasta_parsed' -type d |";
while (<IN>){
	chomp;
	my $folder_adr = $_;
	open INN, "grep '>' $folder_adr/*.fasta |";
	while (<INN>){
		chomp; # for example, a line looks like: 
	    # /storage1/data11/TYMEFLIES_phage/3300044855/vRhyme_best_bins_fasta_parsed/3300044855__vRhyme_100.fasta:>3300044855__vRhyme_100__Ga0453711_005625
		my ($bin) = $_ =~ /vRhyme_best_bins_fasta_parsed\/(.+?)\.fasta\:\>/;
		my ($scf_full) = $_ =~ /\:\>(.+?)$/;
		my ($scf) = $scf_full =~ /^.+?\_\_.+?\_\_(.+?)$/;
		
		$Scf_full2bin{$scf_full} = $bin;
		$Scf_full2scf{$scf_full} = $scf;
		
		if (!exists $Bin{$bin}){
			$Bin{$bin} = $scf_full;
		}else{
			$Bin{$bin} .= "\t".$scf_full;
		}
		
		$Scf2scf_full{$scf} = $scf_full;
	}
	close INN;
}
close IN;

## Step 2.3 Store bin length info to %Bin2length
my %Bin2length = (); # $bin => $length
foreach my $bin (sort keys %Bin){
	my $length = 0;
	
	my @Scf_full_collection = split (/\t/, $Bin{$bin});
	foreach my $scf_full (@Scf_full_collection){
		my $scf = $Scf_full2scf{$scf_full};
		my $scf_length = $Scf2length_n_completeness{$scf}[0];
		$length += $scf_length;
	}

	$Bin2length{$bin} = $length;
}

# Step 3 Store bin to tax
my %Bin2tax = (); # $tax => [0] tax (based on NCBI RefSeq viral protein database) [1] tax (based on VOG marker HMMs)
open IN, "/storage1/data11/TYMEFLIES_phage/Taxonomic_classification/Each_bin_consensus_tax_by_NCBI_RefSeq_viral_protein_searching.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $bin = $tmp[0];
	my $tax = $tmp[1];
	$Bin2tax{$bin}[0] = $tax;
}
close IN;

open IN, "/storage1/data11/TYMEFLIES_phage/Taxonomic_classification/Each_bin_consensus_tax_by_VOG_marker_HMM_searching.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $bin = $tmp[0];
	my $tax = $tmp[1];
	$Bin2tax{$bin}[1] = $tax;
}
close IN;

# Step 4 Store the top n (here, 20) virues info
my %Top_viruses_info = (); # $bin => [0] $length [1] $bin_member_num [2] $tax (based on NCBI RefSeq viral protein database) [3] $tax (based on VOG marker HMMs)

my @Top_viruses = (); 
foreach my $bin (reverse sort {$Bin2length{$a} <=> $Bin2length{$b}} keys %Bin2length){
	push @Top_viruses, $bin;
}

splice @Top_viruses, $n_top; # Get the first $n_top elements

foreach my $bin (@Top_viruses){
	my $length = $Bin2length{$bin};
	my $bin_member_num = $Bin2bin_member_num{$bin};
	my $tax1 = "NA";
	if (exists $Bin2tax{$bin}[0]){
		$tax1 = $Bin2tax{$bin}[0];
	}
	my $tax2 = "NA";
	if (exists $Bin2tax{$bin}[1]){
		$tax2 = $Bin2tax{$bin}[1];
	}	
	$Top_viruses_info{$bin}[0] = $length;
	$Top_viruses_info{$bin}[1] = $bin_member_num;
	$Top_viruses_info{$bin}[2] = $tax1;
	$Top_viruses_info{$bin}[3] = $tax2;
}

# Step 5 Write down the result
open OUT, ">Top_${n_top}_long_virues_info.txt";
print OUT "Bin\tLength\tBin member num\tTax (based on NCBI RefSeq viral protein database)\tTax (based on VOG marker HMMs)\n";
foreach my $bin (@Top_viruses){
	print OUT "$bin\t$Top_viruses_info{$bin}[0]\t$Top_viruses_info{$bin}[1]\t$Top_viruses_info{$bin}[2]\t$Top_viruses_info{$bin}[3]\n";
}
close OUT;


