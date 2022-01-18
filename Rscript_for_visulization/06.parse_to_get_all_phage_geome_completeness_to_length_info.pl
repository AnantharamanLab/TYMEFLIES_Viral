#!/usr/bin/perl

use strict;
use warnings;

# AIM: Get all phage genome completeness and length info for plotting length boxplot for each completeness category

# Step 1. Get all scaffolds length and completeness information
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

# Step 2. Get scf 2 bin information
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

# Step 3. Store CheckV result for bin completeness
my %Bin2completeness = (); # $bin => [0] $checkv_quality [1]$completeness
open IN, "/storage1/data11/TYMEFLIES_phage/CheckV_phage_bin_all/quality_summary.tsv";
while (<IN>){
	chomp;
	if (!/^contig_id/){
		my @tmp = split (/\t/);
		my $bin = $tmp[0];
		my $checkv_quality = $tmp[7];
		my $completeness = $tmp[9];
		$Bin2completeness{$bin}[0] = $checkv_quality;
		$Bin2completeness{$bin}[1] = $completeness;
	}
}
close IN;

# Step 4. Get checkv_quality to completeness range hash and print to check manually
my %Checkv_quality2completeness_values = (); 
my %Checkv_quality2completeness_range = (); 
foreach my $bin (sort keys %Bin2completeness){
	my $checkv_quality = $Bin2completeness{$bin}[0];
	my $completeness = $Bin2completeness{$bin}[1];
	if (! exists $Checkv_quality2completeness_values{$checkv_quality}){
		$Checkv_quality2completeness_values{$checkv_quality} = $completeness;
	}else{
		$Checkv_quality2completeness_values{$checkv_quality} .= "\t".$completeness;
	}
}

foreach my $checkv_quality (sort keys %Checkv_quality2completeness_values){
	my @Values = split (/\t/,$Checkv_quality2completeness_values{$checkv_quality});
	if ($Checkv_quality2completeness_values{$checkv_quality} =~ /NA|completeness/){
		@Values = sort { $a cmp $b } @Values; # Sort the array
	}else{
		@Values = sort { $a <=> $b } @Values; # Sort the array
	}
	my $range = $Values[0]."\-".$Values[-1];
	$Checkv_quality2completeness_range{$checkv_quality} = $range;
}

foreach my $checkv_quality (sort keys %Checkv_quality2completeness_range){
	my $range = $Checkv_quality2completeness_range{$checkv_quality};
	print "$checkv_quality\t$range\n";
}

# Step 5. Store bin length info to %Bin2length
my %Bin2length = (); # $bin => $length
foreach my $bin (sort keys %Bin2completeness){
	my $length = 0;
	
	my @Scf_full_collection = split (/\t/, $Bin{$bin});
	foreach my $scf_full (@Scf_full_collection){
		my $scf = $Scf_full2scf{$scf_full};
		my $scf_length = $Scf2length_n_completeness{$scf}[0];
		$length += $scf_length;
	}

	$Bin2length{$bin} = $length;
}

# Step 6. Write down "Bin2checkv_quality_and_length.txt"
open OUT, ">Bin2checkv_quality_and_length.txt";
foreach my $bin (sort keys %Bin2completeness){
	my $checkv_quality = $Bin2completeness{$bin}[0];
	my $length = $Bin2length{$bin};
	print OUT "$bin\t$checkv_quality\t$length\n";
}
close OUT;



