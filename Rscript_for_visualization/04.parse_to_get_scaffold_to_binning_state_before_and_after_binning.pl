#!/usr/bin/perl

use strict;
use warnings;

# AIM: Get scaffold to binned or unbinned state before and after binning

# Step 1. Get scf 2 bin information
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

# Step 2. Write the result
# Before binning => unbinned => Scaffold percentage
# After binning => binned or unbinned => Scaffold percentage

my $total_scaffold_no = 0;
my $unbinned_before_binning = 0;
my $binned_before_binning = 0;
my $binned_after_binning = 0;
my $unbinned_after_binning = 0;

foreach my $scf_full (sort keys %Scf_full2scf){
	if ($scf_full =~ /unbinned/){
		$unbinned_after_binning++;
	}else{
		$binned_after_binning++;
	}
	$total_scaffold_no++;
}

$unbinned_before_binning = $total_scaffold_no;

my $unbinned_before_binning_percentage = $unbinned_before_binning / $total_scaffold_no * 100;
my $binned_before_binning_percentage = $binned_before_binning / $total_scaffold_no * 100;
my $binned_after_binning_percentage = $binned_after_binning / $total_scaffold_no * 100;
my $unbinned_after_binning_percentage = $unbinned_after_binning / $total_scaffold_no * 100;

open OUT, ">Scaffold_to_binning_state_before_and_after_binning.txt";
print OUT "Before binning\tunbinned\t$unbinned_before_binning_percentage\n";
print OUT "Before binning\tbinned\t$binned_before_binning_percentage\n";
print OUT "After binning\tunbinned\t$binned_after_binning_percentage\n";
print OUT "After binning\tbinned\t$unbinned_after_binning_percentage\n";
close OUT;


