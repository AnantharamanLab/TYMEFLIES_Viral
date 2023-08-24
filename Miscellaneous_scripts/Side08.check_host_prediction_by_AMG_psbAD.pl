#!/usr/bin/env perl

use strict;
use warnings;

# Aim: Check the host prediction result by searching the host for viruses containing psbAD AMG genes
# Note: psbA: K02703 and psbD: K02706

# Step 1 Store the viral gn to host tax hash
my %Viral_gn2host_tax = (); # $gn => [0] $tax, [1] $method;
open IN, "/storage1/data11/TYMEFLIES_phage/Host_prediction/Viral_gn2host_tax_final.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	$Viral_gn2host_tax{$tmp[0]}[0] = $tmp[1];
	$Viral_gn2host_tax{$tmp[0]}[1] = $tmp[2];
}
close IN;

# Step 2 Store AMG summary hash
my %AMG_summary = (); # $pro => $ko
open IN, "/storage1/data11/TYMEFLIES_phage/AMG_analysis/AMG_summary.txt";
while (<IN>){
	chomp;
	if (!/^Pro/){
		my @tmp = split (/\t/);
		my $pro = $tmp[0];
		my $ko = $tmp[2];
		$AMG_summary{$pro} = $ko;
	}
}
close IN;

# Step 3 Get the psbA to host tax hash
my %PsbA2host_tax = (); # $pro => $host_tax
foreach my $pro (sort keys %AMG_summary){
	my $ko = $AMG_summary{$pro};
	if ($ko eq "K02703"){ # Check if it is psbA
		my $host_tax = "NA\tNA";
		my ($gn) = $pro =~ /^(.+?\_\_.+?)\_\_/;
		if (exists $Viral_gn2host_tax{$gn}){
			$host_tax = $Viral_gn2host_tax{$gn}[0]."\t".$Viral_gn2host_tax{$gn}[1];
		}
		$PsbA2host_tax{$pro} = $host_tax;
	}
}

# Step 4 Write down the result
open OUT, ">PsbA2host_tax.txt";
foreach my $pro (sort keys %PsbA2host_tax){
	print OUT "$pro\t$PsbA2host_tax{$pro}\n";
}
close OUT;

# Step 5 Get the psbD to host tax hash
my %PsbD2host_tax = (); # $pro => $host_tax
foreach my $pro (sort keys %AMG_summary){
	my $ko = $AMG_summary{$pro};
	if ($ko eq "K02706"){ # Check if it is psbD
		my $host_tax = "NA\tNA";
		my ($gn) = $pro =~ /^(.+?\_\_.+?)\_\_/;
		if (exists $Viral_gn2host_tax{$gn}){
			$host_tax = $Viral_gn2host_tax{$gn}[0]."\t".$Viral_gn2host_tax{$gn}[1];
		}
		$PsbD2host_tax{$pro} = $host_tax;
	}
}

# Step 6 Write down the result
open OUT, ">PsbD2host_tax.txt";
foreach my $pro (sort keys %PsbD2host_tax){
	print OUT "$pro\t$PsbD2host_tax{$pro}\n";
}
close OUT;


