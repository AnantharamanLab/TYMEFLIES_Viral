#!/usr/bin/perl

use strict;
use warnings;

# Aim: Compare the genus cluster to family cluster

# Step 1 Store genus cluster
my %Gn2genus_cluster = (); # $gn => $cluster
my %Genus_cluster2gn = (); # $cluster => $gn collection separated by "\t"

my $i = 1; # The number of cluster
open IN, "/storage1/data11/TYMEFLIES_phage/Cluster_phage_genomes/genus_clusters.txt";
while (<IN>){
	chomp;
	my $line = $_;
	my $cluster = "Genus_cluster${i}";
	$Genus_cluster2gn{$cluster} = $line;
	
	my @tmp = split (/\t/, $line);
	foreach my $gn (@tmp){
		$Gn2genus_cluster{$gn} = $cluster;
	}
	$i++;
}
close IN;

# Step 2 Store family cluster
my %Gn2family_cluster = (); # $gn => $cluster
my %Family_cluster2gn = (); # $cluster => $gn collection separated by "\t"

my $j = 1; # The number of cluster
open IN, "/storage1/data11/TYMEFLIES_phage/Cluster_phage_genomes/family_clusters.txt";
while (<IN>){
	chomp;
	my $line = $_;
	my $cluster = "Family_cluster${j}";
	$Family_cluster2gn{$cluster} = $line;
	
	my @tmp = split (/\t/, $line);
	foreach my $gn (@tmp){
		$Gn2family_cluster{$gn} = $cluster;
	}
	$j++;
}
close IN;

# Step 3 Find singleton genus cluster
my %Genus_cluster_singleton = (); # $cluster (singleton genus cluster) => 1
foreach my $cluster (sort keys %Genus_cluster2gn){
	my @Gns = split (/\t/, $Genus_cluster2gn{$cluster});
	my $gn_num = scalar @Gns;
	if ($gn_num == 1){
		$Genus_cluster_singleton{$cluster} = 1;
	}
}

# Step 4 Compare genus cluter to family cluster
my @Prec_same_family_cluster = (); # The percentage of genomes from each genus fall into the same family
foreach my $cluster (sort keys %Genus_cluster2gn){
	my @Gns = split (/\t/, $Genus_cluster2gn{$cluster});
	
	my %Family_cluster_freq = (); # $family_cluster => number of frequency of family cluster appearence
	foreach my $gn (@Gns){
		my $family_cluster = $Gn2family_cluster{$gn};
		$Family_cluster_freq{$family_cluster}++;
	}
	
	my @Family_cluster_freq_array = sort { $Family_cluster_freq{$a} <=> $Family_cluster_freq{$b} } keys %Family_cluster_freq;
	
	my $prec_same_family_cluster = $Family_cluster_freq{$Family_cluster_freq_array[-1]} / (scalar @Gns);
	push @Prec_same_family_cluster, $prec_same_family_cluster;
}

my @Prec_same_family_cluster_exclude_singleton = (); # The percentage of genomes from each genus fall into the same family, excluding singleton genus clusters
foreach my $cluster (sort keys %Genus_cluster2gn){
	if (!exists $Genus_cluster_singleton{$cluster}){
		my @Gns = split (/\t/, $Genus_cluster2gn{$cluster});
		
		my %Family_cluster_freq = (); # $family_cluster => number of frequency of family cluster appearence
		foreach my $gn (@Gns){
			my $family_cluster = $Gn2family_cluster{$gn};
			$Family_cluster_freq{$family_cluster}++;
		}
		
		my @Family_cluster_freq_array = sort { $Family_cluster_freq{$a} <=> $Family_cluster_freq{$b} } keys %Family_cluster_freq;
		
		my $prec_same_family_cluster_exclude_singleton = $Family_cluster_freq{$Family_cluster_freq_array[-1]} / (scalar @Gns);
		push @Prec_same_family_cluster_exclude_singleton, $prec_same_family_cluster_exclude_singleton;
	}
}

# Step 5 Write result
my $gn_num_by_genus_cluster = scalar (keys %Gn2genus_cluster);
my $gn_num_by_family_cluster = scalar (keys %Gn2family_cluster);
my $genus_cluster_num = scalar (keys %Genus_cluster2gn);
my $family_cluster_num = scalar (keys %Family_cluster2gn);
my $singleton_genus_cluster_num = scalar %Genus_cluster_singleton;
my $prec_same_family_cluster_average = _avg(@Prec_same_family_cluster);
my $prec_same_family_cluster_exclude_singleton_average = _avg(@Prec_same_family_cluster_exclude_singleton);

open OUT, ">Cluster_phage_genomes/Compare_genus_cluster_and_family_cluster_result.txt";
print OUT "Total number of viral genomes (in genus clusters): $gn_num_by_genus_cluster\n";
print OUT "Total number of viral genomes (in family clusters): $gn_num_by_family_cluster\n";
print OUT "Total number of genus cluster: $genus_cluster_num\n";
print OUT "Total number of family cluster: $family_cluster_num\n";
print OUT "Total number of singleton genus cluster: $singleton_genus_cluster_num\n";
print OUT "The average percentage of genomes from each genus cluster fall into the same family cluster: $prec_same_family_cluster_average\n";
print OUT "The average percentage of genomes from each genus cluster fall into the same family cluster (excluding singleton genus clusters): $prec_same_family_cluster_exclude_singleton_average\n";
close OUT;

# Subroutine
sub _avg { # Input is an array
    my $total;
    $total += $_ foreach @_;
    # sum divided by number of components.
    return $total / @_;
}