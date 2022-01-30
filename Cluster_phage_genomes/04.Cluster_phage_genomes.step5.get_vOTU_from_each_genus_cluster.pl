#!/usr/bin/perl

use strict;
use warnings;

# Aim: Get vOTU (species level) from each genus cluster

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

# Step 2 Find singleton genus cluster
my %Genus_cluster_singleton = (); # $cluster (singleton genus cluster) => 1
foreach my $cluster (sort keys %Genus_cluster2gn){
	my @Gns = split (/\t/, $Genus_cluster2gn{$cluster});
	my $gn_num = scalar @Gns;
	if ($gn_num == 1){
		$Genus_cluster_singleton{$cluster} = 1;
	}
}

# Step 3 Make phage genome lists for each genus cluster (excluding singletons)
my $dRep_dir = "/storage1/data11/TYMEFLIES_phage/dRep_working_dir";
`mkdir $dRep_dir/phage_genome_list`;

foreach my $genus_cluster (sort keys %Genus_cluster2gn){
	if (! exists $Genus_cluster_singleton{$genus_cluster}){ # Not singletons
		my @Gns = split(/\t/, $Genus_cluster2gn{$genus_cluster});
		open OUT, ">$dRep_dir/phage_genome_list/phage_genome_list.$genus_cluster.txt";
		foreach my $gn (@Gns){
			print OUT "$dRep_dir/phage_genomes/$gn.fasta\n";
		}
		close OUT;
	}
}
=pod
# Step 4 Make all viral genome to viral sequences hash and find viral genomes that are not included in genus cluster
my %Viral_gn2viral_seq = (); # $viral_gn => $viral_seq collection separated by "\t"
open IN, "/storage1/data11/TYMEFLIES_phage/Host_prediction/All_phage_genomes.fasta";
while (<IN>){
	chomp;
	if (/^>/){
		my $line = $_;
		my ($viral_seq) = $line =~ /^>(.+?)$/;
		my ($viral_gn) = $viral_seq =~ /^(.+?\_\_.+?)\_\_.+?$/;
		if (! exists $Viral_gn2viral_seq{$viral_gn}){
			$Viral_gn2viral_seq{$viral_gn} = $viral_seq;
		}else{
			$Viral_gn2viral_seq{$viral_gn} .= "\t".$viral_seq;
		}
	}
}

my %Viral_gn_not_in_genus_cluster = (); # $viral_gn => $viral_seq collection separated by "\t"
foreach my $viral_gn (sort keys %Viral_gn2viral_seq){
	if (! exists $Gn2genus_cluster{$viral_gn}){
		$Viral_gn_not_in_genus_cluster{$viral_gn} = $Viral_gn2viral_seq{$viral_gn};
	}
}
=cut
# Step 5 Run dRep for all genus clusters
open OUT, ">tmp.run_dRep.sh";
open IN, "find $dRep_dir/phage_genome_list/ -name 'phage_genome_list.*.txt' |";
while (<IN>){
	chomp;
	my $line = $_;
	my ($genus_cluster) = $line =~ /^.+\/phage_genome_list\.(.+?)\.txt/;
	print OUT "dRep dereplicate dRep_working_dir/Output.$genus_cluster -p 8 -g $dRep_dir/phage_genome_list/phage_genome_list.$genus_cluster.txt -l 2000 --ignoreGenomeQuality -pa 0.8 -sa 0.95 -nc 0.85 -comW 0 -conW 0 -strW 0 -N50W 0 -sizeW 1 -centW 0\n";
}
close IN;

`cat tmp.run_dRep.sh | parallel -j 2`;

`rm tmp.run_dRep.sh`;


