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
=pod
# Step 3 Make phage genome lists for each genus cluster (excluding singletons)
my $dRep_dir = "/storage1/data11/TYMEFLIES_phage/dRep_working_dir";
`mkdir $dRep_dir`;
`mkdir $dRep_dir/phage_genome_list`;
`mkdir $dRep_dir/phage_genomes`;

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

## Copy all virus genomes into phage_genome folder
my $source_dir = "/storage1/data11/TYMEFLIES_phage/33*/vRhyme_best_bins_fasta_parsed/";
my $target_dir = "/storage1/data11/TYMEFLIES_phage/dRep_working_dir/phage_genomes/";

my $find_output = `find $source_dir -type f -name '*.fasta'`;

my @files = split(/\n/, $find_output);

foreach my $file (@files) {
    system("cp $file $target_dir") == 0 or die "Error copying file $file: $!";
}

# Step 4 Run dRep for all genus clusters
open OUT, ">tmp.run_dRep.sh";
open IN, "find $dRep_dir/phage_genome_list/ -name 'phage_genome_list.*.txt' |";
while (<IN>){
	chomp;
	my $line = $_;
	my ($genus_cluster) = $line =~ /^.+\/phage_genome_list\.(.+?)\.txt/;
	print OUT "dRep dereplicate dRep_working_dir/Output.$genus_cluster -p 1 -g $dRep_dir/phage_genome_list/phage_genome_list.$genus_cluster.txt -l 2000 --ignoreGenomeQuality -pa 0.8 -sa 0.95 -nc 0.85 -comW 0 -conW 0 -strW 0 -N50W 0 -sizeW 1 -centW 0\n";
}
close IN;

`cat tmp.run_dRep.sh | parallel -j 20`;

`rm tmp.run_dRep.sh`;
=cut
# Step 5 Make all viral genome to viral sequences hash and find viral genomes that are not included in genus cluster
## Concatenate all virus genomes
`find /storage1/data11/TYMEFLIES_phage/33*/vRhyme_best_bins_fasta_parsed/ -name '*.fasta' -exec cat {} + > TYMEFLIES_all_phages.fasta`;
my %Viral_gn2viral_seq = (); # $viral_gn => $viral_seq collection separated by "\t"
open IN, "/storage1/data11/TYMEFLIES_phage/TYMEFLIES_all_phages.fasta";
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

## Remove the "TYMEFLIES_all_phages.fasta"
`rm TYMEFLIES_all_phages.fasta`; 

my %Viral_gn_not_in_genus_cluster = (); # $viral_gn => $viral_seq collection separated by "\t"
foreach my $viral_gn (sort keys %Viral_gn2viral_seq){
	if (! exists $Gn2genus_cluster{$viral_gn}){
		$Viral_gn_not_in_genus_cluster{$viral_gn} = $Viral_gn2viral_seq{$viral_gn};
	}
}

# Step 6 Store species-level vOTUs
# Species-level vOTUs contains:
# 1) each represenative from genus cluster
# 2) singleton genus cluster genome
# 3) viral genome that is not included in genus cluster

my %Species_level_vOTUs = (); # $gn => [0] $gns_all (all genomes within each vOTUs) [1] $genus_cluster (the name the genus cluster or just not included in any genus cluster) 
## Step 6.1 Grep representiative from each genus cluster
open IN, "find dRep_working_dir/ -name 'Output.Genus_cluster*' | ";
while (<IN>){
	chomp;
	my $folder = $_;
	my ($genus_cluster) = $folder =~ /^.+\/Output\.(.+?)$/;
	
	# Store cluster to representive genome and all genome hash
	my %Hash = (); # $cluster => [0] $rep_gn [1] $gn collection separated by "\,"
	open INN, "$folder/data_tables/Cdb.csv" or warn "$folder/data_tables/Cdb.csv is not present";
	while (<INN>){
		chomp;
		if (!/^genome/){
			my @tmp = split (/\,/);
			my ($gn) = $tmp[0] =~ /^(.+?)\.fasta/;
			my $cluster = $tmp[1];
			if (!exists $Hash{$cluster}[1]){
				$Hash{$cluster}[1] = $gn;
			}else{
				$Hash{$cluster}[1] .= "\,".$gn;
			}
		}
	}
	close INN;
	
	open INN, "$folder/data_tables/Wdb.csv" or warn "$folder/data_tables/Wdb.csv is not present";
	while (<INN>){
		chomp;
		if (!/^genome/){
			my @tmp = split (/\,/);
			my ($gn) = $tmp[0] =~ /^(.+?)\.fasta/;
			my $cluster = $tmp[1];
			$Hash{$cluster}[0] = $gn;

		}
	}
	close INN;	
	
	# Add content of %Hash into the %Species_level_vOTUs
	foreach my $cluster (sort keys %Hash){
		my $rep_gn = $Hash{$cluster}[0];
		my $gns_all = $Hash{$cluster}[1];
		$Species_level_vOTUs{$rep_gn}[0] = $gns_all;
		$Species_level_vOTUs{$rep_gn}[1] = $genus_cluster;
	}
}
close IN;

## Step 6.2 Store singleton genus cluster
foreach my $genus_cluster (sort keys %Genus_cluster_singleton){
	my $gn = $Genus_cluster2gn{$genus_cluster};
	$Species_level_vOTUs{$gn}[0] = $gn;
	$Species_level_vOTUs{$gn}[1] = $genus_cluster;
}

## Step 6.3 Store viral genomes that are not included in genus cluster
foreach my $viral_gn (sort keys %Viral_gn_not_in_genus_cluster){
	$Species_level_vOTUs{$viral_gn}[0] = $viral_gn;
	$Species_level_vOTUs{$viral_gn}[1] = "not included in any genus cluster";
}

# Step 7 Write down species-level vOTU result
open OUT, ">Cluster_phage_genomes/Species_level_vOTUs_cluster.txt";
foreach my $gn (sort keys %Species_level_vOTUs){
	print OUT "$gn\t$Species_level_vOTUs{$gn}[0]\t$Species_level_vOTUs{$gn}[1]\n";
}
close OUT;

# Step 8 Write down species-level vOTU each cluster genome numbers
open OUT, ">Cluster_phage_genomes/Species_level_vOTUs_cluster2genome_number.txt";
foreach my $gn (sort keys %Species_level_vOTUs){
	my $gns = $Species_level_vOTUs{$gn}[0];
	my @Gns = split (/\,/, $gns);
	my $gn_num = scalar @Gns;
	print OUT "$gn\t$gn_num\t$Species_level_vOTUs{$gn}[1]\n";
}
close OUT;

`cat Cluster_phage_genomes/Species_level_vOTUs_cluster2genome_number.txt | sort -k 2 -n > tmp`;
`rm Cluster_phage_genomes/Species_level_vOTUs_cluster2genome_number.txt`;
`mv tmp Cluster_phage_genomes/Species_level_vOTUs_cluster2genome_number.txt`;