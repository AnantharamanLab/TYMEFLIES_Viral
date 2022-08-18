#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw(sum);

# Aim: Calculate gene coverage for AMG-containing viral gn (species representatives) for each year

# Step 1 Get AMG-containing viral gn list (species representatives)
## Step 1.1 Store species info
my %Species = (); # $gn_rep => $gns 
open IN, "/storage1/data11/TYMEFLIES_phage/Cluster_phage_genomes/Species_level_vOTUs_cluster.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $gn_rep = $tmp[0];
	my $gns = $tmp[1];
	$Species{$gn_rep} = $gns;
}
close IN;

# Step 2 Calculate gene depth
## Step 2.1 Store the hash %Viral_gn2viral_scf and hash %Viral_scf2length
my %All_phage_species_rep_gn_containing_AMG_seq = ();
%All_phage_species_rep_gn_containing_AMG_seq = _store_seq("reference_fasta_for_metapop/All_phage_species_rep_gn_containing_AMG.fasta");

my %Viral_gn2viral_scf = (); # $viral_gn => collection of $viral_scf, separated by "\t"
my %Viral_scf2length = (); # $viral_scf => $length
foreach my $key (sort keys %All_phage_species_rep_gn_containing_AMG_seq){
	my ($viral_scf) = $key =~ /^>(.+?)$/;
	my ($viral_gn) = $viral_scf =~ /^(.+?\_\_vRhyme.+?)\_\_/;
	if (!exists $Viral_gn2viral_scf{$viral_gn}){
		$Viral_gn2viral_scf{$viral_gn} = $viral_scf;
	}else{
		$Viral_gn2viral_scf{$viral_gn} .= "\t".$viral_scf;
	}
	
	my $seq = $All_phage_species_rep_gn_containing_AMG_seq{$key}; # The sequence of the scaffold
	
	my $length = length($seq);
	$Viral_scf2length{$viral_scf} = $length;
}

## Step 2.2 Store and write down all gene coordinates 
my %All_gene_coordinates = (); # $gene => $gene, $viral_scf, $viral_scf_length, $start, $stop, $gene_size connected by "\t"
my %Viral_scf2genes = (); # $viral_scf => collection of $gene\|$start\|$stop, separated by "\,"
open IN, "All_phage_species_rep_gn_containing_AMG.mdfed.genes"; # Note: Using mdfed gene files here 
while (<IN>){
	chomp;
	if (/^>/){
		my $line = $_;
		my ($gene, $start, $stop) = $line =~ /^>(.+?) \# (.+?) \# (.+?) \#/;
		my ($viral_scf) = $gene =~ /^(.+)\_/;
		my $viral_scf_length = $Viral_scf2length{$viral_scf};
		my $gene_size = $stop - $start;
		$All_gene_coordinates{$gene} = "$gene\t$viral_scf\t$viral_scf_length\t$start\t$stop\t$gene_size";
		
		my $gene_start_stop = "$gene\|$start\|$stop";
		if (!exists $Viral_scf2genes{$viral_scf}){
			$Viral_scf2genes{$viral_scf} = $gene_start_stop;
		}else{
			$Viral_scf2genes{$viral_scf} .= "\,".$gene_start_stop;
		}
	}
}
close IN;

open OUT, ">MetaPop.for_each_year/MetaPop/All_phage_species_rep_gn_gene_coordinates.txt";
print OUT "gene\tscaffold\tscaffold length\tstart\tstop\tgene size\n";
foreach my $gene (sort keys %All_gene_coordinates){
	print OUT "$All_gene_coordinates{$gene}\n";
}
close OUT;

## Step 2.3 Store the depth per pos info
### Firstly assign all pos depth as 0 to set the range of positions
my %Viral_scf2pos2depth = (); # $viral_scf => $pos => 0; firstly assign all pos depth as 0
foreach my $viral_scf (sort keys %Viral_scf2length){
	my $length = $Viral_scf2length{$viral_scf};
	# Only take sub-region of a given scaffold, cut the start and stop 150 bp regions
	my $start_pos = 151;
	my $end_pos = $length - 150;
	my @Pos = ($start_pos..$end_pos);
	foreach my $pos (@Pos){
		$Viral_scf2pos2depth{$viral_scf}{$pos} = 0;
	}
}

### Process each *depth_per_pos.tsv to get viral gene 2 year 2 depth file
`mkdir MetaPop.for_each_year/MetaPop/04.Depth_per_Pos/Gene_coverage_for_each_year`;

open IN, "ls MetaPop.for_each_year/MetaPop/04.Depth_per_Pos/*.tsv |";
while (<IN>){
	chomp;
	my $file = $_;
	my ($year) = $file =~ /04\.Depth\_per\_Pos\/(\d+?)\./;
	
	# Store the depth for each position
	my %Viral_scf2pos2depth_for_year = %Viral_scf2pos2depth;
	open INN, $file;
	while (<INN>){
		chomp;
		my @tmp = split (/\t/);
		my $viral_scf = $tmp[0];
		my $pos = $tmp[1];
		my $depth = $tmp[2];
		
		if (exists $Viral_scf2pos2depth_for_year{$viral_scf}{$pos}){
			$Viral_scf2pos2depth_for_year{$viral_scf}{$pos} = $depth;
		}
		
	}
	close INN;
	
	my %Viral_gene2depth = (); # $gene => $depth
	# Calculate each gene depth
	foreach my $viral_scf (sort keys %Viral_scf2genes){
		my @Gene_start_stop = split (/\,/, $Viral_scf2genes{$viral_scf});
		foreach my $gene_start_stop (@Gene_start_stop){
			my ($gene, $start, $stop) = $gene_start_stop =~ /^(.+?)\|(.+?)\|(.+?)$/;
			my $gene_length = $stop - $start;
			if ($gene_length < 450){ # If the gene length is less than 450 bp, we do not use it
				$Viral_gene2depth{$gene} = "NA";
			}else{
				my $gene_depth = 0;
				my @Gene_range = ($start..$stop);
				my @Gene_depth = (); # Store all the depth value within the gene region
				foreach my $pos (@Gene_range){
					if (exists $Viral_scf2pos2depth_for_year{$viral_scf}{$pos}){
						push @Gene_depth, $Viral_scf2pos2depth_for_year{$viral_scf}{$pos};
					}
				}
				
				$gene_depth = _mean(@Gene_depth);
				$Viral_gene2depth{$gene} = $gene_depth;
			}
		}
	}
	
	open OUT, ">MetaPop.for_each_year/MetaPop/04.Depth_per_Pos/Gene_coverage_for_each_year/$year.gene_coverage.txt";
	foreach my $gene (sort keys %Viral_gene2depth){
		print OUT "$gene\t$Viral_gene2depth{$gene}\n";
	}
	close OUT;
}
close IN;



## Subroutine
sub _store_seq{ # Store seq; the head will be truncated from the first space if space(s) are present
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

sub _mean { # Get mean value within an array
    return sum(@_)/@_;
}