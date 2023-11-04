#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw(sum);

# Aim: Get gene gain and loss for AMG-containing viral gn (species representatives) for each year

# Step 1. Store Viral_gene2year2gene_freq.txt result
my %Viral_gene2year2gene_freq = (); # $viral_gene => $year => $gene_freq
my @Head = (); # Store the header
open IN, "MetaPop.for_each_year/MetaPop/Viral_gene2year2gene_freq.txt";
while (<IN>){
	chomp;
	if (/^Viral/){
		my @tmp = split (/\t/);
		@Head = @tmp;
	}else{
		my @tmp = split (/\t/);
		my $viral_gene = $tmp[0];
		for(my $i=1; $i<=$#tmp; $i++){
			my $gene_freq = $tmp[$i];
			my $year = $Head[$i];
			$Viral_gene2year2gene_freq{$viral_gene}{$year} = $gene_freq;
		}
	}
}
close IN;

# Step 2 Get the gene gain and loss result based on gene frequency change between 2000-2003 and year 2016-2019 periods
my %Viral_gene2gain_or_loss = (); # $viral_gn => "gain" or "loss"
foreach my $viral_gene (sort keys %Viral_gene2year2gene_freq){
	my $year_2000_2003_gene_freq_mean = "NA";
	my $year_2016_2019_gene_freq_mean = "NA";
	
	my @Year_2000_2003 = (2000..2003);
	my @Year_2016_2019 = (2016..2019);
	
	my @Year_2000_2003_gene_freq = ();
	my @Year_2016_2019_gene_freq = ();
	
	for(my $i=0; $i<=$#Year_2000_2003; $i++){
		my $year = $Year_2000_2003[$i];
		my $gene_freq = $Viral_gene2year2gene_freq{$viral_gene}{$year};
		if ($gene_freq ne "NA"){
			push @Year_2000_2003_gene_freq, $gene_freq;
		}
	}

	for(my $i=0; $i<=$#Year_2016_2019; $i++){
		my $year = $Year_2016_2019[$i];
		my $gene_freq = $Viral_gene2year2gene_freq{$viral_gene}{$year};
		if ($gene_freq ne "NA"){
			push @Year_2016_2019_gene_freq, $gene_freq;
		}
	}	
	
	if (@Year_2000_2003_gene_freq){
		$year_2000_2003_gene_freq_mean = _mean(@Year_2000_2003_gene_freq);
	}
	
	if (@Year_2016_2019_gene_freq){
		$year_2016_2019_gene_freq_mean = _mean(@Year_2016_2019_gene_freq);
	}	
	
	if ($year_2000_2003_gene_freq_mean ne "NA" and $year_2016_2019_gene_freq_mean ne "NA" and abs($year_2000_2003_gene_freq_mean - $year_2016_2019_gene_freq_mean) >= 1){
		if ($year_2000_2003_gene_freq_mean > $year_2016_2019_gene_freq_mean){
			$Viral_gene2gain_or_loss{$viral_gene} = "loss";
		}else{
			$Viral_gene2gain_or_loss{$viral_gene} = "gain";
		}
	}
}

# Step 3 Write down gene gain and loss result
open OUT, ">MetaPop.for_each_year/MetaPop/Viral_gene2gain_or_loss.txt";
foreach my $viral_gene (sort keys %Viral_gene2gain_or_loss){
	print OUT "$viral_gene\t$Viral_gene2gain_or_loss{$viral_gene}\n";
}
close OUT;

# Step 4 Generate gene gain and loss result table for each viral genome and write it down
## Step 4.1 Get all the viral genomes that have either gain or loss genes 
my %Viral_gn_having_gain_or_loss_genes = (); # $viral_gn => the collection of $viral_genes (either gain or loss genes)
foreach my $viral_gene (sort keys %Viral_gene2gain_or_loss){
	my ($viral_gn_from_gene) = $viral_gene =~ /^(.+?)\_\_Ga/; 
	if (!exists $Viral_gn_having_gain_or_loss_genes{$viral_gn_from_gene}){
		$Viral_gn_having_gain_or_loss_genes{$viral_gn_from_gene} = $viral_gene;
	}else{
		$Viral_gn_having_gain_or_loss_genes{$viral_gn_from_gene} .= "\t".$viral_gene;
	}
}

## Step 4.2 Generate the dict of Viral_gn2gain_or_loss_info
my %Viral_gn2gain_or_loss_info = (); # $viral_gn => $num_of_gene_gain \t $num_of_gene_loss \t $genes_gain \t $genes_loss
for my $viral_gn (sort keys %Viral_gn_having_gain_or_loss_genes){	
	my $num_of_gene_gain = 0;
	my $num_of_gene_loss = 0;
	my $genes_gain = "";
	my $genes_loss = "";
	my @Genes_gain = ();
	my @Genes_loss = ();	
	
	my @Viral_genes_either_gain_or_loss = split("\t", $Viral_gn_having_gain_or_loss_genes{$viral_gn});
	
	foreach my $viral_gene (sort @Viral_genes_either_gain_or_loss){
		if ($Viral_gene2gain_or_loss{$viral_gene} eq "gain"){
			$num_of_gene_gain++;
			push @Genes_gain, $viral_gene;
		}elsif($Viral_gene2gain_or_loss{$viral_gene} eq "loss"){
			$num_of_gene_loss++;
			push @Genes_loss, $viral_gene;
		}	
	}
	
	if (@Genes_gain){
		$genes_gain = join("\,", @Genes_gain);
	}

	if (@Genes_loss){
		$genes_loss = join("\,", @Genes_loss);
	}

	$Viral_gn2gain_or_loss_info{$viral_gn} = $num_of_gene_gain."\t".$num_of_gene_loss."\t".$genes_gain."\t".$genes_loss;
}
	
## Step 4.3 Write down the Viral_gn2gain_or_loss_info.txt table
open OUT, ">MetaPop.for_each_year/MetaPop/Viral_gn2gain_or_loss_info.txt";
foreach my $viral_gn (sort keys %Viral_gn2gain_or_loss_info){
	print OUT "$viral_gn\t$Viral_gn2gain_or_loss_info{$viral_gn}\n";
}
close OUT;

# Step 5 Write down gene gain and loss result for four AMGs
## Step 5.1 Store the viral species containing four AMGs
my %Viral_species_containing_four_AMGs = (); # $viral_gn => $info (the corresponding information for each viral species gn)
my @Viral_species_containing_four_AMGs = (); # Store the $viral_gn order
open IN, "viral_species_containing_four_AMGs.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $info = $tmp[0];
	my $viral_gn = $tmp[1];
	$Viral_species_containing_four_AMGs{$viral_gn} = $info;
	
	push @Viral_species_containing_four_AMGs, $viral_gn;
}
close IN;

## Step 5.2 Write down gene gain and loss result for four AMGs
open OUT, ">MetaPop.for_each_year/MetaPop/Viral_gene2gain_or_loss.for_four_AMGs.txt";
print OUT "Viral gn\tViral gn characteristics\tNum of gene gain\tNum of gene loss\tGenes gained\tGenes lost\n";
foreach my $viral_gn (@Viral_species_containing_four_AMGs){
	my $num_of_gene_gain = 0;
	my $num_of_gene_loss = 0;
	my @Genes_gain = ();
	my @Genes_loss = ();
	my $genes_gain = "";
	my $genes_loss = "";
	
	foreach my $viral_gene (sort keys %Viral_gene2gain_or_loss){
		my ($viral_gn_from_gene) = $viral_gene =~ /^(.+?)\_\_Ga/; 
		if ($viral_gn_from_gene eq $viral_gn){
			if ($Viral_gene2gain_or_loss{$viral_gene} eq "gain"){
				$num_of_gene_gain++;
				push @Genes_gain, $viral_gene;
			}elsif($Viral_gene2gain_or_loss{$viral_gene} eq "loss"){
				$num_of_gene_loss++;
				push @Genes_loss, $viral_gene;
			}
		}
	}	
	
	if (@Genes_gain){
		$genes_gain = join("\,", @Genes_gain);
	}

	if (@Genes_loss){
		$genes_loss = join("\,", @Genes_loss);
	}	
	
	print OUT "$viral_gn\t$Viral_species_containing_four_AMGs{$viral_gn}\t$num_of_gene_gain\t$genes_gain\t$num_of_gene_loss\t$genes_loss\n";
}
close OUT;



## Subroutine

sub _mean { # Get mean value within an array
    return sum(@_)/@_;
}