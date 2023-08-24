#!/usr/bin/perl
  
use strict;
use warnings;

# Aim: Get the viral species sequences for four AMGs

# Step 1 Store the viral species containing four AMGs
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

# Step 2 Store the viral gn sequences
my %All_viral_gn_seq = _store_seq("All_phage_species_rep_gn_containing_AMG.fasta");

# Step 3 Write down viral species sequences containing four AMGs
open OUT, ">Summer_vs_Winter_Fst_analysis/Viral_species_containing_four_AMGs.fasta";
foreach my $key (sort keys %All_viral_gn_seq){
	my ($viral_gn) = $key =~ /^>(.+?\_\_.+?)\_\_/; # Extract the viral gn from scaffold ID
	if (exists $Viral_species_containing_four_AMGs{$viral_gn}){
		print OUT "$key\n$All_viral_gn_seq{$key}\n";
	}
}
close OUT;

# Step 4 Store the viral gn gene sequences for each viral genome
my %All_viral_gn_gene_seq = _store_seq("All_phage_species_rep_gn_containing_AMG.mdfed.genes");

# Step 5 Write down viral species gene sequences containing four AMGs
`mkdir Summer_vs_Winter_Fst_analysis/Genes`;

foreach my $key (sort keys %All_viral_gn_gene_seq){
	my ($viral_gn) = $key =~ /^>(.+?\_\_.+?)\_\_/; # Extract the viral gn from scaffold ID
	if (exists $Viral_species_containing_four_AMGs{$viral_gn}){
		open OUT, ">>Summer_vs_Winter_Fst_analysis/Genes/$viral_gn.genes";
		print OUT "$key\n$All_viral_gn_gene_seq{$key}\n";
		close OUT;
	}
}



## Subroutine

sub _store_seq{ # Store the entire header line
	my $file = $_[0];
	my %Seq = (); my $head = "";
	open _IN, "$file";
	while (<_IN>){
		chomp;
		if (/>/){
			($head) = $_ =~ /^(>.+?)$/;
			$Seq{$head} = "";
		}else{
			$Seq{$head} .= $_;
		}
	}
	close _IN;
	return %Seq;
}


