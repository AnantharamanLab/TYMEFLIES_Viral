#!/usr/bin/perl
  
use strict;
use warnings;

# Aim: Get all the viral species gene sequences

# Step 1 Store the viral gn sequences
## Store allthe viral species sequences (rep gn)
my %All_viral_gn_seq = _store_seq("reference_fasta_for_metapop/All_phage_species_rep_gn.fasta");

# Step 2 Get all the viral gn hash
my %Viral_rep_gn = (); # $viral_rep_gn => 1
foreach my $key (sort keys %All_viral_gn_seq){
	my ($key_wo_arrow) = $key =~ /^>(.+?)$/;
	my ($viral_rep_gn) = $key_wo_arrow =~ /^(33.+?\_\_.+?)\_\_/;
	$Viral_rep_gn{$viral_rep_gn} = 1;	
}

# Step 3 Store the viral gn gene sequences for each viral genome
my %All_viral_gn_gene_seq = _store_seq("All_phage_species_rep_gn.mdfed.genes");

# Step 5 Write down viral species gene sequences containing four AMGs
`mkdir Summer_vs_Winter_Fst_analysis/Genes`;

foreach my $key (sort keys %All_viral_gn_gene_seq){
	my ($viral_rep_gn) = $key =~ /^>(33.+?\_\_.+?)\_\_/;
	if (exists $Viral_rep_gn{$viral_rep_gn}){
		open OUT, ">>Summer_vs_Winter_Fst_analysis/Genes/$viral_rep_gn.genes";
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


