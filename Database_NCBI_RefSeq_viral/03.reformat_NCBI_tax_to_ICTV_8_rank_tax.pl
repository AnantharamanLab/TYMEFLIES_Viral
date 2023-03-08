#!/usr/bin/perl

use strict;
use warnings;

# AIM: reformat NCBI tax of viruses to ICTV 8 rank tax

# Step 1. Store NCBI tax and dereplicated it
my %ID2NCBI_tax = (); # $id (protein ID; for example, YP_009854875) => $ncbi_tax
my %NCBI_tax = (); # $ncbi_tax => 1
open IN, "viral.protein.tax.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $id = $tmp[0];
	my $ncbi_tax = $tmp[1];
	my $ncbi_species = $tmp[2];
	$ID2NCBI_tax{$id} = $ncbi_tax;
	$NCBI_tax{$ncbi_tax} = $ncbi_species;
}
close IN;

# Step 2. Store ICTV tax
my %ICTV_tax = (); # $sort_id => [0] $tax_full (separated by ";") [1] $tax_8_rank [2] $genome_composition (dsDNA or ssDNA or other types)
my %Realm = (); # $realm = 1;
my %Subrealm = (); # $subrealm = 1;
my %Kingdom = (); # $kingdom = 1;
my %Subkingdom = (); # $subkingdom = 1;
my %Phylum = (); # $phylum = 1;
my %Subphylum = (); # $subphylum = 1;
my %Class = (); # $class = 1;
my %Subclass = (); # $subclass = 1;
my %Order = (); # $order = 1;
my %Suborder = (); # $suborder = 1;
my %Family = (); # $family = 1;
my %Subfamily = (); # $subfamily = 1;
my %Genus = (); # $genus = 1;
my %Subgenus = (); # $subgenus = 1;
my %Species = (); # $species = 1;
my %String2rank = (); # $string => $rank (realm, subrealm ... )

open IN, "ICTV_Master_Species_List_2021.v3.txt";
while (<IN>){
	chomp;
	if (!/^Sort/){
		my @tmp = split (/\t/);
		my $sort_id = $tmp[0];
		my $realm = $tmp[1]; my $subrealm = $tmp[2];		
		my $kingdom = $tmp[3]; my $subkingdom = $tmp[4];
		my $phylum = $tmp[5]; my $subphylum = $tmp[6];
		my $class = $tmp[7]; my $subclass = $tmp[8];
		my $order = $tmp[9]; my $suborder = $tmp[10];
		my $family = $tmp[11]; my $subfamily = $tmp[12];
		my $genus = $tmp[13]; my $subgenus = $tmp[14];
		my $species = $tmp[15];
		my $genome_composition = $tmp[16];
		
		$Realm{$realm} = 1; $Subrealm{$subrealm} = 1;
		$Kingdom{$kingdom} = 1; $Subkingdom{$subkingdom} = 1;
		$Phylum{$phylum} = 1; $Subphylum{$subphylum} = 1;
		$Class{$class} = 1; $Subclass{$subclass} = 1;
		$Order{$order} = 1; $Suborder{$suborder} = 1;
		$Family{$family} = 1; $Subfamily{$subfamily} = 1;
		$Genus{$genus} = 1; $Subgenus{$subgenus} = 1;
		$Species{$species} = 1;
		
		$String2rank{$realm} = "Realm";
		$String2rank{$subrealm} = "Subrealm";
		$String2rank{$kingdom} = "Kingdom";
		$String2rank{$subkingdom} = "Subkingdom";
		$String2rank{$phylum} = "Phylum";
		$String2rank{$subphylum} = "Subphylum";		
		$String2rank{$class} = "Class";
		$String2rank{$subclass} = "Subclass";	
		$String2rank{$order} = "Order";
		$String2rank{$suborder} = "Suborder";	
		$String2rank{$family} = "Family";
		$String2rank{$subfamily} = "Subfamily";
		$String2rank{$genus} = "Genus";
		$String2rank{$subgenus} = "Subgenus";
		$String2rank{$species} = "Species";	
		
		my $tax_full = "$realm;$subrealm;$kingdom;$subkingdom;$phylum;$subphylum;$class;$subclass;$order;$suborder;$family;$subfamily;$genus;$subgenus;$species";
		my $tax_8_rank = "$realm;$kingdom;$phylum;$class;$order;$family;$genus;$species";
		
		$ICTV_tax{$sort_id}[0] = $tax_full;
		$ICTV_tax{$sort_id}[1] = $tax_8_rank;
		$ICTV_tax{$sort_id}[2] = $genome_composition;
	}
}
close IN;

my @Full_rank = qw/Realm Subrealm Kingdom Subkingdom Phylum Subphylum Class Subclass Order Suborder Family Subfamily Genus Subgenus Species/;

my %NCBI_tax2ictv_8_rank_tax = (); # $ncbi_tax => $ictv_8_rank_tax
my %NCBI_tax_novel_rank = (); # Store the $novel_rank => $ncbi_tax

# Step 3. Reformat the NCBI tax to 8-rank ICTV tax
foreach my $ncbi_tax (sort keys %NCBI_tax){
	my $species = $NCBI_tax{$ncbi_tax}; # The "species" of the full tax
	$ncbi_tax =~ s/^Viruses\;//g; # Delete the "Viruses;" in the front
	my $ictv_8_rank_tax = "";
	my $realm = "";
	my $kingdom = "";
	my $phylum = "";
	my $class = "";
	my $order = "";
	my $family = "";
	my $genus = "";
	
	my @NCBI_tax = split (/\;/,$ncbi_tax);
	for(my $i=0; $i<=$#NCBI_tax; $i++){
		if ($String2rank{$NCBI_tax[$i]}){
			if ($String2rank{$NCBI_tax[$i]} eq "Realm"){
				$realm = $NCBI_tax[$i];
			}elsif ($String2rank{$NCBI_tax[$i]} eq "Kingdom"){
				$kingdom = $NCBI_tax[$i];
			}elsif ($String2rank{$NCBI_tax[$i]} eq "Phylum"){
				$phylum = $NCBI_tax[$i];
			}elsif ($String2rank{$NCBI_tax[$i]} eq "Class"){
				$class = $NCBI_tax[$i];
			}elsif ($String2rank{$NCBI_tax[$i]} eq "Order"){
				$order = $NCBI_tax[$i];
			}elsif ($String2rank{$NCBI_tax[$i]} eq "Family"){
				$family = $NCBI_tax[$i];
			}elsif ($String2rank{$NCBI_tax[$i]} eq "Genus"){
				$genus = $NCBI_tax[$i];
			}elsif ($String2rank{$NCBI_tax[$i]} eq "Species"){
				$species = $NCBI_tax[$i];
			}
		}else{
			my $novel_rank = $NCBI_tax[$i];
			$NCBI_tax_novel_rank{$novel_rank} = $ncbi_tax;
		}
	}
	
	$ictv_8_rank_tax = "$realm;$kingdom;$phylum;$class;$order;$family;$genus;$species";
	if ($ictv_8_rank_tax ne ";;;;;;;"){
		$NCBI_tax2ictv_8_rank_tax{$ncbi_tax} = $ictv_8_rank_tax;  # Note that here $ncbi_tax does not contain "Viruses;" in the front
	}	
}

# Step 4. Write down result
open OUT, ">viral.protein.ictv_8_rank_tax.txt";
foreach my $id (sort keys %ID2NCBI_tax){
	my $ncbi_tax = $ID2NCBI_tax{$id};
	$ncbi_tax =~ s/^Viruses\;//g; # Delete the "Viruses;" in the front
	my $ictv_8_rank_tax = $NCBI_tax2ictv_8_rank_tax{$ncbi_tax};
	print OUT "$id\t$ictv_8_rank_tax\n";
}
close OUT;

# Step 5. Write down novel rank
open OUT, ">NCBI_tax_novel_rank.txt";  # Store any novel ranks that are not present in the ICTV table
foreach my $novel_rank (sort keys %NCBI_tax_novel_rank){
	my $ncbi_tax = $NCBI_tax_novel_rank{$novel_rank};
	print OUT "$novel_rank\t$ncbi_tax\n";
}
close OUT;
