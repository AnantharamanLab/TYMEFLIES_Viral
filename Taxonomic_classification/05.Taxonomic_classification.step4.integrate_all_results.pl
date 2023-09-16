#!/usr/bin/perl

use strict;
use warnings;

# AIM: Integrate all results and provide the final taxonomic classification result

# Step 1 Store taxonomic classification result from three methods
my %Viral_gn2tax = (); # $viral_gn => [0] $tax [1] $method 
open IN, "/storage1/data11/TYMEFLIES_phage/Taxonomic_classification/Each_bin_consensus_tax_by_NCBI_RefSeq_viral_protein_searching.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $viral_gn = $tmp[0];
	my $tax = $tmp[1];
	my $method = "NCBI RefSeq viral protein searching";
	# Change empty ranks into "NA"
	my @Tax = split (/\;/, $tax);
	for(my $i=0; $i<=$#Tax; $i++){
		if ($Tax[$i] eq ""){
			$Tax[$i] = "NA";
		}
	}		
	$tax = join("\;", @Tax);
	
	$Viral_gn2tax{$viral_gn}[0] = $tax;
	$Viral_gn2tax{$viral_gn}[1] = $method;
}
close IN;

open IN, "/storage1/data11/TYMEFLIES_phage/Taxonomic_classification/Each_bin_consensus_tax_by_VOG_marker_HMM_searching.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $viral_gn = $tmp[0];
	my $tax = $tmp[1];
	my $method = "VOG marker HMM searching";
	if (! exists $Viral_gn2tax{$viral_gn}){ # NCBI RefSeq viral protein searching method has the higher priority
		# Add genus and species into the tax (of course both are "NA")
		my @Tax = split (/\;/, $tax); 
		push @Tax, "NA";
		push @Tax, "NA";
		$tax = join("\;", @Tax);
		
		$Viral_gn2tax{$viral_gn}[0] = $tax;
		$Viral_gn2tax{$viral_gn}[1] = $method;
	}
}
close IN;

open IN, "/storage1/data11/TYMEFLIES_phage/Taxonomic_classification/genomad_output/all_virus_genome_fastas_Nlinked_annotate/all_virus_genome_fastas_Nlinked_taxonomy.tsv";
while (<IN>){
	chomp;
	if (!/^seq/){
		my @tmp = split (/\t/);
		my $viral_gn = $tmp[0];
		my $tax = $tmp[4];
		my $method = "geNomad classifying";
		if (! exists $Viral_gn2tax{$viral_gn}){ # NCBI RefSeq viral protein searching and VOG marker HMM searching methods have the higher priority
			# Add the rest ranks into the tax (of course both are "NA"), the total number of ranks is 8 (excluding "Virus")
			my @Tax = split (/\;/, $tax); 
			shift @Tax; # Delete the first element
			# Add "NA" elements
			while (@Tax < 8) {
				push @Tax, "NA";
			}
			$tax = join("\;", @Tax);
			
			$Viral_gn2tax{$viral_gn}[0] = $tax;
			$Viral_gn2tax{$viral_gn}[1] = $method;
		}
	}
}
close IN;

# Step 2 Store genus cluster map
my %Genus_cluster_map = (); # $gn_rep => $gns (all the genomes separated by ",")
open IN, "/storage1/data11/TYMEFLIES_phage/Cluster_phage_genomes/genus_clusters.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $gn_rep = $tmp[0];
	my $gns = join("\,", @tmp);
	$Genus_cluster_map{$gn_rep} = $gns;
}
close IN;

# Step 3 Get taxonomy based on other members' taxonomy from each genus
# Get into each genus cluster to see if any genomes have already got hits (only the hits of "NCBI RefSeq viral protein searching" will be counted), 
# then expand the tax to all the members within this genus cluster
my %Viral_gn2tax_by_genus_cluster = (); # $viral_gn => [0] $tax [1] $method (This hash stores taxonomy result assigned by the "Genus-level vOTU LCA assigning")
foreach my $gn_rep (sort keys %Genus_cluster_map){
	my @Gns = split (/\,/,$Genus_cluster_map{$gn_rep});
	
	my @Tax = (); # The taxonomy collection for all the genomes within this genus cluster
	foreach my $gn (@Gns){
		if (exists $Viral_gn2tax{$gn}){
			if ($Viral_gn2tax{$gn}[1] eq "NCBI RefSeq viral protein searching"){
				push @Tax, $Viral_gn2tax{$gn}[0];
			}
		}
	}
	
	my $lca = 'NA;NA;NA;NA;NA;NA;NA;NA'; # The LCA for all taxonomic hits within this genus (should be above genus level)
	if (@Tax){
		$lca = _get_LCA_from_vOTU(@Tax);
		if ($lca ne 'NA;NA;NA;NA;NA;NA;NA;NA'){
			foreach my $gn (@Gns){
				if (!exists $Viral_gn2tax{$gn}){
					$Viral_gn2tax_by_genus_cluster{$gn}[0] = $lca;
					$Viral_gn2tax_by_genus_cluster{$gn}[1] = "Genus-level vOTU LCA assigning";
				}
			}
		}
	}
}

# Step 4 Write down final taxonomic classification result
open OUT, ">/storage1/data11/TYMEFLIES_phage/Taxonomic_classification/Each_bin_tax_combined_result.txt";
foreach my $gn (sort keys %Viral_gn2tax){
	print OUT "$gn\t$Viral_gn2tax{$gn}[0]\t$Viral_gn2tax{$gn}[1]\n";
}
foreach my $gn (sort keys %Viral_gn2tax_by_genus_cluster){
	print OUT "$gn\t$Viral_gn2tax_by_genus_cluster{$gn}[0]\t$Viral_gn2tax_by_genus_cluster{$gn}[1]\n";
}
close OUT;



# Subroutine

sub _get_LCA_from_vOTU{
	my @Taxonomy = @_; # Get the passed array
	
	my %Realm = (); # $realm => the times of this realm appears
	my %Kingdom = (); # $kingdom => the times of this kingdom appears
	my %Phylum = (); # $phylum => the times of this phylum appears
	my %Class = (); # $class => the times of this class appears
	my %Order = (); # $order => the times of this order appears
	my %Family = (); # $family => the times of this family appears
	my %Genus = (); # $genus => the times of this genus appears
	my %Species = (); # $species => the times of this species appears
	
	foreach my $tax (@Taxonomy){
		my @tmp = split (/\;/, $tax);
		$Realm{$tmp[0]}++;
		$Kingdom{$tmp[1]}++;
		$Phylum{$tmp[2]}++;
		$Class{$tmp[3]}++;
		$Order{$tmp[4]}++;
		$Family{$tmp[5]}++;
		$Genus{$tmp[6]}++;
		$Species{$tmp[7]}++;
	}
	
	my @LCA_final = ();
	OUTTER: {if (scalar (keys %Realm) == 1){
				my @Realm = keys %Realm;
				push @LCA_final, $Realm[0];
				if (scalar (keys %Kingdom) == 1){
					my @Kingdom = keys %Kingdom;
					push @LCA_final, $Kingdom[0];	
					if (scalar (keys %Phylum) == 1){
						my @Phylum = keys %Phylum;
						push @LCA_final, $Phylum[0];
						if (scalar (keys %Class) == 1){
							my @Class = keys %Class;
							push @LCA_final, $Class[0];
							if (scalar (keys %Order) == 1){
								my @Order = keys %Order;
								push @LCA_final, $Order[0];
								if (scalar (keys %Family) == 1){
									my @Family = keys %Family;
									push @LCA_final, $Family[0];
									if (scalar (keys %Genus) == 1){
										my @Genus = keys %Genus;
										push @LCA_final, $Genus[0];
									}else{
										last OUTTER;
									}
								}else{
									last OUTTER;
								}
							}else{
								last OUTTER;
							}
						}else{
							last OUTTER;
						}
					}else{
						last OUTTER;
					}
				}else{
					last OUTTER;
				}
			}else{
				last OUTTER;
			}
	}
	
	my $lca_final_full = 'NA;NA;NA;NA;NA;NA;NA;NA';
	my @LCA_final_full = split (/\;/, $lca_final_full);
	
	if (@LCA_final){
		for(my $i=0; $i<=$#LCA_final; $i++){
			$LCA_final_full[$i] = $LCA_final[$i];
		}
	}	
	
	$lca_final_full = join("\;", @LCA_final_full);
	
	return $lca_final_full;
}




