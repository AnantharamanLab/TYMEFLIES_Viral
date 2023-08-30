#!/usr/bin/perl

use strict;
use warnings;

# AIM: Integrate all host prediction results and write down the final result

# Step 1 Store host prediction results from three methods
## Step 1.1 Store results from prophage
my %Prophage_gn2host = (); # $gn => $host_tax
open IN, "/storage1/data11/TYMEFLIES_phage/Host_prediction/Prophage_gn2host.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $gn = $tmp[0];
	my $host_tax = $tmp[2];
	$host_tax = _make_up_full_ranks($host_tax);
	if ($host_tax !~ /\;f\_\_\;/){ # If the host tax is above family level then do not use
		$Prophage_gn2host{$gn} = $host_tax; 
	}
}
close IN;

## Step 1.2 Store results from Phage_gn2host_tax_based_on_AMG.txt
my %Phage_gn2host_tax_based_on_AMG = (); # $gn => $host_tax
open IN, "/storage1/data11/TYMEFLIES_phage/Host_prediction/Phage_gn2host_tax_based_on_AMG.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $gn = $tmp[0];
	my $host_tax = $tmp[1];
	$host_tax = _make_up_full_ranks($host_tax);
	if ($host_tax !~ /\;f\_\_\;/){ # If the host tax is above family level then do not use
		$Phage_gn2host_tax_based_on_AMG{$gn} = $host_tax; 
	}
}
close IN;

## Step 1.3 Store results from iPHoP prediction
my %Phage_gn2host_tax_based_on_iPHoP = (); # $gn => $host_tax
open IN, "/storage1/data11/TYMEFLIES_phage/Host_prediction/iPHoP_result/iPHoP_result_parsed.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $gn = $tmp[0];
	my $host_tax = $tmp[1];
	$host_tax = _make_up_full_ranks($host_tax);
	if ($host_tax !~ /\;f\_\_\;/){ # If the host tax is above family level then do not use
		$Phage_gn2host_tax_based_on_iPHoP{$gn} = $host_tax; 
	}	
}
close IN;

# Step 2 Store species cluster map
my %Species_cluster_map = (); # $gn_rep => $gns (all the genomes separated by ",")
open IN, "/storage1/data11/TYMEFLIES_phage/Cluster_phage_genomes/Species_level_vOTUs_cluster.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $gn_rep = $tmp[0];
	my $gns = $tmp[1];
	$Species_cluster_map{$gn_rep} = $gns;
}
close IN;

# Step 3 Get host prediction based on other members' host prediction from each species
# Get into each species cluster to see if any genomes have already got hits, then expand the host prediction to all the members within this species cluster
# Use prophage, AMG host prediction, and iPHoP prediction results to get other vOTU member host prediction
# Only use the result that is down to the family level
my %Viral_gn2host_tax_by_species_cluster = (); # $viral_gn => $host_tax
foreach my $gn_rep (sort keys %Species_cluster_map){
	my @Gns = split (/\,/,$Species_cluster_map{$gn_rep});
	
	my @Tax = (); # The host tax collection for all the genomes within this species cluster
	foreach my $gn (@Gns){
		if (exists $Prophage_gn2host{$gn}){
			push @Tax, $Prophage_gn2host{$gn};
		}elsif(exists $Phage_gn2host_tax_based_on_AMG{$gn}){
			push @Tax, $Phage_gn2host_tax_based_on_AMG{$gn};
		}elsif(exists $Phage_gn2host_tax_based_on_iPHoP{$gn}){
			push @Tax, $Phage_gn2host_tax_based_on_iPHoP{$gn};
		}			
	}
	
	my $lca = ""; # The LCA for all host tax hits within this species cluster
	if (@Tax){
		$lca = _get_LCA_from_vOTU(@Tax);
		if ($lca and $lca !~ /\;f\_\_\;/){ # If the $lca is above family level then do not use
			foreach my $gn (@Gns){
				if (!exists $Prophage_gn2host{$gn} and !exists $Phage_gn2host_tax_based_on_AMG{$gn} and !exists $Phage_gn2host_tax_based_on_iPHoP{$gn}){
					$Viral_gn2host_tax_by_species_cluster{$gn} = $lca;
				}
			}
		}
	}
}

# Step 4 Integrate all results
# The overlapped host taxonomies were solved based on the following priority: 
# 1) prophage within a host genome; 
# 2) AMG match to host genome;
# 3) iPHoP result;
# 4) derived from vOTU host taxonomy. 

my %Viral_gn2host_tax_all_combined = (%Prophage_gn2host, %Phage_gn2host_tax_based_on_AMG, %Phage_gn2host_tax_based_on_iPHoP, %Viral_gn2host_tax_by_species_cluster);

my %Viral_gn2host_tax_final = (); # $viral_gn => [0] $host_tax [1] $method
foreach my $gn (sort keys %Viral_gn2host_tax_all_combined){
	my $array_host_tax = 'NA;NA;NA;NA';
	my @Array_host_tax = split (/\;/, $array_host_tax);
	my @Methods = ();
	$Methods[0] = "prophage within a host genome";
	$Methods[1] = "AMG match to a host genome";
	$Methods[2] = "iPHoP result";
	$Methods[3] = "derived from vOTU host taxonomy";
	
	if (exists $Prophage_gn2host{$gn}){
		$Array_host_tax[0] = $Prophage_gn2host{$gn};
	}
	
	if (exists $Phage_gn2host_tax_based_on_AMG{$gn}){
		$Array_host_tax[1] = $Phage_gn2host_tax_based_on_AMG{$gn};
	}

	if (exists $Phage_gn2host_tax_based_on_iPHoP{$gn}){
		$Array_host_tax[2] = $Phage_gn2host_tax_based_on_iPHoP{$gn};
	}	
	
	if (exists $Viral_gn2host_tax_by_species_cluster{$gn}){
		$Array_host_tax[3] = $Viral_gn2host_tax_by_species_cluster{$gn};
	}
	
	my $host_tax_final = "";
	my $method = "";
	
	for(my $i=0; $i<=$#Array_host_tax; $i++){
		if ($Array_host_tax[$i] ne "NA"){
			$host_tax_final = $Array_host_tax[$i];
			$method = $Methods[$i];
			last;
		}
	}
	
	$Viral_gn2host_tax_final{$gn}[0] = $host_tax_final;
	$Viral_gn2host_tax_final{$gn}[1] = $method;
}

# Step 5 Write down result
open OUT, ">Host_prediction/Viral_gn2host_tax_final.txt";
foreach my $gn (sort keys %Viral_gn2host_tax_final){
	print OUT "$gn\t$Viral_gn2host_tax_final{$gn}[0]\t$Viral_gn2host_tax_final{$gn}[1]\n";
}
close OUT;

# Subroutine

sub _get_LCA_from_vOTU{
	my @Taxonomy = @_; # Get the passed array

	my %Domain = (); # $domain => the times of this domain appears
	my %Phylum = (); # $phylum => the times of this phylum appears
	my %Class = (); # $class => the times of this class appears
	my %Order = (); # $order => the times of this order appears
	my %Family = (); # $family => the times of this family appears
	my %Genus = (); # $genus => the times of this genus appears
	my %Species = (); # $species => the times of this species appears
	
	foreach my $tax (@Taxonomy){
		my @tmp = split (/\;/, $tax);
		$Domain{$tmp[0]}++;
		$Phylum{$tmp[1]}++;
		$Class{$tmp[2]}++;
		$Order{$tmp[3]}++;
		$Family{$tmp[4]}++;
		$Genus{$tmp[5]}++;
		$Species{$tmp[6]}++;
	}
	
	my @LCA_final = ();
	OUTTER: {if (scalar (keys %Domain) == 1){
				my @Domain = keys %Domain;
				push @LCA_final, $Domain[0];
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
									if (scalar (keys %Species) == 1){
										my @Species = keys %Species;
										push @LCA_final, $Species[0];
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
	
	my $lca_final_full = 'd__;p__;c__;o__;f__;g__;s__';
	my @LCA_final_full = split (/\;/, $lca_final_full);
	
	if (@LCA_final){
		for(my $i=0; $i<=$#LCA_final; $i++){
			$LCA_final_full[$i] = $LCA_final[$i];
		}
	}	
	
	$lca_final_full = join("\;", @LCA_final_full);
	
	return $lca_final_full;
}

sub _make_up_full_ranks{
	my $input_tax = $_[0];
	my $output_tax = "";
	
	my @Input_tax_ranks = split (/\;/, $input_tax);
	my $domain = "d__";
	my $phylum = "p__";
	my $class = "c__";
	my $order = "o__";
	my $family = "f__";
	my $genus = "g__";
	my $species = "s__";
	for(my $i=0; $i<=$#Input_tax_ranks; $i++){
		if ($Input_tax_ranks[$i] =~ /d\_\_/){
			$domain = $Input_tax_ranks[$i]; 
		}elsif ($Input_tax_ranks[$i] =~ /p\_\_/){
			$phylum = $Input_tax_ranks[$i]; 
		}elsif ($Input_tax_ranks[$i] =~ /c\_\_/){
			$class = $Input_tax_ranks[$i]; 
		}elsif ($Input_tax_ranks[$i] =~ /o\_\_/){
			$order = $Input_tax_ranks[$i]; 
		}elsif ($Input_tax_ranks[$i] =~ /f\_\_/){
			$family = $Input_tax_ranks[$i]; 
		}elsif ($Input_tax_ranks[$i] =~ /g\_\_/){
			$genus = $Input_tax_ranks[$i]; 
		}elsif ($Input_tax_ranks[$i] =~ /s\_\_/){
			$species = $Input_tax_ranks[$i]; 
		}
	}
	
	$output_tax = "$domain;$phylum;$class;$order;$family;$genus;$species";
	return $output_tax;
}