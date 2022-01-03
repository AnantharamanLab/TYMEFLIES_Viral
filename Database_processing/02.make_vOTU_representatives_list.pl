#!/usr/bin/perl

use strict;
use warnings;

# Aim: make vOTU representative list; pick the high quality and long phage genome

# Store all genome property
my %Gn_property = (); # $gn => [0] $vOTU, [1] $length, [2] $completeness, [3] $genome_quality, [4] $tax
my %VOTU2gn = (); # $vOTU => $gn separated by "\t";
open IN, "IMGVR_all_Sequence_information.tsv";
while (<IN>){
	chomp;
	if (!/^#/){
		my @tmp = split (/\t/);
		my $gn = $tmp[0];
		my $vOTU = $tmp[5];
		my $length = $tmp[6];
		my $completeness = $tmp[8];
		my $genome_quality = $tmp[9];
		my $tax = $tmp[11];
		$Gn_property{$gn}[0] = $vOTU;
		$Gn_property{$gn}[1] = $length;
		$Gn_property{$gn}[2] = $completeness;
		$Gn_property{$gn}[3] = $genome_quality;
		$Gn_property{$gn}[4] = $tax;
		if (!exists $VOTU2gn{$vOTU}){
			$VOTU2gn{$vOTU} = $gn;
		}else{
			$VOTU2gn{$vOTU} .= "\t".$gn;
		}
	}
}
close IN;

# Pick vOTU representative 
my %VOTU2rep = (); # $vOTU => $gn_rep (representative genome)
foreach my $vOTU (sort keys %VOTU2gn){
	my @Gn = split (/\t/, $VOTU2gn{$vOTU});
	
	# Divide genomes into 4 collections by genome quality
	my @Gn2 = (); # 2nd level genome collection for picking, containing "Reference"
	my @Gn3 = (); # 3rd level genome collection for picking, containing "High-quality"
	my @Gn4 = (); # 4th level genome collection for picking, containing "Genome fragment"
	my @Gn5 = (); # 5th level genome collection for picking, containing "Unsure (completeness > 120%)" or "Unsure (Ns)" or "Concatemer_artifact" 
	foreach my $gn (@Gn){
		if ($Gn_property{$gn}[3] eq "Reference"){
			push @Gn2, $gn;
		}elsif ($Gn_property{$gn}[3] eq "High-quality"){
			push @Gn3, $gn;
		}elsif ($Gn_property{$gn}[3] eq "Genome fragment"){
			push @Gn4, $gn;
		}elsif ($Gn_property{$gn}[3] =~ /Unsure/ or $Gn_property{$gn}[3] eq "Concatemer_artifact"){
			push @Gn5, $gn;
		}
	}
	
	if (@Gn2){ # If the 2nd level genome collection is not empty
		my $gn_rep = $Gn2[0]; 
		for(my $i=1; $i<=$#Gn2; $i++){
			if ($Gn_property{$Gn2[$i]}[1] > $Gn_property{$gn_rep}[1]){ # Pick the longest one from the 2nd level genome collection
				$gn_rep = $Gn2[$i];
			}
		}
		$VOTU2rep{$vOTU} = $gn_rep;
	}elsif (!@Gn2 and @Gn3){ # If the 2nd level genome collection is empty, and 3rd level genome collection is not empty
		my $gn_rep = $Gn3[0]; 
		for(my $i=1; $i<=$#Gn3; $i++){
			if ($Gn_property{$Gn3[$i]}[1] > $Gn_property{$gn_rep}[1]){ # Pick the longest one from the 3nd level genome collection
				$gn_rep = $Gn3[$i];
			}
		}
		$VOTU2rep{$vOTU} = $gn_rep;
	}elsif (!@Gn2 and !@Gn3 and @Gn4){ # If the 3rd level genome collection is empty, and 4th level genome collection is not empty
		my $gn_rep = $Gn4[0]; 
		for(my $i=1; $i<=$#Gn4; $i++){
			if ($Gn_property{$Gn4[$i]}[1] > $Gn_property{$gn_rep}[1]){ # Pick the longest one from the 4th level genome collection
				$gn_rep = $Gn4[$i];
			}
		}
		$VOTU2rep{$vOTU} = $gn_rep;
	}elsif (!@Gn2 and !@Gn3 and !@Gn4 and @Gn5){ # If the 4th level genome collection is empty, and 5th level genome collection is not empty
		my $gn_rep = $Gn5[0]; 
		for(my $i=1; $i<=$#Gn5; $i++){
			if ($Gn_property{$Gn5[$i]}[1] > $Gn_property{$gn_rep}[1]){ # Pick the longest one from the 5th level genome collection
				$gn_rep = $Gn5[$i];
			}
		}
		$VOTU2rep{$vOTU} = $gn_rep;		
	}
}

# Write down vOTU representatives table 

open OUT, ">IMGVR_all_phage_vOTU_representatives.txt";
print OUT "vOTU\tvOTU representative\tvOTU representative taxonomy\n";
foreach my $key (sort keys %VOTU2rep){
	print OUT "$key\t$VOTU2rep{$key}\t$Gn_property{$VOTU2rep{$key}}[4]\n";
}
close OUT;
close OUT;


# Make a new folder containing vOTU representatives
`mkdir vOTU_representatives_phage_genomes`;

open OUT, ">tmp.copy_genomes.sh";
foreach my $vOTU (sort keys %VOTU2rep){
	my $gn_rep = $VOTU2rep{$vOTU};
	
	# Replace "|" to "__" in the gn name
	my $gn_rep_new = $gn_rep; $gn_rep_new =~ s/\|/\_\_/g;
	print OUT "cp Individual_phage_genomes/$gn_rep_new.fasta vOTU_representatives_phage_genomes/$vOTU\_\_\_$gn_rep_new.fasta\n";
}
close OUT;

`cat tmp.copy_genomes.sh | parallel -j 20`;

`rm tmp.copy_genomes.sh`;




