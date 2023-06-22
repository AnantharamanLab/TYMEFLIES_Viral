#!/usr/bin/env perl

use strict;
use warnings;
use Array::Split qw( split_into );

# AIM: Investigate AMG variation in each species and genus (non-singleton)

# Step 1 Grep all AMG proteins
my %AMG_pro2ko = (); # $amg_pro => $ko
my %Gn2ko2num = (); # $gn => $ko => number of hits
my %Gn2ko = (); # $gn => $ko collections
my %KO_all = (); # $ko => 1
my $header = ""; # Store the header of AMG_summary.txt
open IN, "AMG_analysis/AMG_summary.txt";
while (<IN>){
	chomp;
	if (/^Pro/){
		$header = $_;
	}else{	
		my $line = $_;
		my @tmp = split (/\t/, $line);
		my $amg_pro = $tmp[0];
		my $ko = $tmp[2]; 
		$AMG_pro2ko{$amg_pro} = $ko;
		
		my ($gn) = $amg_pro =~ /^(.+?\_\_.+?)\_\_/;
		$Gn2ko2num{$gn}{$ko}++;
		
		if (!exists $Gn2ko{$gn}){
			$Gn2ko{$gn} = $ko;
		}else{
			$Gn2ko{$gn} .= "\t".$ko;
		}
		
		$KO_all{$ko} = 1;
	}
}
close IN;

# Step 2 Store and write down non-singleton species-level vOTU and check AMG variation
## Step 2.1 Store AMG variation
my %Species_level_vOTU = (); # $gn_rep => [0] $gns [1] $num_gns
open IN, "Cluster_phage_genomes/Species_level_vOTUs_cluster.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $gns = $tmp[1];
	my @Gns = split (/\,/, $gns);
	if ((scalar @Gns) >= 2){
		my $gn_rep = $tmp[0];
		$Species_level_vOTU{$gn_rep}[0] = $gns;
		$Species_level_vOTU{$gn_rep}[1] = scalar @Gns;
	}
}
close IN;

my %Species_level_vOTU2ko2freq = (); # $gn_rep => $ko => $freq_ko_in_all_gns
                                     # species-level vOTUs that do not have any ko hits are not included in this hash
									 
my %Species_level_vOTU_with_ko = (); # $gn_rep => 1
my %Species_level_vOTU_with_high_freq_in_all_gns = (); # $gn_rep => 1
my %Species_level_vOTU2kos = (); # $gn_rep => $ko collection
foreach my $gn_rep (sort keys %Species_level_vOTU){
	my @Gns = split (/\,/, $Species_level_vOTU{$gn_rep}[0]);
	
	my %KOs = (); # Store all the AMG KOs in this species-level vOTU
	foreach my $gn (@Gns){
		if (exists $Gn2ko{$gn}){
			my $kos = $Gn2ko{$gn};
			my @tmp = split (/\t/, $kos);
			foreach my $key (@tmp){
				$KOs{$key} = 1;
			}
		}
	}
	
	if (%KOs){
		$Species_level_vOTU_with_ko{$gn_rep} = 1;
		$Species_level_vOTU2kos{$gn_rep} = join("\t", (keys %KOs));
	}
	
	foreach my $ko (sort keys %KOs){
		my $freq_count_for_each_ko = 0;
		foreach my $gn (@Gns){
			if ($Gn2ko2num{$gn}{$ko}){
				$freq_count_for_each_ko++;
			}
		}
		my $freq_ko_in_all_gns = $freq_count_for_each_ko / (scalar @Gns);
		
		if ($freq_ko_in_all_gns){
			$Species_level_vOTU2ko2freq{$gn_rep}{$ko} = $freq_ko_in_all_gns;
		}
	}
	
	# Store vOTU that has high frequency of AMG KOs in all its genomes
	# cutoff =  0.75
	my $logic = 1;
	foreach my $ko (sort keys %KOs){
		if ($Species_level_vOTU2ko2freq{$gn_rep}{$ko} <= 0.75){ # If any ko freq is below the cutoff - 0.75
			$logic = 0;
		}
	}
	
	if ($logic and %KOs){
		$Species_level_vOTU_with_high_freq_in_all_gns{$gn_rep} = 1;
	}
}

## Step 2.2 Write down AMG variation result
open OUT, ">AMG_analysis/Species_level_vOTU_AMG_variation.txt";
my $num_nonsingleton_vOTU = scalar (keys %Species_level_vOTU);
print OUT "The number of non-singleton vOTU: $num_nonsingleton_vOTU\n";
my $num_vOTU_with_ko = scalar (keys %Species_level_vOTU_with_ko);
print OUT "The number of non-singleton vOTU with any AMG KO hit(s): $num_vOTU_with_ko\n";
my $num_vOTU_with_high_freq_in_all_gns = scalar (keys %Species_level_vOTU_with_high_freq_in_all_gns);
print OUT "The number of non-singleton vOTU with high AMG frequency (> 75%) in its all genomes: $num_vOTU_with_high_freq_in_all_gns\n";
close OUT;

## Step 2.3 Make statistics for AMG variation
### Step 2.3.1 Four quartile abundance 
my %VOTU_n_ko_freq_category2abundance = (); # $vOTU_n_ko_freq_category (four quartile: 1st quartile, 2nd quartile, 3rd quartile, and 4th quartile) => $abundance (in percentage)
my @VOTU_n_ko_freq_collection = (); # Store all the freq collections of vOTU and ko combinations
foreach my $gn_rep (sort keys %Species_level_vOTU_with_ko){
	my $kos = $Species_level_vOTU2kos{$gn_rep};
	my @KOs = split (/\t/, $kos);
	foreach my $ko (@KOs){
		my $vOTU_n_ko_freq = $Species_level_vOTU2ko2freq{$gn_rep}{$ko};
		push @VOTU_n_ko_freq_collection,$vOTU_n_ko_freq;
	}
}

%VOTU_n_ko_freq_category2abundance = _get_four_quartiles_abundance(@VOTU_n_ko_freq_collection);

#### Write down result
open OUT, ">AMG_analysis/Species_level_vOTU_AMG_variation_statistics.txt";
print OUT "The total number of vOTU and ko combinations (for vOTUs with any AMG KO hits(s)): ".(scalar @VOTU_n_ko_freq_collection)."\n";
foreach my $key (sort keys %VOTU_n_ko_freq_category2abundance){
	print OUT "$key\t$VOTU_n_ko_freq_category2abundance{$key}\n";
}
close OUT;

### Step 2.3.2 Four quartiles abundance for different species-level vOTU size
#### Split into four vOTU size
my %VOTU_w_ko2size = (); # $vOTU => $size (non-singleton vOTU with any AMG KO hit(s))
foreach my $vOTU (sort keys %Species_level_vOTU_with_ko){
	$VOTU_w_ko2size{$vOTU} = $Species_level_vOTU{$vOTU}[1]; 
}

my @VOTU_w_ko = (); # Sorted species-level vOTU
foreach my $key ( sort { $VOTU_w_ko2size{$a} <=> $VOTU_w_ko2size{$b} } keys %VOTU_w_ko2size ) {
 	push @VOTU_w_ko, $key;
}

my @Array_refs = split_into(4,@VOTU_w_ko);

my $array_ref_1 = $Array_refs[0];my $array_ref_2 = $Array_refs[1];
my $array_ref_3 = $Array_refs[2];my $array_ref_4 = $Array_refs[3];

my @VOTU_w_ko_1st_quartile = @$array_ref_1;my @VOTU_w_ko_2nd_quartile = @$array_ref_2;
my @VOTU_w_ko_3rd_quartile = @$array_ref_3;my @VOTU_w_ko_4th_quartile = @$array_ref_4;

my $vOTU_w_ko_1st_quartile_size_range = $Species_level_vOTU{$VOTU_w_ko_1st_quartile[0]}[1]."\-".$Species_level_vOTU{$VOTU_w_ko_1st_quartile[-1]}[1];
my $vOTU_w_ko_2nd_quartile_size_range = $Species_level_vOTU{$VOTU_w_ko_2nd_quartile[0]}[1]."\-".$Species_level_vOTU{$VOTU_w_ko_2nd_quartile[-1]}[1];
my $vOTU_w_ko_3rd_quartile_size_range = $Species_level_vOTU{$VOTU_w_ko_3rd_quartile[0]}[1]."\-".$Species_level_vOTU{$VOTU_w_ko_3rd_quartile[-1]}[1];
my $vOTU_w_ko_4th_quartile_size_range = $Species_level_vOTU{$VOTU_w_ko_4th_quartile[0]}[1]."\-".$Species_level_vOTU{$VOTU_w_ko_4th_quartile[-1]}[1];

#### Calculate four quartile abundance for each array (species-level vOTUs with different size)
my %VOTU_n_ko_freq_category2abundance_for_1st_quartile = ();
my @VOTU_n_ko_freq_collection_for_1st_quartile = (); # Store all the freq collections of vOTU and ko combinations for 1st quartile size species-level vOTUs
foreach my $gn_rep (@VOTU_w_ko_1st_quartile){
	my $kos = $Species_level_vOTU2kos{$gn_rep};
	my @KOs = split (/\t/, $kos);
	foreach my $ko (@KOs){
		my $vOTU_n_ko_freq = $Species_level_vOTU2ko2freq{$gn_rep}{$ko};
		push @VOTU_n_ko_freq_collection_for_1st_quartile,$vOTU_n_ko_freq;
	}
}
%VOTU_n_ko_freq_category2abundance_for_1st_quartile = _get_four_quartiles_abundance(@VOTU_n_ko_freq_collection_for_1st_quartile);

my %VOTU_n_ko_freq_category2abundance_for_2nd_quartile = ();
my @VOTU_n_ko_freq_collection_for_2nd_quartile = (); # Store all the freq collections of vOTU and ko combinations for 2nd quartile size species-level vOTUs
foreach my $gn_rep (@VOTU_w_ko_2nd_quartile){
	my $kos = $Species_level_vOTU2kos{$gn_rep};
	my @KOs = split (/\t/, $kos);
	foreach my $ko (@KOs){
		my $vOTU_n_ko_freq = $Species_level_vOTU2ko2freq{$gn_rep}{$ko};
		push @VOTU_n_ko_freq_collection_for_2nd_quartile,$vOTU_n_ko_freq;
	}
}
%VOTU_n_ko_freq_category2abundance_for_2nd_quartile = _get_four_quartiles_abundance(@VOTU_n_ko_freq_collection_for_2nd_quartile);

my %VOTU_n_ko_freq_category2abundance_for_3rd_quartile = ();
my @VOTU_n_ko_freq_collection_for_3rd_quartile = (); # Store all the freq collections of vOTU and ko combinations for 3rd quartile size species-level vOTUs
foreach my $gn_rep (@VOTU_w_ko_3rd_quartile){
	my $kos = $Species_level_vOTU2kos{$gn_rep};
	my @KOs = split (/\t/, $kos);
	foreach my $ko (@KOs){
		my $vOTU_n_ko_freq = $Species_level_vOTU2ko2freq{$gn_rep}{$ko};
		push @VOTU_n_ko_freq_collection_for_3rd_quartile,$vOTU_n_ko_freq;
	}
}
%VOTU_n_ko_freq_category2abundance_for_3rd_quartile = _get_four_quartiles_abundance(@VOTU_n_ko_freq_collection_for_3rd_quartile);

my %VOTU_n_ko_freq_category2abundance_for_4th_quartile = ();
my @VOTU_n_ko_freq_collection_for_4th_quartile = (); # Store all the freq collections of vOTU and ko combinations for 4th quartile size species-level vOTUs
foreach my $gn_rep (@VOTU_w_ko_4th_quartile){
	my $kos = $Species_level_vOTU2kos{$gn_rep};
	my @KOs = split (/\t/, $kos);
	foreach my $ko (@KOs){
		my $vOTU_n_ko_freq = $Species_level_vOTU2ko2freq{$gn_rep}{$ko};
		push @VOTU_n_ko_freq_collection_for_4th_quartile,$vOTU_n_ko_freq;
	}
}
%VOTU_n_ko_freq_category2abundance_for_4th_quartile = _get_four_quartiles_abundance(@VOTU_n_ko_freq_collection_for_4th_quartile);

#### Write down result
open OUT, ">>AMG_analysis/Species_level_vOTU_AMG_variation_statistics.txt";
print OUT "The total number of vOTU and ko combinations (for vOTUs with any AMG KO hits(s)) of the 1st quartile size species-level vOTUs: ".(scalar @VOTU_n_ko_freq_collection_for_1st_quartile)."\n";
print OUT "The size range of 1st quartile size species-level vOTUs: $vOTU_w_ko_1st_quartile_size_range\n";
foreach my $key (sort keys %VOTU_n_ko_freq_category2abundance_for_1st_quartile){
	print OUT "$key\t$VOTU_n_ko_freq_category2abundance_for_1st_quartile{$key}\n";
}

print OUT "The total number of vOTU and ko combinations (for vOTUs with any AMG KO hits(s)) of the 2nd quartile size species-level vOTUs: ".(scalar @VOTU_n_ko_freq_collection_for_2nd_quartile)."\n";
print OUT "The size range of 2nd quartile size species-level vOTUs: $vOTU_w_ko_2nd_quartile_size_range\n";
foreach my $key (sort keys %VOTU_n_ko_freq_category2abundance_for_2nd_quartile){
	print OUT "$key\t$VOTU_n_ko_freq_category2abundance_for_2nd_quartile{$key}\n";
}

print OUT "The total number of vOTU and ko combinations (for vOTUs with any AMG KO hits(s)) of the 3rd quartile size species-level vOTUs: ".(scalar @VOTU_n_ko_freq_collection_for_3rd_quartile)."\n";
print OUT "The size range of 3rd quartile size species-level vOTUs: $vOTU_w_ko_3rd_quartile_size_range\n";
foreach my $key (sort keys %VOTU_n_ko_freq_category2abundance_for_3rd_quartile){
	print OUT "$key\t$VOTU_n_ko_freq_category2abundance_for_3rd_quartile{$key}\n";
}

print OUT "The total number of vOTU and ko combinations (for vOTUs with any AMG KO hits(s)) of the 4th quartile size species-level vOTUs: ".(scalar @VOTU_n_ko_freq_collection_for_4th_quartile)."\n";
print OUT "The size range of 4th quartile size species-level vOTUs: $vOTU_w_ko_4th_quartile_size_range\n";
foreach my $key (sort keys %VOTU_n_ko_freq_category2abundance_for_4th_quartile){
	print OUT "$key\t$VOTU_n_ko_freq_category2abundance_for_4th_quartile{$key}\n";
}
close OUT;

### Step 2.3.3 Get KO distribution from the 1st quartile ko frequency (> 75%) and the 4th quartile in bin size (@VOTU_w_ko_4th_quartile)
### and make the KO2ko_abun_n_mean_ko_freq hash
my %KO_distribution_from_1st_quartile_ko_frequency = (); # $ko => $abundance (in percentage)
my $ko_total_num = 0;
my %KO2vOTU_collection = (); # $ko => $vOTU collection separated by "\t"
foreach my $gn_rep (@VOTU_w_ko_4th_quartile){
	my $kos = $Species_level_vOTU2kos{$gn_rep};
	my @KOs = split (/\t/, $kos);
	foreach my $ko (@KOs){
		my $vOTU_n_ko_freq = $Species_level_vOTU2ko2freq{$gn_rep}{$ko};
		if ($vOTU_n_ko_freq > 0.75){
			$ko_total_num++;
			$KO_distribution_from_1st_quartile_ko_frequency{$ko}++;
			if (!exists $KO2vOTU_collection{$ko}){ 
				$KO2vOTU_collection{$ko} = $gn_rep;
			}else{
				$KO2vOTU_collection{$ko} .= "\t".$gn_rep;
			}
		}
	}
}

foreach my $ko (sort keys %KO_distribution_from_1st_quartile_ko_frequency){
	my $abun = $KO_distribution_from_1st_quartile_ko_frequency{$ko};
	$abun = $abun / $ko_total_num;
	$KO_distribution_from_1st_quartile_ko_frequency{$ko} = $abun;
}

my @KO_list_of_KO_distribution_from_1st_quartile_ko_frequency = ();
open OUT, ">AMG_analysis/Species_level_vOTU_AMG_variation_statistics.2.txt";
print OUT "The total number of KO for AMG KO frequency > 75%: ".$ko_total_num."\n";
my $ko_distribtion_abun_sum = 0;
foreach my $ko (sort keys %KO_distribution_from_1st_quartile_ko_frequency){
	if ($KO_distribution_from_1st_quartile_ko_frequency{$ko} >= 0){
		print OUT "$ko\t$KO_distribution_from_1st_quartile_ko_frequency{$ko}\n";
		$ko_distribtion_abun_sum += $KO_distribution_from_1st_quartile_ko_frequency{$ko};
		push @KO_list_of_KO_distribution_from_1st_quartile_ko_frequency, $ko;
	}
}
print OUT "The rest KO\t".(1 - $ko_distribtion_abun_sum)."\n";
close OUT;

my %KO2ko_abun_n_mean_ko_freq = (); # $ko => $ko_abun and $mean_ko_freq combined by "\t"
foreach my $ko (sort keys %KO_distribution_from_1st_quartile_ko_frequency){
	my $ko_abun = $KO_distribution_from_1st_quartile_ko_frequency{$ko};
	
	my $vOTU_num = 0;
	
	my @VOTU_collection = split (/\t/, $KO2vOTU_collection{$ko});
	my @KO_freq = ();
	foreach my $vOTU (@VOTU_collection){
		push @KO_freq, $Species_level_vOTU2ko2freq{$vOTU}{$ko};
	}
	
	my $mean_ko_freq = 0;
	$mean_ko_freq = _avg(@KO_freq);
	
	$KO2ko_abun_n_mean_ko_freq{$ko} = "$ko_abun\t$mean_ko_freq";
}

open OUT, ">AMG_analysis/KO2ko_abun_n_mean_ko_freq.txt";
foreach my $ko (sort keys %KO2ko_abun_n_mean_ko_freq){
	print OUT "$ko\t$KO2ko_abun_n_mean_ko_freq{$ko}\n";
}
close OUT;

### Step 2.3.4 Get the sample date distribution for KO-carrying viruses from the 1st quartile ko frequency (> 75%) and the 4th quartile in bin size (@VOTU_w_ko_4th_quartile)
#### Store the sample date (date in a year) to each metagenome
my %IMG_ID2date_in_a_year = ();
open IN, "TYMEFLIES_metagenome_info.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		my @tmp = split (/\t/);
		my $img_id = $tmp[0];
		my $date_in_a_year = $tmp[9];
		$IMG_ID2date_in_a_year{$img_id} = $date_in_a_year;
	}
}	
close IN;

my @Dates_in_a_year = (1..365); # Store all the dates in a year

#### Store KO to dates in a year collection
my %KO2dates_in_a_year = (); # $ko => $dates_in_a_year (collection of $date_in_a_year, separated by "\t")
foreach my $ko (sort keys %KO2vOTU_collection){
	my %Dates_in_a_year_for_this_ko = (); # Store all $date_in_a_year
	
	my @VOTUs = split (/\t/, $KO2vOTU_collection{$ko});
	foreach my $gn_rep (@VOTUs){
		my @Gns = split (/\,/, $Species_level_vOTU{$gn_rep}[0]);
		foreach my $gn (@Gns){
			my ($img_id) = $gn =~ /^(.+?)\_\_/;
			my $date_in_a_year = $IMG_ID2date_in_a_year{$img_id};
			$Dates_in_a_year_for_this_ko{$date_in_a_year} = 1;
		}
	}
	
	my $dates_in_a_year = join ("\t", (sort keys %Dates_in_a_year_for_this_ko));
	$KO2dates_in_a_year{$ko} = $dates_in_a_year;
}

#### Make KO to date in a year to presence hash
my %KO2date_in_a_year2presence = (); # $ko => $date_in_a_year => $presence
foreach my $ko (sort keys %KO2vOTU_collection){
	foreach my $date_in_a_year (@Dates_in_a_year){
		my $presence = 0;
		my @Dates_in_a_year =  split (/\t/, $KO2dates_in_a_year{$ko});
		my %Dates_in_a_year = map { $_ => 1 } @Dates_in_a_year;
		if (exists $Dates_in_a_year{$date_in_a_year}){
			$presence = 1;
		}
		$KO2date_in_a_year2presence{$ko}{$date_in_a_year} = $presence;
	}
}

#### Write down %KO2date_in_a_year2presence hash
open OUT, ">AMG_analysis/KO2dates_in_a_year.txt";
my $row=join("\t", @Dates_in_a_year);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %KO2date_in_a_year2presence){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (@Dates_in_a_year) {       
                if (exists $KO2date_in_a_year2presence{$tmp1}{$tmp2}){
                        push @tmp, $KO2date_in_a_year2presence{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;




# Subroutine
sub _store_seq{
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

sub _get_four_quartiles_abundance{
	my @list = @_;
	my %Quartile2abundance = (); 
	foreach my $key (@list){
		if ($key > 0.75 and $key <= 1){
			$Quartile2abundance{"1st quartile"}++;
		}elsif($key > 0.5 and $key <= 0.75){
			$Quartile2abundance{"2nd quartile"}++;
		}elsif($key > 0.25 and $key <= 0.5){
			$Quartile2abundance{"3rd quartile"}++;
		}elsif($key > 0 and $key <= 0.25){
			$Quartile2abundance{"4th quartile"}++;
		}	
	}

	foreach my $quartile (sort keys %Quartile2abundance){
		my $abun = $Quartile2abundance{$quartile};
		$abun = $abun / (scalar @list);
		$Quartile2abundance{$quartile} = $abun;
	}
	
	return %Quartile2abundance;
}	

sub _avg {
    my $total;
    $total += $_ foreach @_;
    # sum divided by number of components.
    return $total / @_;
}



