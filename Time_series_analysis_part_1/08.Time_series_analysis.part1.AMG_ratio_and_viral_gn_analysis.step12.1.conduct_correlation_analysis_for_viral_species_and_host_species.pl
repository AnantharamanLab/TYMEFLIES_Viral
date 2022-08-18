#!/usr/bin/perl

use strict;
use warnings;

# Aim: Сonduct correlation analysis for viral species and the corresponding host species

# Note: This script should be run under env vRhyme (conda activate vRhyme)

# Step 1 Store viral host pair information
my %Viral_host_pair = ();
$Viral_host_pair{"3300042345__vRhyme_unbinned588|3300042540__vRhyme_272"} = "d__Bacteria;p__Cyanobacteria;c__Cyanobacteriia;o__Cyanobacteriales;f__Microcoleaceae;g__Planktothrix;s__Planktothrix agardhii";
$Viral_host_pair{"3300042525__vRhyme_unbinned1285"} = "d__Bacteria;p__Cyanobacteria;c__Cyanobacteriia;o__Pseudanabaenales;f__Pseudanabaenaceae;g__Pseudanabaena;s__Pseudanabaena sp014696925";
$Viral_host_pair{"3300042865__vRhyme_115|3300042433__vRhyme_19|3300044641__vRhyme_unbinned881|3300043567__vRhyme_122|3300042891__vRhyme_523|3300034021__vRhyme_113|3300034061__vRhyme_59|3300042312__vRhyme_4|3300042474__vRhyme_200"} = "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Methylococcales;f__Methylomonadaceae;g__Methylomonas;s__Methylomonas sp903931355";
$Viral_host_pair{"3300042906__vRhyme_379|3300042327__vRhyme_unbinned1|3300042373__vRhyme_unbinned1|3300043772__vRhyme_626|3300042515__vRhyme_37|3300033993__vRhyme_676|3300042349__vRhyme_122"} = "d__Bacteria;p__Cyanobacteria;c__Cyanobacteriia;o__PCC-6307;f__Cyanobiaceae;g__Cyanobium;s__Cyanobium sp014191755";
$Viral_host_pair{"3300043438__vRhyme_unbinned495|3300042862__vRhyme_118|3300042503__vRhyme_196"} = "d__Bacteria;p__Cyanobacteria;c__Cyanobacteriia;o__Cyanobacteriales;f__Nostocaceae;g__Dolichospermum;s__Dolichospermum sp000312705";
$Viral_host_pair{"3300043467__vRhyme_122"} = "d__Bacteria;p__Cyanobacteria;c__Cyanobacteriia;o__Cyanobacteriales;f__Nostocaceae;g__Dolichospermum;s__Dolichospermum sp001277295";
$Viral_host_pair{"3300044627__vRhyme_unbinned4646"} = "d__Bacteria;p__Cyanobacteria;c__Cyanobacteriia;o__Cyanobacteriales;f__Microcoleaceae;g__Lyngbya;s__Lyngbya robusta";
$Viral_host_pair{"3300047095__vRhyme_323"} = "d__Bacteria;p__Cyanobacteria;c__Cyanobacteriia;o__Cyanobacteriales;f__Microcystaceae;g__Microcystis;s__Microcystis aeruginosa";

my %Host_species = (); 
foreach my $viral_species (sort keys %Viral_host_pair){
	my $host_species = $Viral_host_pair{$viral_species};
	$Host_species{$host_species} = 1;
}

# Step 2 Store the the viral species cov
## Step 2.1 Store Viral_gn2IMG2cov_norm_filtered.txt
my %Viral_gn2IMG2cov_norm_filtered = (); # $viral_gn => $img => $cov_norm
my @Header = (); # Store the header
open IN, "MetaPop/Viral_gn2IMG2cov_norm_filtered.txt";
while (<IN>){
	chomp;
	if (/^Head/){
		my @tmp = split (/\t/);
		@Header = @tmp;
	}else{
		my @tmp = split (/\t/);
		my $viral_gn = $tmp[0];
		for(my $i=1; $i<=$#tmp; $i++){
			my $img = $Header[$i];
			my $cov_norm = $tmp[$i];
			$Viral_gn2IMG2cov_norm_filtered{$viral_gn}{$img} = $cov_norm;
		}
	}
}
close IN;

my @IMG = (); # Store the array of IMG ID
@IMG = @Header; 
shift @IMG; # Delete the first element of the array

## Step 2.2 Store %Viral_species2cov_array hash
my %Viral_species2cov_array = (); # $viral_gn => $cov_array (include the gn name as the first element)
foreach my $viral_species (sort keys %Viral_host_pair){
	my @Cov_array = ();
	push @Cov_array, $viral_species; # The first element should be $viral_species
	
	for(my $i=0; $i<=$#IMG; $i++){
		my $img = $IMG[$i];
		if ($viral_species !~ /\|/){ # If the viral_species is just one species
			push @Cov_array, $Viral_gn2IMG2cov_norm_filtered{$viral_species}{$img};
		}else{ # If the viral_species is made by two or more species
			my @Viral_species = split (/\|/, $viral_species);
			my $cov = 0;
			foreach my $key (@Viral_species){
				if ($Viral_gn2IMG2cov_norm_filtered{$key}{$img} ne "NA"){
					$cov += $Viral_gn2IMG2cov_norm_filtered{$key}{$img};
				}
			}
			
			if ($cov == 0){
				$cov = "NA";
			}
			
			push @Cov_array, $cov;
		}
	}
	
	$Viral_species2cov_array{$viral_species} = join("\t", @Cov_array);
}

# Step 3 Store all MAG abundance (normalized)
## Step 3.1 Store TYMEFLIES_all_MAGs_stat.txt
my %MAG2scfs = (); # $mag => $scfs (collection of $scf, separated by "\,")
my %MAG2tax = (); # $mag => $tax
my %IMG2mag = (); # $img_id => $mags (collection of $mag, separated by "\,")
my %Host_species2mag = (); # $host_species => $mags (collection of $mag, separated by "\,")
open IN, "cat /storage1/data11/TYMEFLIES_phage/Binning_Data/TYMEFLIES_all_MAGs_stat.txt /storage1/data11/TYMEFLIES_phage/Binning_Data/TYMEFLIES_all_MAGs_stat.additional_49_metagenomes.txt |";
while (<IN>){
	chomp;
	if (!/^IMG/){
		my @tmp = split (/\t/);
		my $mag = $tmp[0];
		my $scf = $tmp[4];
		my $tax = $tmp[1];
		$MAG2scfs{$mag} = $scf;
		$MAG2tax{$mag} = $tax;
		
		my ($img_id) = $mag =~ /^(.+?)\_/;
		if (!exists $IMG2mag{$img_id}){
			$IMG2mag{$img_id} = $mag;
		}else{
			$IMG2mag{$img_id} .= "\,".$mag;
		}
		
		if (exists $Host_species{$tax}){
			if (! exists $Host_species2mag{$tax}){
				$Host_species2mag{$tax} = $mag;
			}else{
				$Host_species2mag{$tax} .= "\,".$mag;
			}
		}
	}
}
close IN;

## Step 3.2 Get read numbers for each metagenome
my %IMG_ID2read_num = ();
open IN, "TYMEFLIES_metagenome_info.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		my @tmp = split (/\t/);
		my $img_id = $tmp[0];
		my $read_num = $tmp[13];
		$IMG_ID2read_num{$img_id} = $read_num;
	}
}	
close IN;

## Step 3.3 Make the %MAG2IMG2abun hash
my %MAG2IMG2abun = (); # $mag => $img_id => $abun (normalized)
open IN, "ls /storage1/data11/TYMEFLIES_phage/33*/33*.id97.coverm_depth.txt |";
while (<IN>){
	chomp;
	my $file = $_;
	my ($img_id) = $file =~ /^.+\/(.+?)\.id97\.coverm_depth\.txt/;
	my %Scf2depth = (); # Store all scaffold to depth information
	open INN, "$file";
	while (<INN>){
		chomp;
		if (/^Ga/){
			my @tmp = split (/\t/);
			my $scf = $tmp[0];
			my $depth = $tmp[1];
			$Scf2depth{$scf} = $depth;
		}
	}
	close INN;
	
	if (exists $IMG2mag{$img_id}){
		my @MAGs = split (/\,/, $IMG2mag{$img_id});
		
		foreach my $mag (@MAGs){
			my @Scfs = split (/\,/, $MAG2scfs{$mag});
			my @Depths = ();
			foreach my $scf (@Scfs){
				push @Depths, $Scf2depth{$scf};
			}
			my $depth_mean = _avg(@Depths);
			
			my $read_num = $IMG_ID2read_num{$img_id};
			my $abun_normalized = $depth_mean / ($read_num / 100000000); # Normalized abun (or depth) for the $mag
			$MAG2IMG2abun{$mag}{$img_id} = $abun_normalized;
		}
	}
}
close IN;

## Step 3.4 Store %Host_species2IMG2abun hash
my %Host_species2IMG2abun = (); # $host_species => $img => $abun
foreach my $host_species (sort keys %Host_species){
	foreach my $img (sort keys %IMG_ID2read_num){
		my $abun = 0;
		my @MAGs = split (/\,/, $Host_species2mag{$host_species});
		foreach my $mag (@MAGs){
			if ($MAG2IMG2abun{$mag}{$img}){ # To check if the abundance value is 0 or not
				$abun += $MAG2IMG2abun{$mag}{$img};
			}
		}
		$Host_species2IMG2abun{$host_species}{$img} = $abun;
	}
}

# Step 3.5 Store %Host_species2cov_array hash
my %Host_species2cov_array = (); # $host_species => $cov_array (include the host_species as the first element)
foreach my $host_species (sort keys %Host_species2IMG2abun){
	my @Cov_array = ();
	push @Cov_array, $host_species; # The first element should be $host_species
	
	for(my $i=0; $i<=$#IMG; $i++){
		my $img = $IMG[$i];
		push @Cov_array, $Host_species2IMG2abun{$host_species}{$img};
	}
	
	$Host_species2cov_array{$host_species} = join("\t", @Cov_array);
}

# Step 4 Conduct Spearman's correlation test
my %Corr_result = (); # $viral_species => [0] $viral_host [1] $corr_coeff [2] $pvalue
my $i = 1; # Store the viral species to host species pair
`mkdir MetaPop/tmp_folder_for_correlation_test`;
foreach my $viral_species (sort keys %Viral_species2cov_array){
	open OUT, ">MetaPop/tmp_folder_for_correlation_test/$i.corr_test.input.txt";
	print OUT "$Viral_species2cov_array{$viral_species}\n";
	my $host_species = $Viral_host_pair{$viral_species};
	print OUT "$Host_species2cov_array{$host_species}\n";
	close OUT;
	
	`python3 /storage1/data14/for_chao/calc_spearman_correlation.py -i MetaPop/tmp_folder_for_correlation_test/$i.corr_test.input.txt -o MetaPop/tmp_folder_for_correlation_test/$i.corr_test.output.txt`;

	open IN, "MetaPop/tmp_folder_for_correlation_test/$i.corr_test.output.txt";
	while (<IN>){
		chomp;
		if (!/^var1/){
			my @tmp = split (/\t/);
			my $viral_species = $tmp[0];
			my $host_species = $tmp[1];
			my $pvalue = $tmp[2];
			my $corr_coeff = $tmp[3];
			$Corr_result{$viral_species}[0] = $host_species;
			$Corr_result{$viral_species}[1] = $corr_coeff;
			$Corr_result{$viral_species}[2] = $pvalue;
		}
	}
	close IN;
	$i++;
}
#`rm -r MetaPop/tmp_folder_for_correlation_test`;

## Step 5 Write down the correlation test result 
open OUT, ">MetaPop/Correlation_between_viral_species_and_host_species.txt";
print OUT "Viral species\tHost species\tcorr coeff\tp value\n";
foreach my $viral_species (sort keys %Corr_result){
	print OUT "$viral_species\t$Corr_result{$viral_species}[0]\t$Corr_result{$viral_species}[1]\t$Corr_result{$viral_species}[2]\n";
}
close OUT;



# Subroutine 
sub _avg {
    my $total;
    $total += $_ foreach @_;
    # sum divided by number of components.
    return $total / @_;
}