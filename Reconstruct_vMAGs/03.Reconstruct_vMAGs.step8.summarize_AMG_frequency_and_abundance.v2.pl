#!/usr/bin/perl

use strict;
use warnings;

# AIM: 1) Summarize AMG frequencies for each season
#      2) Calculate AMG abundances for each season
#      3) Calculate AMG abundances for each year_season

# Step 1. Store AMG summary table
my %AMG_summary = (); # $pro => $ko
my %KOs= (); # $ko => 1;
my %IMG2date = (); # $img_id => $date_n_season
open IN, "AMG_analysis/AMG_summary.txt";
while (<IN>){
	chomp;
	if (!/^Pro/){
		my @tmp = split (/\t/);
		my $pro = $tmp[0];
		my $ko = $tmp[2];
		my $ko_detail = $tmp[3];
		my $date_n_season = $tmp[1];
		$AMG_summary{$pro} = $ko;
		my ($img_id) = $pro =~ /^(33.+?)\_/;
		$IMG2date{$img_id} = $date_n_season;
		$KOs{$ko} = $ko_detail;
	}
}
close IN;
 
# Step 2. Summarize AMG frequencies (the metagenome size was normalized by 100M reads) in each metagenome
## Step 2.1 Get read numbers for each metagenome
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

## Step 2.2 Get AMG frequencies 
my %KO2metagenome2freq = (); # $ko => $img_id => $freq
foreach my $ko (sort keys %KOs){
	foreach my $pro (sort keys %AMG_summary){
		if ($AMG_summary{$pro} eq $ko){
			my ($img_id) = $pro =~ /^(33.+?)\_/;
			$KO2metagenome2freq{$ko}{$img_id}++;
		}
	}
}

## Step 2.3 Get AMG frequencies normalized by each metagenome read number as 100M
my %KO2metagenome2freq_normalized = (); # $ko => $img_id => $freq
foreach my $ko (sort keys %KO2metagenome2freq){
	foreach my $img_id (sort keys %IMG_ID2read_num){
		my $freq = $KO2metagenome2freq{$ko}{$img_id};
		my $read_num = $IMG_ID2read_num{$img_id};
		my $freq_normalized = 0;
		if ($freq){
			$freq_normalized = $freq / ($read_num / 100000000);
		}
		$KO2metagenome2freq_normalized{$ko}{$img_id} = $freq_normalized;
	}
}

## Step 2.4 Get AMG frequencies for each season (the number of metagenome each season was normalized)
my %Season2num = (); # $season => $num_metagenome; Store how many samples (metagenomes) are there in each season for the whole datasets
my %Season2img_id = (); # $season => collection of $img_id separeted by "\t"
my %Year_season2num = (); # $year_season => $num_metagenome
my %Year_season2img_id = (); # $year_season (for example, "2010-Ice-on") => collection of $img_id separeted by "\t"
open IN, "TYMEFLIES_metagenome_info.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		my @tmp = split (/\t/);
		my $img_id = $tmp[0];
		my $date = $tmp[8];
		my ($year) = $date =~ /(\d\d\d\d)-\d\d-\d\d/; 
		my $season = $tmp[10];
		$Season2num{$season}++;
		if (!exists $Season2img_id{$season}){
			$Season2img_id{$season} = $img_id;
		}else{
			$Season2img_id{$season} .= "\t".$img_id;
		}
		
		my $year_season = $year.'-'.$season;
		$Year_season2num{$year_season}++;
		if (!exists $Year_season2img_id{$year_season}){
			$Year_season2img_id{$year_season} = $img_id;
		}else{
			$Year_season2img_id{$year_season} .= "\t".$img_id;
		}		
	}
}
close IN;

my %KO2season2freq = (); # $ko => $season => $freq 
foreach my $ko (sort keys %KO2metagenome2freq){
	foreach my $season (sort keys %Season2num){
		my $num_metagenome = $Season2num{$season}; # num of metagenomes in this season
		my $freq_sum_for_season = 0; # Store the sum of KO frequencies from all metagenomes from this season
		my $freq_for_season_normalized = 0; # Store the normlized KO frequency for this season
			
		my @IMG_ID = split (/\t/,$Season2img_id{$season});
		foreach my $img_id (@IMG_ID){
			my $freq = $KO2metagenome2freq_normalized{$ko}{$img_id}; 
			$freq_sum_for_season += $freq;
		}
		
		if ($freq_sum_for_season){
			$freq_for_season_normalized = $freq_sum_for_season / $num_metagenome;
		}
		
		$KO2season2freq{$ko}{$season} = $freq_for_season_normalized; 
	}
}

## Step 2.5 Write down AMG frequencies (normalized by read number per metagenome and metagenome number per season) for each season
my @Season = ("Spring", "Clearwater", "Early Summer", "Late Summer", "Fall", "Ice-on");
open OUT, ">AMG_analysis/KO2season2freq.txt";
my $row=join("\t", @Season);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %KO2season2freq){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (@Season) {       
                if (exists $KO2season2freq{$tmp1}{$tmp2}){
                        push @tmp, $KO2season2freq{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

# Step 3. Summarize AMG abundances (the metagenome size was normalized by 100M reads) in each metagenome
## Step 3.1 Store scaffold coverage 
my %Scf2cov = (); # $scf => $cov (within individual metagenome)
open IN, "ls /storage1/data11/TYMEFLIES_phage/*/*.id97.coverm_depth.txt | ";
while (<IN>){
	chomp;
	my $file = $_;
	open INN, "$file";
	while (<INN>){
		chomp;
		if (!/^contigName/){
			my @tmp = split (/\t/);
			my $scf = $tmp[0];
			my $cov = $tmp[2];
			$Scf2cov{$scf} = $cov;
		}
	}
	close INN;
}
close IN;

## Step 3.2 Store prophage scaffold coverage (containing 'fragment' in scaffold header)
my %Prophage_scf2cov = (); # $scf => $cov (within individual metagenome)
open IN, "ls /storage1/data11/TYMEFLIES_phage/*/PropagAtE_result/PropagAtE_result.tsv | ";
while (<IN>){
	chomp;
	my $file = $_;
	open INN, "$file";
	while (<INN>){
		chomp;
		if (!/^prophage/){
			my @tmp = split (/\t/);
			my $prophage_scf = $tmp[0];
			my $cov = $tmp[7];
			if ($cov eq 'NA'){
				$cov = 0;
			}			
			$Prophage_scf2cov{$prophage_scf} = $cov;
		}
	}
	close INN;
}
close IN;

## Step 3.3 Get AMG abundances (normalzied by 100M reads per metagenome)
my %KO2metagenome2abun_normalized = (); # $ko => $img_id => $abun_normalized_sum (normalized abundance sum of all AMG)
foreach my $ko (sort keys %KOs){
	foreach my $pro (sort keys %AMG_summary){
		if ($AMG_summary{$pro} eq $ko){
			my $abun_normalized = 0;
			my ($img_id,$scf) = $pro =~ /^(33\d+?)\_\_vRhyme\_.+?\_\_(Ga.+)\_\d+?$/;
			my $abun = 0;
			if ($scf =~ /fragment/){ # If this scaffold is from prophage 
				$abun = $Prophage_scf2cov{$scf}; 
				print "$scf\n"; # Print out the $scf to make a check
			}else{
				$abun = $Scf2cov{$scf}; # If this scaffold is not from prophage 
			}
			my $read_num = $IMG_ID2read_num{$img_id}; 
			
			if ($abun){
				$abun_normalized = $abun / ($read_num / 100000000);
			}
			
			$KO2metagenome2abun_normalized{$ko}{$img_id} += $abun_normalized;

		}
	}
}

## Step 3.4 Get AMG abundances for each season (the number of metagenome each season was normalized)
my %KO2season2abun = (); # $ko => $season => $abun (normalized abundance sum of all AMG) 
foreach my $ko (sort keys %KO2metagenome2abun_normalized){
	foreach my $season (sort keys %Season2num){
		my $num_metagenome = $Season2num{$season}; # num of metagenomes in this season
		my $abun_sum_for_season = 0; # Store the sum of KO abun from all metagenomes from this season
		my $abun_for_season_normalized = 0; # Store the normlized KO abun for this season
			
		my @IMG_ID = split (/\t/,$Season2img_id{$season});
		foreach my $img_id (@IMG_ID){
			my $abun = $KO2metagenome2abun_normalized{$ko}{$img_id}; 
			if ($abun){
				$abun_sum_for_season += $abun;
			}
		}
		
		if ($abun_sum_for_season){
			$abun_for_season_normalized = $abun_sum_for_season / $num_metagenome;
		}
		
		$KO2season2abun{$ko}{$season} = $abun_for_season_normalized; 
	}
}

## Step 3.5 Write down AMG abundances (normalized by read number per metagenome and metagenome number per season) for each season
open OUT, ">AMG_analysis/KO2season2abun.txt";
my $row2=join("\t", @Season);
print OUT "Head\t$row2\n";
foreach my $tmp1 (sort keys %KO2season2abun){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (@Season) {       
                if (exists $KO2season2abun{$tmp1}{$tmp2}){
                        push @tmp, $KO2season2abun{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

## Step 3.6 Get AMG abundances for each year_season (the number of metagenome each year_season was normalized)
my %KO2year_season2abun = (); # $ko => $year_season => $abun (normalized abundance sum of all AMG) 
foreach my $ko (sort keys %KO2metagenome2abun_normalized){
	foreach my $year_season (sort keys %Year_season2num){
		my $num_metagenome = $Year_season2num{$year_season}; # num of metagenomes in this year_season
		my $abun_sum_for_year_season = 0; # Store the sum of KO abun from all metagenomes from this year_season
		my $abun_for_year_season_normalized = 0; # Store the normlized KO abun for this year_season
			
		my @IMG_ID = split (/\t/,$Year_season2img_id{$year_season});
		foreach my $img_id (@IMG_ID){
			my $abun = $KO2metagenome2abun_normalized{$ko}{$img_id}; 
			if ($abun){
				$abun_sum_for_year_season += $abun;
			}
		}
		
		if ($abun_sum_for_year_season){
			$abun_for_year_season_normalized = $abun_sum_for_year_season / $num_metagenome;
		}
		
		$KO2year_season2abun{$ko}{$year_season} = $abun_for_year_season_normalized; 
	}
}

## Step 3.7 Write down AMG abundances (normalized by read number per metagenome and metagenome number per year_season) for each year_season
my @Year = ('2000','2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019');
my @Year_season = ();
foreach my $year (@Year){
	foreach my $season (@Season){
		my $year_season = "$year\-$season";
		push @Year_season, $year_season;
	}
}

open OUT, ">AMG_analysis/KO2year_season2abun.txt";
my $row3=join("\t", @Year_season);
print OUT "Head\t$row3\n";
foreach my $tmp1 (sort keys %KO2year_season2abun){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (@Year_season) {       
                if (exists $KO2year_season2abun{$tmp1}{$tmp2}){
                        push @tmp, $KO2year_season2abun{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

# Step 4. Write down the KO details
open OUT, ">AMG_analysis/KO2season_ko_details.txt";
foreach my $ko (sort keys %KO2season2abun){
	my $ko_details = $KOs{$ko};
	$ko_details =~ s/^'(.*)'$/$1/;
	$ko_details =~ s/'/(prime)/g
	print OUT "$ko\t$ko_details\n";
}
close OUT;

# Step 5. Write down the number of metagenomes in each season and year_season
open OUT, ">AMG_analysis/Season2num_of_metagenomes.txt";
foreach my $season (@Season){
	print OUT "$season\t$Season2num{$season}\n";
}
close OUT;

open OUT, ">AMG_analysis/Year_season2num_of_metagenomes.txt";
foreach my $year_season (@Year_season){
	my $num = 0;
	if (exists $Year_season2num{$year_season}){
		print OUT "$year_season\t$Year_season2num{$year_season}\n";
	}else{
		print OUT "$year_season\t$num\n";
	}
}
close OUT;
