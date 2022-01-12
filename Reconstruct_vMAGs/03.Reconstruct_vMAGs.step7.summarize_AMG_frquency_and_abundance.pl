#!/usr/bin/perl

use strict;
use warnings;

# AIM: 1) Summarize AMG frequencies monthly
#      2) Calculate AMG abundances monthly
#      3) Calculate AMG abundances for each year_month

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

## Step 2.4 Get AMG frequencies for each month (the number of metagenome each month was normalized)
my %Month2num = (); # $month => $num_metagenome; Store how many samples (metagenomes) are there in each month for the whole datasets
my %Month2img_id = (); # $month => collection of $img_id separeted by "\t"
my %Year_month2num = (); # $year_month => $num_metagenome
my %Year_month2img_id = (); # $year_month (for example, "2010-06") => collection of $img_id separeted by "\t"
open IN, "TYMEFLIES_metagenome_info.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		my @tmp = split (/\t/);
		my $img_id = $tmp[0];
		my $date = $tmp[8];
		my ($month) = $date =~ /\d\d\d\d-(\d\d)-\d\d/; 
		$Month2num{$month}++;
		if (!exists $Month2img_id{$month}){
			$Month2img_id{$month} = $img_id;
		}else{
			$Month2img_id{$month} .= "\t".$img_id;
		}
		
		my ($year_month) = $date =~ /(\d\d\d\d-\d\d)-\d\d/;
		$Year_month2num{$year_month}++;
		if (!exists $Year_month2img_id{$year_month}){
			$Year_month2img_id{$year_month} = $img_id;
		}else{
			$Year_month2img_id{$year_month} .= "\t".$img_id;
		}		
	}
}
close IN;

my %KO2month2freq = (); # $ko => $month => $freq 
foreach my $ko (sort keys %KO2metagenome2freq){
	foreach my $month (sort keys %Month2num){
		my $num_metagenome = $Month2num{$month}; # num of metagenomes in this month
		my $freq_sum_for_month = 0; # Store the sum of KO frequencies from all metagenomes from this month
		my $freq_for_month_normalized = 0; # Store the normlized KO frequency for this month
			
		my @IMG_ID = split (/\t/,$Month2img_id{$month});
		foreach my $img_id (@IMG_ID){
			my $freq = $KO2metagenome2freq_normalized{$ko}{$img_id}; 
			$freq_sum_for_month += $freq;
		}
		
		if ($freq_sum_for_month){
			$freq_for_month_normalized = $freq_sum_for_month / $num_metagenome;
		}
		
		$KO2month2freq{$ko}{$month} = $freq_for_month_normalized; 
	}
}

## Step 2.5 Write down AMG frequencies (normalized by read number per metagenome and metagenome number per month) for each month
open OUT, ">AMG_analysis/KO2month2freq.txt";
my $row=join("\t", sort keys %Month2num);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %KO2month2freq){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %Month2num) {       
                if (exists $KO2month2freq{$tmp1}{$tmp2}){
                        push @tmp, $KO2month2freq{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

# Step 3. Summarize AMG abundances (the metagenome size was normalized by 100M reads) in each metagenome
## Step 3.1 Store scaffold coverage for each phage genome
my %Scf2cov = (); # $scf => $cov (within individual metagenome)
open IN, "ls 33**/vRhyme_result/vRhyme_coverage_files/vRhyme_coverage_values.tsv | ";
while (<IN>){
	chomp;
	my $file = $_;
	my ($img_id) = $file =~ /(33\d+?)\//;
	open INN, "$file";
	while (<INN>){
		chomp;
		if (!/^scaffold/){
			my @tmp = split (/\t/);
			my $scf = $tmp[0];
			my $cov = $tmp[1];
			$Scf2cov{$scf} = $cov;
		}
	}
	close INN;
}
close IN;

## Step 3.2 Get AMG abundances (normalzied by 100M reads per metagenome)
my %KO2metagenome2abun_normalized = (); # $ko => $img_id => $abun_normalized_sum (normalized abundance sum of all AMG)
foreach my $ko (sort keys %KOs){
	foreach my $pro (sort keys %AMG_summary){
		if ($AMG_summary{$pro} eq $ko){
			my $abun_normalized = 0;
			my ($img_id,$scf) = $pro =~ /^(33\d+?)\_\_vRhyme\_.+?\_\_(Ga.+)\_\d+?$/;
			my $read_num = $IMG_ID2read_num{$img_id}; 
			my $abun = $Scf2cov{$scf}; 
			if ($abun){
				$abun_normalized = $abun / ($read_num / 100000000);
			}
			
			$KO2metagenome2abun_normalized{$ko}{$img_id} += $abun_normalized;

		}
	}
}

## Step 3.3 Get AMG abundances for each month (the number of metagenome each month was normalized)
my %KO2month2abun = (); # $ko => $month => $abun (normalized abundance sum of all AMG) 
foreach my $ko (sort keys %KO2metagenome2abun_normalized){
	foreach my $month (sort keys %Month2num){
		my $num_metagenome = $Month2num{$month}; # num of metagenomes in this month
		my $abun_sum_for_month = 0; # Store the sum of KO abun from all metagenomes from this month
		my $abun_for_month_normalized = 0; # Store the normlized KO abun for this month
			
		my @IMG_ID = split (/\t/,$Month2img_id{$month});
		foreach my $img_id (@IMG_ID){
			my $abun = $KO2metagenome2abun_normalized{$ko}{$img_id}; 
			if ($abun){
				$abun_sum_for_month += $abun;
			}
		}
		
		if ($abun_sum_for_month){
			$abun_for_month_normalized = $abun_sum_for_month / $num_metagenome;
		}
		
		$KO2month2abun{$ko}{$month} = $abun_for_month_normalized; 
	}
}

## Step 3.4 Write down AMG abundances (normalized by read number per metagenome and metagenome number per month) for each month
open OUT, ">AMG_analysis/KO2month2abun.txt";
my $row2=join("\t", sort keys %Month2num);
print OUT "Head\t$row2\n";
foreach my $tmp1 (sort keys %KO2month2abun){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %Month2num) {       
                if (exists $KO2month2abun{$tmp1}{$tmp2}){
                        push @tmp, $KO2month2abun{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

## Step 3.5 Get AMG abundances for each year_month (the number of metagenome each year_month was normalized)
my %KO2year_month2abun = (); # $ko => $year_month => $abun (normalized abundance sum of all AMG) 
foreach my $ko (sort keys %KO2metagenome2abun_normalized){
	foreach my $year_month (sort keys %Year_month2num){
		my $num_metagenome = $Year_month2num{$year_month}; # num of metagenomes in this year_month
		my $abun_sum_for_year_month = 0; # Store the sum of KO abun from all metagenomes from this year_month
		my $abun_for_year_month_normalized = 0; # Store the normlized KO abun for this year_month
			
		my @IMG_ID = split (/\t/,$Year_month2img_id{$year_month});
		foreach my $img_id (@IMG_ID){
			my $abun = $KO2metagenome2abun_normalized{$ko}{$img_id}; 
			if ($abun){
				$abun_sum_for_year_month += $abun;
			}
		}
		
		if ($abun_sum_for_year_month){
			$abun_for_year_month_normalized = $abun_sum_for_year_month / $num_metagenome;
		}
		
		$KO2year_month2abun{$ko}{$year_month} = $abun_for_year_month_normalized; 
	}
}

## Step 3.6 Write down AMG abundances (normalized by read number per metagenome and metagenome number per year_month) for each year_month
my @Month = ('01','02','03','04','05','06','07','08','09','10','11','12');
my @Year = ('2000','2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019');
my @Year_month = ();
foreach my $year (@Year){
	foreach my $month (@Month){
		my $year_month = "$year\-$month";
		push @Year_month, $year_month;
	}
}

open OUT, ">AMG_analysis/KO2year_month2abun.txt";
my $row3=join("\t", @Year_month);
print OUT "Head\t$row3\n";
foreach my $tmp1 (sort keys %KO2year_month2abun){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (@Year_month) {       
                if (exists $KO2year_month2abun{$tmp1}{$tmp2}){
                        push @tmp, $KO2year_month2abun{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

# Step 4. Write down the KO details
open OUT, ">AMG_analysis/KO2month_ko_details.txt";
foreach my $ko (sort keys %KO2month2abun){
	my $ko_details = $KOs{$ko};
	$ko_details =~ s/\'//g;
	print OUT "$ko\t$ko_details\n";
}
close OUT;

# Step 5. Write down the number of metagenomes in each month and year_month
open OUT, ">AMG_analysis/Month2num_of_metagenomes.txt";
foreach my $month (@Month){
	print OUT "$month\t$Month2num{$month}\n";
}
close OUT;

open OUT, ">AMG_analysis/Year_month2num_of_metagenomes.txt";
foreach my $year_month (@Year_month){
	my $num = 0;
	if (exists $Year_month2num{$year_month}){
		print OUT "$year_month\t$Year_month2num{$year_month}\n";
	}else{
		print OUT "$year_month\t$num\n";
	}
}
close OUT;














