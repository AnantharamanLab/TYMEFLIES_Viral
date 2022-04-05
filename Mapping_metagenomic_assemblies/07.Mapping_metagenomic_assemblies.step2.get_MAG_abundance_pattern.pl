#!/usr/bin/perl

use strict;
use warnings;

# Aim: Get MAG abundance pattern, both monthly pattern (aggregated across 20 years) and yearly pattern (for each year)

# Step 1 Store all MAG abundance (normalized)
## Step 1.1 Store TYMEFLIES_all_MAGs_stat.txt
my %MAG2scfs = (); # $mag => $scfs (collection of $scf, separated by "\,")
my %MAG2tax = (); # $mag => $tax
my %IMG2mag = (); # $img_id => $mags (collection of $mag, separated by "\,")
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
	}
}
close IN;

## Step 1.2 Get read numbers for each metagenome
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

## Step 1.3 Make the %MAG2IMG2abun hash
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

# Step 2 Get family (contains the order too) to $img to $abun
my %Family2mag = (); # $family (contains the order in the front too) => $mags (collection of $mag, separated by "\,")
foreach my $mag (sort keys %MAG2tax){
	my $tax = $MAG2tax{$mag};
	my @tmp = split (/\;/, $tax);
	
	my $family = $tmp[3]."\;".$tmp[4];
	
	if (!exists $Family2mag{$family}){
		$Family2mag{$family} = $mag;
	}else{
		$Family2mag{$family} .= "\,".$mag;
	}
}

my %Family2IMG2abun = (); # $family => $img_id => $abun (normalized)
foreach my $family (sort keys %Family2mag){
	my @MAGs = split (/\,/, $Family2mag{$family});
	
	foreach my $mag (@MAGs){
		my ($img_id) = $mag =~ /^(.+?)\_/;
		if (exists $MAG2IMG2abun{$mag}{$img_id}){
			$Family2IMG2abun{$family}{$img_id} += $MAG2IMG2abun{$mag}{$img_id};
		}
	}
}

# Step 3 Store and write down %Family2month2abun hash
## Step 3.1 Store month/year_month to img id info
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

## Step 3.2 Get family abundances for each month (the number of metagenome each month was normalized)
my %Family2month2abun = (); # $family => $month => $abun (normalized abundance) 
foreach my $family (sort keys %Family2IMG2abun){
	foreach my $month (sort keys %Month2num){
		my $num_metagenome = $Month2num{$month}; # num of metagenomes in this month
		my $abun_sum_for_month = 0; # Store the sum of family abun from all metagenomes from this month
		my $abun_for_month_normalized = 0; # Store the normlized family abun for this month
			
		my @IMG_ID = split (/\t/,$Month2img_id{$month});
		foreach my $img_id (@IMG_ID){
			my $abun = $Family2IMG2abun{$family}{$img_id}; 
			if ($abun){
				$abun_sum_for_month += $abun;
			}
		}
		
		if ($abun_sum_for_month){
			$abun_for_month_normalized = $abun_sum_for_month / $num_metagenome;
		}
		
		$Family2month2abun{$family}{$month} = $abun_for_month_normalized; 
	}
}

## Step 3.3 Write down Family abundances (normalized by read number per metagenome and metagenome number per month) for each month
`mkdir MAG_abundance`;
open OUT, ">MAG_abundance/Family2month2abun.txt";
my $row=join("\t", sort keys %Month2num);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %Family2month2abun){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %Month2num) {       
                if (exists $Family2month2abun{$tmp1}{$tmp2}){
                        push @tmp, $Family2month2abun{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

# Step 4 Store and write down %Family2year_month2abun hash
## Step 4.1 Get family abundances for each month (the number of metagenome each month was normalized)
my %Family2year_month2abun = (); # $family => $year_month => $abun (normalized abundance sum of all AMG within the family) 
foreach my $family (sort keys %Family2IMG2abun){
	foreach my $year_month (sort keys %Year_month2num){
		my $num_metagenome = $Year_month2num{$year_month}; # num of metagenomes in this year_month
		my $abun_sum_for_year_month = 0; # Store the sum of family abun from all metagenomes from this year_month
		my $abun_for_year_month_normalized = 0; # Store the normlized family abun for this year_month
			
		my @IMG_ID = split (/\t/,$Year_month2img_id{$year_month});
		foreach my $img_id (@IMG_ID){
			my $abun = $Family2IMG2abun{$family}{$img_id}; 
			if ($abun){
				$abun_sum_for_year_month += $abun;
			}
		}
		
		if ($abun_sum_for_year_month){
			$abun_for_year_month_normalized = $abun_sum_for_year_month / $num_metagenome;
		}
		
		$Family2year_month2abun{$family}{$year_month} = $abun_for_year_month_normalized; 
	}
}

## Step 4.2 Write down family abundances (normalized by read number per metagenome and metagenome number per year_month) for each year_month
my @Month = ('01','02','03','04','05','06','07','08','09','10','11','12');
my @Year = ('2000','2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019');
my @Year_month = ();
foreach my $year (@Year){
	foreach my $month (@Month){
		my $year_month = "$year\-$month";
		push @Year_month, $year_month;
	}
}

open OUT, ">MAG_abundance/Family2year_month2abun.txt";
my $row2=join("\t", @Year_month);
print OUT "Head\t$row2\n";
foreach my $tmp1 (sort keys %Family2year_month2abun){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (@Year_month) {       
                if (exists $Family2year_month2abun{$tmp1}{$tmp2}){
                        push @tmp, $Family2year_month2abun{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

# Step 5. Write down the number of metagenomes in each month and year_month
open OUT, ">MAG_abundance/Month2num_of_metagenomes.txt";
foreach my $month (@Month){
	print OUT "$month\t$Month2num{$month}\n";
}
close OUT;

open OUT, ">MAG_abundance/Year_month2num_of_metagenomes.txt";
foreach my $year_month (@Year_month){
	my $num = 0;
	if (exists $Year_month2num{$year_month}){
		print OUT "$year_month\t$Year_month2num{$year_month}\n";
	}else{
		print OUT "$year_month\t$num\n";
	}
}
close OUT;



# Subroutine 
sub _avg {
    my $total;
    $total += $_ foreach @_;
    # sum divided by number of components.
    return $total / @_;
}
