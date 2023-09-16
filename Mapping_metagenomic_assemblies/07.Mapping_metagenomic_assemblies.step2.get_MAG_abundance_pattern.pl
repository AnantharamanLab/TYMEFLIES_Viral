#!/usr/bin/perl

use strict;
use warnings;

# Aim: Get MAG abundance pattern, including seasonly pattern (aggregated across 20 years), year-season pattern, and yearly pattern (for each year)

# Step 1 Store all MAG abundance (normalized)
## Step 1.1 Store Robin_MAG_stat.txt
my %MAG2scfs = (); # $mag => $scfs (collection of $scf, separated by "\,")
my %MAG2tax = (); # $mag => $tax
my %IMG2mag = (); # $img_id => $mags (collection of $mag, separated by "\,")
open IN, "/storage1/data11/TYMEFLIES_phage/Robin_MAGs/Robin_MAG_stat.txt";
while (<IN>){
	chomp;
	if (!/^tymeflies/){
		my @tmp = split (/\t/);
		my $mag = $tmp[5];		
		my $num_in_cluster = $tmp[15];		
		if ($num_in_cluster ne "NA"){
			my ($img_id) = $mag =~ /_(33\d+?)_/;
			my @Scfs = ();
			my $MAG_addr = "/storage1/data11/TYMEFLIES_phage/Robin_MAGs/".$img_id."/".$mag.".fasta";
			my %MAG_seq = _store_seq("$MAG_addr");
			foreach my $header (sort keys %MAG_seq){
				my ($scf) = $header =~ /^>(.+?)$/;
				push @Scfs, $scf;
			}
			$MAG2scfs{$mag} = join("\,", @Scfs);
			my $tax = join(";", @tmp[16..22]);
			$MAG2tax{$mag} = $tax;
			
			if (!exists $IMG2mag{$img_id}){
				$IMG2mag{$img_id} = $mag;
			}else{
				$IMG2mag{$img_id} .= "\,".$mag;
			}	
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

`mkdir MAG_abundance`;
## Step 1.4 Write down the %MAG2IMG2abun result
open OUT, ">MAG_abundance/MAG2IMG2abun.txt";
my $row0=join("\t", sort keys %IMG_ID2read_num);
print OUT "Head\t$row0\n";
foreach my $tmp1 (sort keys %MAG2IMG2abun){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %IMG_ID2read_num) {       
                if (exists $MAG2IMG2abun{$tmp1}{$tmp2}){
                        push @tmp, $MAG2IMG2abun{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

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
		my ($img_id) = $mag =~ /_(33\d+?)_/;
		if (exists $MAG2IMG2abun{$mag}{$img_id}){
			$Family2IMG2abun{$family}{$img_id} += $MAG2IMG2abun{$mag}{$img_id};
		}
	}
}

# Step 3 Store and write down %Family2season2abun hash
## Step 3.1 Store season/year_season to img id info
my %Season2num = (); # $season => $num_metagenome; Store how many samples (metagenomes) are there in each season for the whole datasets
my %Season2img_id = (); # $season => collection of $img_id separeted by "\t"
my %Year_season2num = (); # $year_season => $num_metagenome
my %Year_season2img_id = (); # $year_season (for example, "2010-Ice-on") => collection of $img_id separeted by "\t"
my %Year2num = (); # $year => $num_metagenome
my %Year2img_id = (); # $year (for example, "2010") => collection of $img_id separeted by "\t"
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
		
		$Year2num{$year}++;
		if (!exists $Year2img_id{$year}){
			$Year2img_id{$year} = $img_id;
		}else{
			$Year2img_id{$year} .= "\t".$img_id;
		}		
	}
}
close IN;

## Step 3.2 Get family abundances for each season (the number of metagenome each season was normalized)
my %Family2season2abun = (); # $family => $season => $abun (normalized abundance) 
foreach my $family (sort keys %Family2IMG2abun){
	foreach my $season (sort keys %Season2num){
		my $num_metagenome = $Season2num{$season}; # num of metagenomes in this season
		my $abun_sum_for_season = 0; # Store the sum of family abun from all metagenomes from this season
		my $abun_for_season_normalized = 0; # Store the normlized family abun for this season
			
		my @IMG_ID = split (/\t/,$Season2img_id{$season});
		foreach my $img_id (@IMG_ID){
			my $abun = $Family2IMG2abun{$family}{$img_id}; 
			if ($abun){
				$abun_sum_for_season += $abun;
			}
		}
		
		if ($abun_sum_for_season){
			$abun_for_season_normalized = $abun_sum_for_season / $num_metagenome;
		}
		
		$Family2season2abun{$family}{$season} = $abun_for_season_normalized; 
	}
}

## Step 3.3 Write down Family abundances (normalized by read number per metagenome and metagenome number per season) for each season
my @Season = ("Spring", "Clearwater", "Early Summer", "Late Summer", "Fall", "Ice-on");
open OUT, ">MAG_abundance/Family2season2abun.txt";
my $row=join("\t", @Season);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %Family2season2abun){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (@Season) {       
                if (exists $Family2season2abun{$tmp1}{$tmp2}){
                        push @tmp, $Family2season2abun{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

# Step 4 Store and write down %Family2year_season2abun hash
## Step 4.1 Get family abundances for each season (the number of metagenome each season was normalized)
my %Family2year_season2abun = (); # $family => $year_season => $abun (normalized abundance sum of all AMG within the family) 
foreach my $family (sort keys %Family2IMG2abun){
	foreach my $year_season (sort keys %Year_season2num){
		my $num_metagenome = $Year_season2num{$year_season}; # num of metagenomes in this year_season
		my $abun_sum_for_year_season = 0; # Store the sum of family abun from all metagenomes from this year_season
		my $abun_for_year_season_normalized = 0; # Store the normlized family abun for this year_season
			
		my @IMG_ID = split (/\t/,$Year_season2img_id{$year_season});
		foreach my $img_id (@IMG_ID){
			my $abun = $Family2IMG2abun{$family}{$img_id}; 
			if ($abun){
				$abun_sum_for_year_season += $abun;
			}
		}
		
		if ($abun_sum_for_year_season){
			$abun_for_year_season_normalized = $abun_sum_for_year_season / $num_metagenome;
		}
		
		$Family2year_season2abun{$family}{$year_season} = $abun_for_year_season_normalized; 
	}
}

## Step 4.2 Write down family abundances (normalized by read number per metagenome and metagenome number per year_season) for each year_season
my @Year = ('2000','2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019');
my @Year_season = ();
foreach my $year (@Year){
	foreach my $season (@Season){
		my $year_season = "$year\-$season";
		push @Year_season, $year_season;
	}
}

open OUT, ">MAG_abundance/Family2year_season2abun.txt";
my $row2=join("\t", @Year_season);
print OUT "Head\t$row2\n";
foreach my $tmp1 (sort keys %Family2year_season2abun){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (@Year_season) {       
                if (exists $Family2year_season2abun{$tmp1}{$tmp2}){
                        push @tmp, $Family2year_season2abun{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

# Step 5 Store and write down %Family2year2abun hash
## Step 5.1 Get family abundances for each year (the number of metagenome each year was normalized)
my %Family2year2abun = (); # $family => $year => $abun (normalized abundance) 
foreach my $family (sort keys %Family2IMG2abun){
	foreach my $year (sort keys %Year2num){
		my $num_metagenome = $Year2num{$year}; # num of metagenomes in this year
		my $abun_sum_for_year = 0; # Store the sum of family abun from all metagenomes from this year
		my $abun_for_year_normalized = 0; # Store the normlized family abun for this year
			
		my @IMG_ID = split (/\t/,$Year2img_id{$year});
		foreach my $img_id (@IMG_ID){
			my $abun = $Family2IMG2abun{$family}{$img_id}; 
			if ($abun){
				$abun_sum_for_year += $abun;
			}
		}
		
		if ($abun_sum_for_year){
			$abun_for_year_normalized = $abun_sum_for_year / $num_metagenome;
		}
		
		$Family2year2abun{$family}{$year} = $abun_for_year_normalized; 
	}
}

## Step 5.2 Write down Family abundances (normalized by read number per metagenome and metagenome number per year) for each year
open OUT, ">MAG_abundance/Family2year2abun.txt";
my $row3=join("\t", sort keys %Year2num);
print OUT "Head\t$row3\n";
foreach my $tmp1 (sort keys %Family2year2abun){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %Year2num) {       
                if (exists $Family2year2abun{$tmp1}{$tmp2}){
                        push @tmp, $Family2year2abun{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

# Step 6. Write down the number of metagenomes in each season, year_season, and year
open OUT, ">MAG_abundance/Season2num_of_metagenomes.txt";
foreach my $season (@Season){
	print OUT "$season\t$Season2num{$season}\n";
}
close OUT;

open OUT, ">MAG_abundance/Year_season2num_of_metagenomes.txt";
foreach my $year_season (@Year_season){
	my $num = 0;
	if (exists $Year_season2num{$year_season}){
		print OUT "$year_season\t$Year_season2num{$year_season}\n";
	}else{
		print OUT "$year_season\t$num\n";
	}
}
close OUT;

open OUT, ">MAG_abundance/Year2num_of_metagenomes.txt";
foreach my $year (@Year){
	my $num = 0;
	if (exists $Year2num{$year}){
		print OUT "$year\t$Year2num{$year}\n";
	}else{
		print OUT "$year\t$num\n";
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
