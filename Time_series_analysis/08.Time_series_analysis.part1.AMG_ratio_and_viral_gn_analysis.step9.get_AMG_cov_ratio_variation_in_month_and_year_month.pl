#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw(sum);

# Aim: Parse to get the AMG variation in month and year_month

# Step 1 Store AMG_gene_cov_ratio_variation_table
my %AMG_gene2IMG2cov_ratio = (); # $amg_gene => $img => $cov_ratio
my %AMG_gene2ko = (); # $amg_gene => $ko
my %AMG_gene2distribution = (); # $amg_gene => $distribution
my @Header = (); # Store the header
open IN, "MetaPop/AMG_gene_cov_ratio_variation_table.txt";
while (<IN>){
	chomp;
	if (/^Head/){
		@Header = split (/\t/);
	}else{
		my @tmp = split (/\t/);
		my $amg_gene = $tmp[0];
		my $ko = $tmp[1];
		my $distribution = $tmp[4];
		$AMG_gene2ko{$amg_gene} = $ko;
		$AMG_gene2distribution{$amg_gene} = $distribution;
		for(my $i=7; $i<=$#tmp; $i++){
			my $img = $Header[$i];
			my $cov_ratio = $tmp[$i];
			$AMG_gene2IMG2cov_ratio{$amg_gene}{$img} = $cov_ratio
		}
	}
}
close IN;

# Step 2 Store AMG KO information and month and year_month metagenome information
## Step 2.1 Store AMG KO information
my %AMG_summary = (); # $pro => $ko
my %KOs= (); # $ko => 1;
my %IMG2date = (); # $img_id => $date_n_season
my %AMG_containing_viral_gn = (); # $gn => 1
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
		
		my ($gn) = $pro =~ /^(.+?\_\_.+?)\_\_/;
		$AMG_containing_viral_gn{$gn} = 1;
	}
}
close IN;

## Step 2.2 Store month metagenome information
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

my @Month = ('01','02','03','04','05','06','07','08','09','10','11','12');
my @Year = ('2000','2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019');
my @Year_month = ();
foreach my $year (@Year){
	foreach my $month (@Month){
		my $year_month = "$year\-$month";
		push @Year_month, $year_month;
	}
}

# Step 3 Make %AMG_gene2month2cov_ratio hash
my %AMG_gene2month2cov_ratio = (); # $amg_gene => $month => $cov_ratio
                                   # Only use $amg_gene with distribution >= 5
foreach my $amg_gene (sort keys %AMG_gene2distribution){
	my $distribution = $AMG_gene2distribution{$amg_gene};
	if ($distribution >= 5){
		foreach my $month (sort keys %Month2img_id){
			my @IMG_ID = split (/\t/,$Month2img_id{$month}); # Store all metagenomes from this month
			my @Cov_ratio_collection = (); # Store all the non-"NA" cov ratio values 
			my $cov_ratio_for_this_month = 0; # Store the cov ratio for this month (the mean value of all non-"NA" values)
			foreach my $img (@IMG_ID){
				my $cov_ratio = $AMG_gene2IMG2cov_ratio{$amg_gene}{$img};
				if ($cov_ratio ne "NA"){
					push @Cov_ratio_collection, $cov_ratio;
				}
			}
			
			if (@Cov_ratio_collection){
				$cov_ratio_for_this_month = mean(@Cov_ratio_collection);
			}
			$AMG_gene2month2cov_ratio{$amg_gene}{$month} = $cov_ratio_for_this_month;
		}
	}
}

## Write down %AMG_gene2month2cov_ratio hash
open OUT, ">MetaPop/AMG_gene2month2cov_ratio.txt";
my $row=join("\t", sort keys %Month2num);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %AMG_gene2month2cov_ratio){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %Month2num) {       
                if (exists $AMG_gene2month2cov_ratio{$tmp1}{$tmp2}){
                        push @tmp, $AMG_gene2month2cov_ratio{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

# Step 4 Make %AMG_gene2year_month2cov_ratio hash and write it down
my %AMG_gene2year_month2cov_ratio = (); # $amg_gene => $year_month => $cov_ratio
                                   # Only use $amg_gene with distribution >= 5
foreach my $amg_gene (sort keys %AMG_gene2distribution){
	my $distribution = $AMG_gene2distribution{$amg_gene};
	if ($distribution >= 5){
		foreach my $year_month (sort keys %Year_month2img_id){
			my @IMG_ID = split (/\t/,$Year_month2img_id{$year_month}); # Store all metagenomes from this year_month
			my @Cov_ratio_collection = (); # Store all the non-"NA" cov ratio values 
			my $cov_ratio_for_this_year_month = 0; # Store the cov ratio for this year_month (the mean value of all non-"NA" values)
			foreach my $img (@IMG_ID){
				my $cov_ratio = $AMG_gene2IMG2cov_ratio{$amg_gene}{$img};
				if ($cov_ratio ne "NA"){
					push @Cov_ratio_collection, $cov_ratio;
				}
			}
			
			if (@Cov_ratio_collection){
				$cov_ratio_for_this_year_month = mean(@Cov_ratio_collection);
			}
			$AMG_gene2year_month2cov_ratio{$amg_gene}{$year_month} = $cov_ratio_for_this_year_month;
		}
	}
}

## Write down %AMG_gene2year_month2cov_ratio hash
open OUT, ">MetaPop/AMG_gene2year_month2cov_ratio.txt";
my $row2=join("\t", @Year_month);
print OUT "Head\t$row2\n";
foreach my $tmp1 (sort keys %AMG_gene2year_month2cov_ratio){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (@Year_month) {       
                if (exists $AMG_gene2year_month2cov_ratio{$tmp1}{$tmp2}){
                        push @tmp, $AMG_gene2year_month2cov_ratio{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

# Step 5 Make %AMG_gene_containing_viral_gn2month2cov_ratio hash and write it down
## Step 5.1 Store Viral_gn2IMG2cov_norm_filtered.txt
my %Viral_gn2IMG2cov_norm_filtered = (); # $viral_gn => $img => $cov_norm
my @Header2 = (); # Store the header
open IN, "MetaPop/Viral_gn2IMG2cov_norm_filtered.txt";
while (<IN>){
	chomp;
	if (/^Head/){
		my @tmp = split (/\t/);
		@Header2 = @tmp;
	}else{
		my @tmp = split (/\t/);
		my $viral_gn = $tmp[0];
		for(my $i=1; $i<=$#tmp; $i++){
			my $img = $Header2[$i];
			my $cov_norm = $tmp[$i];
			$Viral_gn2IMG2cov_norm_filtered{$viral_gn}{$img} = $cov_norm;
		}
	}
}
close IN;

## Step 5.2 Make %AMG_gene_containing_viral_gn2month2cov
my %AMG_gene_containing_viral_gn2month2cov = (); # $viral_gn => $month => $cov
foreach my $amg_gene (sort keys %AMG_gene2month2cov_ratio){
	my ($viral_gn) = $amg_gene =~ /^(.+?\_\_.+?)\_\_/;
	foreach my $month (sort keys %Month2img_id){
		my @IMG_ID = split (/\t/,$Month2img_id{$month}); # Store all metagenomes from this month
		my @Cov_collection = (); # Store all the non-"NA" cov values 
		my $cov_for_this_month = 0; # Store the cov for this month (the mean value of all non-"NA" values)
		foreach my $img (@IMG_ID){
			my $cov = $Viral_gn2IMG2cov_norm_filtered{$viral_gn}{$img};
			if ($cov ne "NA"){
				push @Cov_collection, $cov;
			}
		}
		
		if (@Cov_collection){
			$cov_for_this_month = mean(@Cov_collection);
		}
		$AMG_gene_containing_viral_gn2month2cov{$viral_gn}{$month} = $cov_for_this_month;	
	}
}

## Step 5.3 Write down %AMG_gene_containing_viral_gn2month2cov
open OUT, ">MetaPop/AMG_gene_containing_viral_gn2month2cov.txt";
my $row3=join("\t", sort keys %Month2num);
print OUT "Head\t$row3\n";
foreach my $tmp1 (sort keys %AMG_gene_containing_viral_gn2month2cov){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %Month2num) {       
                if (exists $AMG_gene_containing_viral_gn2month2cov{$tmp1}{$tmp2}){
                        push @tmp, $AMG_gene_containing_viral_gn2month2cov{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

## Step 5.4 Store %AMG_gene_containing_viral_gn2KOs and write it down
my %AMG_gene_containing_viral_gn2KOs = (); # $viral_gn => $kos (collection of $ko, separated by ";")
foreach my $amg_gene (sort keys %AMG_gene2month2cov_ratio){
	my ($viral_gn) = $amg_gene =~ /^(.+?\_\_.+?)\_\_/;
	my $ko = $AMG_summary{$amg_gene};
	if (!exists $AMG_gene_containing_viral_gn2KOs{$viral_gn}){
		$AMG_gene_containing_viral_gn2KOs{$viral_gn} = $ko;
	}else{
		if ($AMG_gene_containing_viral_gn2KOs{$viral_gn} !~ /$ko/){
			$AMG_gene_containing_viral_gn2KOs{$viral_gn} .= "\;".$ko;
		}
	}
}

open OUT, ">MetaPop/AMG_gene_containing_viral_gn2KOs.txt";
foreach my $viral_gn (sort keys %AMG_gene_containing_viral_gn2KOs){
	print OUT "$viral_gn\t$AMG_gene_containing_viral_gn2KOs{$viral_gn}\n";
}
close OUT;

## Step 5.5 Make %AMG_gene2ko_n_detail and write it down
my %KO2detail = (); # $ko => $detail
open IN, "MetaPop/KO2month_ko_details.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $ko = $tmp[0];
	my $detail = $tmp[1];
	$KO2detail{$ko} = $detail;
}
close IN;

my %AMG_gene2ko_n_detail = (); # $amg_gene => $ko_n_detail
foreach my $amg_gene (sort keys %AMG_gene2month2cov_ratio){
	my $ko = $AMG_summary{$amg_gene};
	my $detail = $KO2detail{$ko};
	my $ko_n_detail = "$ko\: $detail";
	$AMG_gene2ko_n_detail{$amg_gene} = $ko_n_detail;
}

open OUT, ">MetaPop/AMG_gene2ko_n_detail.txt";
foreach my $amg_gene (sort keys %AMG_gene2ko_n_detail){
	print OUT "$amg_gene\t$AMG_gene2ko_n_detail{$amg_gene}\n";
}
close OUT;

# Step 6 Make %AMG_gene_containing_viral_gn2year_month2cov hash and write it down
## Step 6.1 Make %AMG_gene_containing_viral_gn2year_month2cov
my %AMG_gene_containing_viral_gn2year_month2cov = (); # $viral_gn => $year_month => $cov
foreach my $amg_gene (sort keys %AMG_gene2year_month2cov_ratio){
	my ($viral_gn) = $amg_gene =~ /^(.+?\_\_.+?)\_\_/;
	foreach my $year_month (sort keys %Year_month2img_id){
		my @IMG_ID = split (/\t/,$Year_month2img_id{$year_month}); # Store all metagenomes from this year_month
		my @Cov_collection = (); # Store all the non-"NA" cov values 
		my $cov_for_this_year_month = 0; # Store the cov for this year_month (the mean value of all non-"NA" values)
		foreach my $img (@IMG_ID){
			my $cov = $Viral_gn2IMG2cov_norm_filtered{$viral_gn}{$img};
			if ($cov ne "NA"){
				push @Cov_collection, $cov;
			}
		}
		
		if (@Cov_collection){
			$cov_for_this_year_month = mean(@Cov_collection);
		}
		$AMG_gene_containing_viral_gn2year_month2cov{$viral_gn}{$year_month} = $cov_for_this_year_month;
	}
}

# Step 6.2 Write down the %AMG_gene_containing_viral_gn2year_month2cov hash 
open OUT, ">MetaPop/AMG_gene_containing_viral_gn2year_month2cov.txt";
my $row4=join("\t", @Year_month);
print OUT "Head\t$row4\n";
foreach my $tmp1 (sort keys %AMG_gene_containing_viral_gn2year_month2cov){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (@Year_month) {       
                if (exists $AMG_gene_containing_viral_gn2year_month2cov{$tmp1}{$tmp2}){
                        push @tmp, $AMG_gene_containing_viral_gn2year_month2cov{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;



# Subroutine
sub mean {
    return sum(@_)/@_;
}






