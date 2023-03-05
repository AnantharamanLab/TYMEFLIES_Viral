#!/usr/bin/perl

use strict;
use warnings;

# Aim: Get the glnA-containing viral gn coverage statistics (the distribution of viral gn should be >= 5 out of 465 metagenomes)

# Step 1 Get glnA-containing viral gn list (species representatives)
## Step 1.1 Store species info
my %Species = (); # $gn_rep => $gns 
open IN, "/storage1/data11/TYMEFLIES_phage/Cluster_phage_genomes/Species_level_vOTUs_cluster.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $gn_rep = $tmp[0];
	my $gns = $tmp[1];
	$Species{$gn_rep} = $gns;
}
close IN;

## Step 1.2 Store AMG KO information
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

### Change the old gene to new gene
my %Old_gene2new_gene_map = (); # $gene_old => $gene_new
open IN, "New_gene2old_gene_map.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $gene_new = $tmp[0]; my $gene_old = $tmp[1];
	$Old_gene2new_gene_map{$gene_old} = $gene_new;
}
close IN;

foreach my $pro (sort keys %AMG_summary){
	if ($Old_gene2new_gene_map{$pro}){
		my $ko = $AMG_summary{$pro};
		my $gene_new = $Old_gene2new_gene_map{$pro};
		delete $AMG_summary{$pro}; # Delete the old gene and its value
		$AMG_summary{$gene_new} = $ko; # Add the new gene and its value
	}
}

## Step 1.3 Get glnA-containing viral genome
my %GlnA_containing_viral_gn = (); # $gn => 1
foreach my $pro (sort keys %AMG_summary){
	my ($gn) = $pro =~ /^(.+?\_\_.+?)\_\_/;
	my $ko = $AMG_summary{$pro};
	if (exists $Species{$gn} and $ko eq "K01915"){ # If both this genome is species representative genome and this AMG is glnA
		$GlnA_containing_viral_gn{$gn} = 1;
	}
}

# Step 2 Store the hash of Viral_gn2distribution_n_cov_mean
my %Viral_gn2distribution_n_cov_mean = (); # viral_gn => [0] distribution [1] $cov_mean
open IN, "MetaPop/Viral_gn2distribution_n_cov_mean.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $viral_gn = $tmp[0];
	my $distribution = $tmp[1];
	my $cov_mean = $tmp[2];
	$Viral_gn2distribution_n_cov_mean{$viral_gn}[0] = $distribution;
	$Viral_gn2distribution_n_cov_mean{$viral_gn}[1] = $cov_mean;
}
close IN;

# Step 3 Get glnA-containing viral gn coverage statistics and write it down
# Step 3.1 Get distribution and coverage mean information
my %GlnA_containing_viral_gn2distribution_n_cov_mean = (); # $gn => [0] distribution [1] $cov_mean
foreach my $gn (sort keys %GlnA_containing_viral_gn){
	my $distribution = $Viral_gn2distribution_n_cov_mean{$gn}[0];
	my $cov_mean = $Viral_gn2distribution_n_cov_mean{$gn}[1];
	$GlnA_containing_viral_gn2distribution_n_cov_mean{$gn}[0] = $distribution;
	$GlnA_containing_viral_gn2distribution_n_cov_mean{$gn}[1] = $cov_mean;
}

# Step 3.2 Get tax and host tax information
my %Viral_gn2tax = (); # $viral_gn => [0] $tax [1] $tax_method
open IN, "/storage1/data11/TYMEFLIES_phage/Taxonomic_classification/Each_bin_tax_combined_result.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $viral_gn = $tmp[0];
	my $tax = $tmp[1];
	my $tax_method = $tmp[2];
	$Viral_gn2tax{$viral_gn}[0] = $tax;
	$Viral_gn2tax{$viral_gn}[1] = $tax_method;
}
close IN;

my %Viral_gn2host_tax = (); # $viral_gn => [0] $host_tax [1] $host_tax_method
open IN, "/storage1/data11/TYMEFLIES_phage/Host_prediction/Viral_gn2host_tax_final.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $viral_gn = $tmp[0];
	my $host_tax = $tmp[1];
	my $host_tax_method = $tmp[2];
	$Viral_gn2host_tax{$viral_gn}[0] = $host_tax;
	$Viral_gn2host_tax{$viral_gn}[1] = $host_tax_method;
}
close IN;

open OUT, ">MetaPop/GlnA_containing_viral_gn2distribution_n_cov_mean_n_tax_n_host_tax.txt";
foreach my $gn (sort keys %GlnA_containing_viral_gn2distribution_n_cov_mean){
	my $distribution = $GlnA_containing_viral_gn2distribution_n_cov_mean{$gn}[0];
	my $cov_mean = $GlnA_containing_viral_gn2distribution_n_cov_mean{$gn}[1];
	my $tax = "NA"; my $tax_method = "NA";
	my $host_tax = "NA"; my $host_tax_method = "NA";
	if (exists $Viral_gn2tax{$gn}[0]){
		$tax = $Viral_gn2tax{$gn}[0]; $tax_method = $Viral_gn2tax{$gn}[1];
	}
	if (exists $Viral_gn2host_tax{$gn}[0]){
		$host_tax = $Viral_gn2host_tax{$gn}[0]; $host_tax_method = $Viral_gn2host_tax{$gn}[1];
	}	
	print OUT "$gn\t$distribution\t$cov_mean\t$tax\t$tax_method\t$host_tax\t$host_tax_method\n";
}
close OUT;

# Step 4 Get glnA-containing viral gn 2 month 2 cov_ratio result and 2 year_month 2 cov_ratio result

## Store month and year month information
my @Month = ('01','02','03','04','05','06','07','08','09','10','11','12');
my @Year = ('2000','2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019');
my @Year_month = ();
foreach my $year (@Year){
	foreach my $month (@Month){
		my $year_month = "$year\-$month";
		push @Year_month, $year_month;
	}
}

## Step 4.1 Get glnA-containing viral gn 2 month 2 cov result
my %GlnA_containing_viral_gn2month2cov = (); # $glnA_containing_viral_gn => $month => $cov
my @Header = (); # Store the header line of AMG_gene_containing_viral_gn2month2cov.txt 
open IN, "MetaPop/AMG_gene_containing_viral_gn2month2cov.txt";
while (<IN>){
	chomp;
	if (/^Head/){
		my @tmp = split (/\t/);
		@Header = @tmp;
	}else{
		my @tmp = split (/\t/);
		my $viral_gn = $tmp[0];
		if (exists $GlnA_containing_viral_gn{$viral_gn} and $GlnA_containing_viral_gn2distribution_n_cov_mean{$viral_gn}[0] >= 5){	
			for(my $i=1; $i<=$#tmp; $i++){
				my $month = $Header[$i];
				my $cov = $tmp[$i];
				$GlnA_containing_viral_gn2month2cov{$viral_gn}{$month} = $cov;
			}
		}	
	}
}
close IN;

open OUT, ">MetaPop/GlnA_containing_viral_gn2month2cov.txt";
my $row=join("\t", @Month);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %GlnA_containing_viral_gn2month2cov){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (@Month) {       
                if (exists $GlnA_containing_viral_gn2month2cov{$tmp1}{$tmp2}){
                        push @tmp, $GlnA_containing_viral_gn2month2cov{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

## Step 4.2 Get glnA-containing viral gn 2 year_month 2 cov result
my %GlnA_containing_viral_gn2year_month2cov = (); # $glnA_containing_viral_gn => $year_month => $cov
my @Header2 = (); # Store the header line of AMG_gene_containing_viral_gn2year_month2cov.txt 
open IN, "MetaPop/AMG_gene_containing_viral_gn2year_month2cov.txt";
while (<IN>){
	chomp;
	if (/^Head/){
		my @tmp = split (/\t/);
		@Header2 = @tmp;
	}else{
		my @tmp = split (/\t/);
		my $viral_gn = $tmp[0];
		if (exists $GlnA_containing_viral_gn{$viral_gn} and $GlnA_containing_viral_gn2distribution_n_cov_mean{$viral_gn}[0] >= 5){	
			for(my $i=1; $i<=$#tmp; $i++){
				my $year_month = $Header2[$i];
				my $cov = $tmp[$i];
				$GlnA_containing_viral_gn2year_month2cov{$viral_gn}{$year_month} = $cov;
			}
		}	
	}
}
close IN;

open OUT, ">MetaPop/GlnA_containing_viral_gn2year_month2cov.txt";
my $row2=join("\t", @Year_month);
print OUT "Head\t$row2\n";
foreach my $tmp1 (sort keys %GlnA_containing_viral_gn2year_month2cov){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (@Year_month) {       
                if (exists $GlnA_containing_viral_gn2year_month2cov{$tmp1}{$tmp2}){
                        push @tmp, $GlnA_containing_viral_gn2year_month2cov{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

# Step 5 Get glnA AMG gene 2 month 2 cov result and 2 year_month 2 cov ratio result
## Step 5.1 Get glnA AMG gene 2 month 2 cov ratio result
my %GlnA_AMG_gene2month2cov_ratio = (); # $glnA_amg_gene => $month => $cov_ratio
my @Header3 = (); # Store the header line of AMG_gene2month2cov_ratio.txt
open IN, "MetaPop/AMG_gene2month2cov_ratio.txt";
while (<IN>){
	chomp;
	if (/^Head/){
		my @tmp = split (/\t/);
		@Header3 = @tmp;
	}else{
		my @tmp = split (/\t/);
		my $amg_gene = $tmp[0];
		my ($viral_gn) = $amg_gene =~ /^(.+?\_\_.+?)\_\_/;
		# Satisfy three conditions:
		# 1) $amg_gene is glnA 2) $amg_gene is in a glnA containing viral gn 3) this viral gn has distribution >= 5 out of 465 metagenomes
		if ($AMG_summary{$amg_gene} eq "K01915" and exists $GlnA_containing_viral_gn{$viral_gn} and $GlnA_containing_viral_gn2distribution_n_cov_mean{$viral_gn}[0] >= 5){	
			for(my $i=1; $i<=$#tmp; $i++){
				my $month = $Header3[$i];
				my $cov_ratio = $tmp[$i];
				$GlnA_AMG_gene2month2cov_ratio{$amg_gene}{$month} = $cov_ratio;
			}
		}	
	}
}
close IN;

open OUT, ">MetaPop/GlnA_AMG_gene2month2cov_ratio.txt";
my $row3=join("\t", @Month);
print OUT "Head\t$row3\n";
foreach my $tmp1 (sort keys %GlnA_AMG_gene2month2cov_ratio){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (@Month) {       
                if (exists $GlnA_AMG_gene2month2cov_ratio{$tmp1}{$tmp2}){
                        push @tmp, $GlnA_AMG_gene2month2cov_ratio{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

## Step 5.2 Get glnA AMG gene 2 year_month 2 cov ratio result
my %GlnA_AMG_gene2year_month2cov_ratio = (); # $glnA_amg_gene => $year_month => $cov_ratio
my @Header4 = (); # Store the header line of AMG_gene2year_month2cov_ratio.txt
open IN, "MetaPop/AMG_gene2year_month2cov_ratio.txt";
while (<IN>){
	chomp;
	if (/^Head/){
		my @tmp = split (/\t/);
		@Header4 = @tmp;
	}else{
		my @tmp = split (/\t/);
		my $amg_gene = $tmp[0];
		my ($viral_gn) = $amg_gene =~ /^(.+?\_\_.+?)\_\_/;
		# Satisfy three conditions:
		# 1) $amg_gene is glnA 2) $amg_gene is in a glnA containing viral gn 3) this viral gn has distribution >= 5 out of 465 metagenomes
		if ($AMG_summary{$amg_gene} eq "K01915" and exists $GlnA_containing_viral_gn{$viral_gn} and $GlnA_containing_viral_gn2distribution_n_cov_mean{$viral_gn}[0] >= 5){	
			for(my $i=1; $i<=$#tmp; $i++){
				my $year_month = $Header4[$i];
				my $cov_ratio = $tmp[$i];
				$GlnA_AMG_gene2year_month2cov_ratio{$amg_gene}{$year_month} = $cov_ratio;
			}
		}	
	}
}
close IN;

open OUT, ">MetaPop/GlnA_AMG_gene2year_month2cov_ratio.txt";
my $row4=join("\t", @Year_month);
print OUT "Head\t$row4\n";
foreach my $tmp1 (sort keys %GlnA_AMG_gene2year_month2cov_ratio){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (@Year_month) {       
                if (exists $GlnA_AMG_gene2year_month2cov_ratio{$tmp1}{$tmp2}){
                        push @tmp, $GlnA_AMG_gene2year_month2cov_ratio{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;