#!/usr/bin/perl

use strict;
use warnings;

# Aim: Get the pmoC-containing viral gn coverage statistics (the distribution of viral gn should be >= 1 out of 20 years) for each year

# Step 1 Get pmoC-containing viral gn list (species representatives)
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

## Step 1.3 Get pmoC-containing viral genome
my %PmoC_containing_viral_gn = (); # $gn => 1
foreach my $pro (sort keys %AMG_summary){
	my ($gn) = $pro =~ /^(.+?\_\_.+?)\_\_/;
	my $ko = $AMG_summary{$pro};
	if (exists $Species{$gn} and $ko eq "K10946"){ # If both this genome is species representative genome and this AMG is pmoC
		$PmoC_containing_viral_gn{$gn} = 1;
	}
}

# Step 2 Store the hash of Viral_gn2distribution_n_cov_mean
my %Viral_gn2distribution_n_cov_mean = (); # viral_gn => [0] distribution [1] $cov_mean
open IN, "MetaPop.for_each_year/MetaPop/Viral_gn2distribution_n_cov_mean.txt";
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

# Step 3 Get pmoC-containing viral gn coverage statistics and write it down
# Step 3.1 Get distribution and coverage mean information
my %PmoC_containing_viral_gn2distribution_n_cov_mean = (); # $gn => [0] distribution [1] $cov_mean
foreach my $gn (sort keys %PmoC_containing_viral_gn){
	my $distribution = $Viral_gn2distribution_n_cov_mean{$gn}[0];
	my $cov_mean = $Viral_gn2distribution_n_cov_mean{$gn}[1];
	$PmoC_containing_viral_gn2distribution_n_cov_mean{$gn}[0] = $distribution;
	$PmoC_containing_viral_gn2distribution_n_cov_mean{$gn}[1] = $cov_mean;
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

open OUT, ">MetaPop.for_each_year/MetaPop/PmoC_containing_viral_gn2distribution_n_cov_mean_n_tax_n_host_tax.txt";
foreach my $gn (sort keys %PmoC_containing_viral_gn2distribution_n_cov_mean){
	my $distribution = $PmoC_containing_viral_gn2distribution_n_cov_mean{$gn}[0];
	my $cov_mean = $PmoC_containing_viral_gn2distribution_n_cov_mean{$gn}[1];
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

# Step 4 Get pmoC-containing viral gn 2 year 2 cov_ratio result 
## Store the year information
my @Year = ('2000','2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019');

my %PmoC_containing_viral_gn2year2cov = (); # $pmoC_containing_viral_gn => $year => $cov
my @Header = (); # Store the header line of Viral_gn2year2cov_norm.txt
open IN, "MetaPop.for_each_year/MetaPop/Viral_gn2year2cov_norm.txt";
while (<IN>){
	chomp;
	if (/^Head/){
		my @tmp = split (/\t/);
		@Header = @tmp;
	}else{
		my @tmp = split (/\t/);
		my $viral_gn = $tmp[0];
		if (exists $PmoC_containing_viral_gn{$viral_gn} and $PmoC_containing_viral_gn2distribution_n_cov_mean{$viral_gn}[0] >= 1){	
			for(my $i=1; $i<=$#tmp; $i++){
				my $year = $Header[$i];
				my $cov = $tmp[$i];
				$PmoC_containing_viral_gn2year2cov{$viral_gn}{$year} = $cov;
			}
		}	
	}
}
close IN;

open OUT, ">MetaPop.for_each_year/MetaPop/PmoC_containing_viral_gn2year2cov.txt";
my $row=join("\t", @Year);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %PmoC_containing_viral_gn2year2cov){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (@Year) {       
                if (exists $PmoC_containing_viral_gn2year2cov{$tmp1}{$tmp2}){
                        push @tmp, $PmoC_containing_viral_gn2year2cov{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

# Step 5 Get pmoC AMG gene 2 year 2 cov ratio result
## Step 5.1 Get pmoC AMG gene 2 year 2 cov ratio result
my %PmoC_AMG_gene2year2cov_ratio = (); # $pmoC_amg_gene => $year => $cov_ratio
my @Header3 = (); # Store the header line of AMG_gene2year2cov_ratio.txt
open IN, "MetaPop.for_each_year/MetaPop/AMG_gene2year2cov_ratio.txt";
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
		# 1) $amg_gene is pmoC 2) $amg_gene is in a pmoC containing viral gn 3) this viral gn has distribution >= 1 out of 20 years
		if ($AMG_summary{$amg_gene} eq "K10946" and exists $PmoC_containing_viral_gn{$viral_gn} and $PmoC_containing_viral_gn2distribution_n_cov_mean{$viral_gn}[0] >= 1){	
			for(my $i=1; $i<=$#tmp; $i++){
				my $year = $Header3[$i];
				my $cov_ratio = $tmp[$i];
				if ($cov_ratio eq "Viral gn absent, AMG present" or $cov_ratio eq "Both viral gn and AMG absent"){
					$cov_ratio = 0;
					$PmoC_AMG_gene2year2cov_ratio{$amg_gene}{$year} = $cov_ratio;
				}else{
					$PmoC_AMG_gene2year2cov_ratio{$amg_gene}{$year} = $cov_ratio;
				}
			}
		}	
	}
}
close IN;

open OUT, ">MetaPop.for_each_year/MetaPop/PmoC_AMG_gene2year2cov_ratio.txt";
my $row3=join("\t", @Year);
print OUT "Head\t$row3\n";
foreach my $tmp1 (sort keys %PmoC_AMG_gene2year2cov_ratio){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (@Year) {       
                if (exists $PmoC_AMG_gene2year2cov_ratio{$tmp1}{$tmp2}){
                        push @tmp, $PmoC_AMG_gene2year2cov_ratio{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;