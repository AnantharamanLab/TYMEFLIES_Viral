#!/usr/bin/perl

use strict;
use warnings;

# Aim: Get the viral species (species rep gn) SNP allele frequencies for each year

# Step 1 Get AMG-containing viral gn list (species representatives)
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
my %New_gene2old_gene_map = (); # $gene_new => $gene_old
open IN, "New_gene2old_gene_map.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $gene_new = $tmp[0]; my $gene_old = $tmp[1];
	$Old_gene2new_gene_map{$gene_old} = $gene_new;
	$New_gene2old_gene_map{$gene_new} = $gene_old
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

# Step 2 Parse the raw SNP result table
my %All_SNP_pos2viral_gn = (); # $snp_pos => $viral_gn
my %Viral_gn2all_SNP_pos = (); # $viral_gn => collection of all SNP position (separated by "\,", and ordered)
my %SNP_pos2year2ATCG_num = (); # $snp_pos => $img => $atcg_num (for instance, 0,14,4,0)
my %Year = (); # $year => 1

open IN, "MetaPop.for_each_year/MetaPop/10.Microdiversity/global_raw_microdiversity_data_snp_loci_only.tsv";
while (<IN>){
	chomp;
	if (!/^contig/){
		my @tmp = split (/\t/);
		my $snp_pos = $tmp[0];
		my ($viral_gn) = $snp_pos =~ /^(.+?\_\_.+?)\_\_/;
		my ($year) = $tmp[9] =~ /^(.+?)\./;
		$Year{$year} = 1;
		$All_SNP_pos2viral_gn{$snp_pos} = $viral_gn;
		
		my $atcg_num = $tmp[5]."\,".$tmp[6]."\,".$tmp[7]."\,".$tmp[8];
		$SNP_pos2year2ATCG_num{$snp_pos}{$year} = $atcg_num;
	}
}
close IN;

## Store %Viral_gn2all_SNP_pos
foreach my $snp_pos (sort keys %All_SNP_pos2viral_gn){
	my $viral_gn = $All_SNP_pos2viral_gn{$snp_pos};
	if (!exists $Viral_gn2all_SNP_pos{$viral_gn}){
		$Viral_gn2all_SNP_pos{$viral_gn} = $snp_pos;
	}else{
		$Viral_gn2all_SNP_pos{$viral_gn} .= "\,".$snp_pos;
	}
}

## Reorder the SNP positions
foreach my $viral_gn (sort keys %Viral_gn2all_SNP_pos){
	my @SNP_pos = split (/\,/, $Viral_gn2all_SNP_pos{$viral_gn});
	@SNP_pos = _reorder_SNP_position(@SNP_pos);
	my $snp_pos_collection = join("\,", @SNP_pos);
	$Viral_gn2all_SNP_pos{$viral_gn} = $snp_pos_collection;
}

## Store the reference allele
my %All_SNP_pos2last_year2ref_allele = (); # Store the reference allele (the majority allele in the last year - 2018)
                                           # 2018 is the last metagenome
foreach my $snp_pos (sort keys %All_SNP_pos2viral_gn){
	if (!exists $SNP_pos2year2ATCG_num{$snp_pos}{"2018"}){
		$All_SNP_pos2last_year2ref_allele{$snp_pos} = "NA";
	}else{
		my $atcg_num = $SNP_pos2year2ATCG_num{$snp_pos}{"2018"};
		$All_SNP_pos2last_year2ref_allele{$snp_pos} = _pick_major_allele($atcg_num);
	}	
}

## Get the SNP allele frequency for each SNP position
my %All_SNP_pos2year2allele_freq = (); # $snp_pos => $year => $allele_freq
foreach my $snp_pos (sort keys %All_SNP_pos2viral_gn){
	foreach my $year (sort keys %Year){
		if (!exists $SNP_pos2year2ATCG_num{$snp_pos}{$year}){
			$All_SNP_pos2year2allele_freq{$snp_pos}{$year} = "NA";
		}else{
			my $atcg_num = $SNP_pos2year2ATCG_num{$snp_pos}{$year};
			my $ref_allele = $All_SNP_pos2last_year2ref_allele{$snp_pos};
			my $allele_freq = _get_allele_freq($ref_allele, $atcg_num);
			$All_SNP_pos2year2allele_freq{$snp_pos}{$year} = $allele_freq;
		}	
	}
}

open OUT, ">MetaPop.for_each_year/MetaPop/All_SNP_pos2year2allele_freq.txt";
my $row=join("\t", sort keys %Year);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %All_SNP_pos2viral_gn){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %Year) {       
                if (exists $All_SNP_pos2year2allele_freq{$tmp1}{$tmp2}){
                        push @tmp, $All_SNP_pos2year2allele_freq{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"NA";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;



# Subroutine
sub _pick_major_allele{
	my @ATCG_num = split (/\,/, $_[0]);
	
	my $sum_num = 0; my $largest_num_index = 0; my $largest_num = 0;
	for(my $i=0; $i<=$#ATCG_num; $i++){
		my $index = $i;
		my $num = $ATCG_num[$i];
		$sum_num += $num;
		
		if ($num > $largest_num){
			$largest_num = $num;
			$largest_num_index = $index;
		}
	}
	
	my $major_allele_returned = "NA";
	if ($sum_num){
		if ($largest_num_index == 0){
			$major_allele_returned = "A";
		}elsif($largest_num_index == 1){
			$major_allele_returned = "T";
		}elsif($largest_num_index == 2){
			$major_allele_returned = "C";
		}elsif($largest_num_index == 3){
			$major_allele_returned = "G";
		}
	}
	
	return $major_allele_returned;
}

sub _reorder_SNP_position{
	my @Origial_array = @_;
	
	my %SNP_pos_map = (); # $snp_pos_new => $snp_pos_original
	my @New_array = ();
	foreach my $key (@Origial_array){
		my ($key_new_1, $key_new_2) = $key =~ /^(.+)\_(\d+?)$/; # Get the first part and second part of a input key
		$key_new_2=(sprintf "%05d", $key_new_2); # Add 0s to the left to make it a 5-digit num
		my $key_new = $key_new_1."\_".$key_new_2;
		$SNP_pos_map{$key_new} = $key;
		push @New_array, $key_new;
	}

	my @New_array_sorted = sort @New_array;
	
	my @Original_array_sorted = ();
	for(my $i=0; $i<=$#New_array_sorted; $i++){
		my $snp_pos_new = $New_array_sorted[$i];
		my $snp_pos_original = $SNP_pos_map{$snp_pos_new};
		push @Original_array_sorted, $snp_pos_original;
	}
	
	return @Original_array_sorted;
}

sub _get_allele_freq{
	my $ref_allele_ = $_[0];
	my $atcg_num_ = $_[1];
	
	my $allele_freq_ = "NA";
	
	if ($ref_allele_ ne "NA"){
		my @ATCG_num_ = split (/\,/, $atcg_num_);
		if ($ref_allele_ eq "A"){
			$allele_freq_ = $ATCG_num_[0] / ($ATCG_num_[0] + $ATCG_num_[1] + $ATCG_num_[2] + $ATCG_num_[3]);
		}elsif ($ref_allele_ eq "T"){
			$allele_freq_ = $ATCG_num_[1] / ($ATCG_num_[0] + $ATCG_num_[1] + $ATCG_num_[2] + $ATCG_num_[3]);
		}elsif ($ref_allele_ eq "C"){
			$allele_freq_ = $ATCG_num_[2] / ($ATCG_num_[0] + $ATCG_num_[1] + $ATCG_num_[2] + $ATCG_num_[3]);
		}elsif ($ref_allele_ eq "G"){
			$allele_freq_ = $ATCG_num_[3] / ($ATCG_num_[0] + $ATCG_num_[1] + $ATCG_num_[2] + $ATCG_num_[3]);
		}
	}
	
	return $allele_freq_;
}

