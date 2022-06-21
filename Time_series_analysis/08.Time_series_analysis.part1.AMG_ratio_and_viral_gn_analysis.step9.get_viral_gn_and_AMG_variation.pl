#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw(sum);

# Aim: Parse to get the viral gn coverage and AMG coverage variation pattern across metagenomes
# Processing the viral gn presence with both coverage (>= 1.5) and breadth (>= 70%)

# Step 1 Store all scaffold length
my %Scf2length = (); # $scf => $length
my %Seq_all_scfs = _store_seq("reference_fasta_for_metapop/All_phage_species_rep_gn_containing_AMG.fasta");
foreach my $key (sort keys %Seq_all_scfs){
	my ($scf) = $key =~ /^>(.+?)$/;
	my $length = length($Seq_all_scfs{$key});
	$Scf2length{$scf} = $length;
}

# Step 2 Store the viral_gn 2 img 2 breadth hash
# Step 2.1 Store %Viral_gn2scfs
my %Viral_gn2scfs = (); # $viral_gn => $scf collection separated by "\t" (Only for rep gn with AMG genes)
foreach my $key (sort keys %Seq_all_scfs){
	my ($scf) = $key =~ /^>(.+?)$/;
	my ($viral_gn) = $scf =~ /^(.+?\_\_.+?)\_\_/;
	if (!exists $Viral_gn2scfs{$viral_gn}){
		$Viral_gn2scfs{$viral_gn} = $scf;
	}else{
		$Viral_gn2scfs{$viral_gn} .= "\t".$scf;
	}
}

# Step 2.2 Store %Scf2IMG2Breadth hash
my %Scf2IMG2Breadth = (); # $scf => $img => $breadth
my %IMG = (); # $img => 1
open IN, "ls MetaPop/03.Breadth_and_Depth/*.viral_species_rep.id90_breadth_and_depth.tsv|";
while (<IN>){
	chomp;
	my $file = $_;
	my ($img) = $file =~ /03\.Breadth_and_Depth\/(.+?)\.viral_species_rep/;
	$IMG{$img} = 1;
	open INN, "$file";
	while (<INN>){
		chomp;
		my @tmp = split (/\t/);
		my $scf = $tmp[0];
		my $breadth = $tmp[2];
		$Scf2IMG2Breadth{$scf}{$img} = $breadth;
	}
	close INN;
}
close IN;

# Step 2.3 Store %Viral_gn2IMG2Breadth hash
my %Viral_gn2IMG2Breadth = (); # $viral_gn => $img => $breadth
foreach my $viral_gn (sort keys %Viral_gn2scfs){
	foreach my $img (sort keys %IMG){
		my @Scfs = split(/\t/, $Viral_gn2scfs{$viral_gn});
		my $total_Breadth_multipled_by_length = 0;
		my $total_length = 0;
		foreach my $scf (@Scfs){
			my $breath = 0;
			if (exists $Scf2IMG2Breadth{$scf}{$img}){
				$breath = $Scf2IMG2Breadth{$scf}{$img};
			}
			my $length = $Scf2length{$scf};
			$total_Breadth_multipled_by_length += $breath * $length;
			$total_length += $length;
		}
		
		my $breadth_for_viral_gn = 0;
		if ($total_Breadth_multipled_by_length){
			$breadth_for_viral_gn = $total_Breadth_multipled_by_length / $total_length;
		}
		$Viral_gn2IMG2Breadth{$viral_gn}{$img} = $breadth_for_viral_gn;
	}
}

# Step 3 Store Viral_gn2IMG2cov_norm.txt and filter low coverage and low breadth viral gn
my %Viral_gn2IMG2cov_norm_filtered = (); # $viral_gn => $img => $cov_norm
my @Header = (); # Store the header line
open IN, "MetaPop/Viral_gn2IMG2cov_norm.txt";
while (<IN>){
	chomp;
	if (/^Head/){
		my $line = $_;
		@Header = split (/\t/, $line);
	}else{
		my $line = $_;
		my @tmp = split (/\t/, $line);
		my $viral_gn = $tmp[0];
		for(my $i=1; $i<=$#tmp; $i++){
			my $img = $Header[$i];
			my $cov_norm = $tmp[$i];
			my $breadth = $Viral_gn2IMG2Breadth{$viral_gn}{$img};
			if ($cov_norm >= 1.5 and $breadth >= 70){
				$Viral_gn2IMG2cov_norm_filtered{$viral_gn}{$img} = $cov_norm;
			}else{
				$Viral_gn2IMG2cov_norm_filtered{$viral_gn}{$img} = "NA";
			}
		}
	}
}
close IN;

## Write down %Viral_gn2IMG2cov_norm_filtered hash
open OUT, ">MetaPop/Viral_gn2IMG2cov_norm_filtered.txt";
my $row=join("\t", sort keys %IMG);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %Viral_gn2IMG2cov_norm_filtered){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %IMG) {       
                if (exists $Viral_gn2IMG2cov_norm_filtered{$tmp1}{$tmp2}){
                        push @tmp, $Viral_gn2IMG2cov_norm_filtered{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

## Get statistics of viral_gn distribution and average coverage across all metagenomes
my %Viral_gn2distribution_n_cov_mean = (); # $viral_gn => [0] $distribution [1] $cov_mean
foreach my $viral_gn (sort keys %Viral_gn2IMG2cov_norm_filtered){
	my $distribution = 0;
	my @Covs = (); # Store all non-zero coverages
	foreach my $img (sort keys %IMG){
		if (exists $Viral_gn2IMG2cov_norm_filtered{$viral_gn}{$img} and $Viral_gn2IMG2cov_norm_filtered{$viral_gn}{$img} ne "NA"){
			$distribution++;
			push @Covs, $Viral_gn2IMG2cov_norm_filtered{$viral_gn}{$img};
		}
	}
	my $cov_mean = 0;
	if (@Covs){
		$cov_mean = mean(@Covs);
	}
	
	$Viral_gn2distribution_n_cov_mean{$viral_gn}[0] = $distribution;
	$Viral_gn2distribution_n_cov_mean{$viral_gn}[1] = $cov_mean;
}

## Write down %Viral_gn2distribution_n_cov_mean hash
open OUT, ">MetaPop/Viral_gn2distribution_n_cov_mean.txt";
foreach my $viral_gn (sort keys %Viral_gn2distribution_n_cov_mean){
	print OUT "$viral_gn\t$Viral_gn2distribution_n_cov_mean{$viral_gn}[0]\t$Viral_gn2distribution_n_cov_mean{$viral_gn}[1]\n";
}
close OUT;

# Step 4 Store AMG_gene2IMG2cov_ratio.txt result and parse to get AMG distribution variation
## Step 4.1 Store AMG_gene2IMG2cov_ratio.txt result
my %AMG_gene2IMG2cov_ratio_filtered = (); # $amg_gene => $img => $cov_ratio; filter those values that are corresponding to filtered viral_gn coverage
my @Header2 = (); # Store the header line
open IN, "MetaPop/AMG_gene2IMG2cov_ratio.txt";
while (<IN>){
	chomp;
	if (/^Head/){
		my $line = $_;
		@Header2 = split (/\t/, $line);
	}else{
		my $line = $_;
		my @tmp = split (/\t/, $line);
		my $amg_gene = $tmp[0];
		for(my $i=1; $i<=$#tmp; $i++){
			my $img = $Header2[$i];
			my $cov_ratio = $tmp[$i];
			
			my ($viral_gn) = $amg_gene =~ /^(.+?\_\_.+?)\_\_/;
			my $viral_gn_cov = $Viral_gn2IMG2cov_norm_filtered{$viral_gn}{$img};
			if ($viral_gn_cov eq "NA"){
				$cov_ratio = "NA";
			}
			
			$AMG_gene2IMG2cov_ratio_filtered{$amg_gene}{$img} = $cov_ratio;
		}
	}
}
close IN;

## Step 4.2 Parse to get AMG cov ration variation
### Store AMG KO information
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

### Store AMG_KO2KO_details.txt
my %AMG_KO2KO_details = ();
open IN, "AMG_KO2KO_details.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $ko = $tmp[0];
	my $ko_detail = $tmp[1];
	$AMG_KO2KO_details{$ko} = $ko_detail;
}
close IN;

### Store %AMG_gene_cov_ratio_variation_table
my %AMG_gene_cov_ratio_variation_table = (); # $amg_gene => [0] KO ($ko) [1] KO detail ($ko_detail) [2] Viral gn ($viral_gn)
														#	[3] Distribution of viral gn ($distribution) [4] Coverage mean across all metagenomes ($cov_mean)
														#   [5] Distribution of AMG gene [6] - [470] 465 metagenomes
my @Header3 = ();														
foreach my $amg_gene (sort keys %AMG_gene2IMG2cov_ratio_filtered){
	my $ko = $AMG_summary{$amg_gene};
	my $ko_detail = $AMG_KO2KO_details{$ko};
	my ($viral_gn) = $amg_gene =~ /^(.+?\_\_.+?)\_\_/;
	my $distribution = $Viral_gn2distribution_n_cov_mean{$viral_gn}[0];
	my $cov_mean = $Viral_gn2distribution_n_cov_mean{$viral_gn}[1];
	$AMG_gene_cov_ratio_variation_table{$amg_gene}[0] = $ko;
	$AMG_gene_cov_ratio_variation_table{$amg_gene}[1] = $ko_detail;
	$AMG_gene_cov_ratio_variation_table{$amg_gene}[2] = $viral_gn;
	$AMG_gene_cov_ratio_variation_table{$amg_gene}[3] = $distribution;
	$AMG_gene_cov_ratio_variation_table{$amg_gene}[4] = $cov_mean;
	
	$Header3[0] = "KO";
	$Header3[1] = "KO detail";
	$Header3[2] = "Viral gn";
	$Header3[3] = "Distribution of viral gn";
	$Header3[4] = "Coverage mean across all metagenomes";
	$Header3[5] = "Distribution of AMG genes";
	
	my $amg_gene_distribution = 0; # The distribution of AMG genes
	
	my $i = 6; # Store the index of both @Header3 and %AMG_gene_cov_ratio_variation_table
	foreach my $img (sort keys %IMG){
		my $cov_ratio = $AMG_gene2IMG2cov_ratio_filtered{$amg_gene}{$img};
		$AMG_gene_cov_ratio_variation_table{$amg_gene}[$i] = $cov_ratio;
		if ($cov_ratio ne "NA"){
			$amg_gene_distribution++;
		}
		$Header3[$i] = $img; 
		$i++;
	}
	$AMG_gene_cov_ratio_variation_table{$amg_gene}[5] = $amg_gene_distribution;
}

### Write down %AMG_gene_cov_ratio_variation_table
open OUT, ">MetaPop/AMG_gene_cov_ratio_variation_table.txt";
my $header3 = join("\t",@Header3);
print OUT "Head\t$header3\n";
foreach my $amg_gene (sort keys %AMG_gene_cov_ratio_variation_table){
	my @Line1 = (); my @Line2 = ();
	push @Line1, $amg_gene;
	for(my $i=0; $i<=5; $i++){
		push @Line1, $AMG_gene_cov_ratio_variation_table{$amg_gene}[$i];
	}
	
	for(my $i=6; $i<=$#Header3; $i++){
		push @Line2, $AMG_gene_cov_ratio_variation_table{$amg_gene}[$i];
	}
	
	my $line1 = join("\t",@Line1);
	my $line2 = join("\t",@Line2);
	print OUT "$line1\t$line2\n";
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

# Subroutine
sub mean {
    return sum(@_)/@_;
}






