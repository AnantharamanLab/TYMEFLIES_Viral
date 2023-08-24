#!/usr/bin/perl
  
use strict;
use warnings;

# Aim: Get the viral gn 2 IMG 2 cov_norm filtered table for four AMGs

# Step 1 Store the viral species containing four AMGs
my %Viral_species_containing_four_AMGs = (); # $viral_gn => $info (the corresponding information for each viral species gn)
my @Viral_species_containing_four_AMGs = (); # Store the $viral_gn order
open IN, "viral_species_containing_four_AMGs.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $info = $tmp[0];
	my $viral_gn = $tmp[1];
	$Viral_species_containing_four_AMGs{$viral_gn} = $info;
	
	push @Viral_species_containing_four_AMGs, $viral_gn;
}
close IN;

# Step 2 Store Viral_gn2IMG2cov_norm_filtered.txt
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

# Step 3 Store AMG KO information
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

## Change the old gene to new gene
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

# Step 4 Write down Viral_gn2IMG2cov_norm_filtered.four_AMGs.txt
open OUT, ">MetaPop/Viral_gn2IMG2cov_norm_filtered.four_AMGs.txt";
my $row=join("\t", sort keys %IMG2date);
print OUT "Head\tViral species\t$row\n";
foreach my $tmp1 (@Viral_species_containing_four_AMGs){
        print OUT $Viral_species_containing_four_AMGs{$tmp1}."\t".$tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %IMG2date) {
                if (exists $Viral_gn2IMG2cov_norm_filtered{$tmp1}{$tmp2}){
                        push @tmp, $Viral_gn2IMG2cov_norm_filtered{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"NA";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;


