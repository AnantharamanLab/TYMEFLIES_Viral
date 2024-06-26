#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw(sum);

# Aim: 1) Get the AMG-containing viral gn microdiversity statistics 
# (the distribution of viral gn should be >= 20 out of 471 metagenomes for calculating season and year_season distribution)
# 2) Get the AMG-containing viral gn microdiversity statistics for four AMGs

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

# Step 3 Store season metagenome information
my %Season2num = (); # $season => $num_metagenome; Store how many samples (metagenomes) are there in each season for the whole datasets
my %Season2img_id = (); # $season => collection of $img_id separeted by "\t"
my %Year_season2num = (); # $year_season => $num_metagenome
my %Year_season2img_id = (); # $year_season (for example, "2010-06") => collection of $img_id separeted by "\t"
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

## Store season and year season information
my @Season = ('Spring', 'Clearwater', 'Early Summer', 'Late Summer', 'Fall', 'Ice-on');
my @Year = ('2000','2001','2002','2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019');
my @Year_season = ();
foreach my $year (@Year){
	foreach my $season (@Season){
		my $year_season = "$year\-$season";
		push @Year_season, $year_season;
	}
}

# Step 4 Store the viral species containing four AMGs
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

# Step 5 Store the nucleotide diversity of each viral genome
## Step 5.1 Store the hash %Viral_gn2viral_scf and hash %Viral_scf2length
my %All_phage_species_rep_gn_containing_AMG_seq = ();
%All_phage_species_rep_gn_containing_AMG_seq = _store_seq("All_phage_species_rep_gn_containing_AMG.fasta");

my %Viral_gn2viral_scf = (); # $viral_gn => collection of $viral_scf, separated by "\t"
my %Viral_scf2length = (); # $viral_scf => $length
foreach my $key (sort keys %All_phage_species_rep_gn_containing_AMG_seq){
	my ($viral_scf) = $key =~ /^>(.+?)$/;
	my ($viral_gn) = $viral_scf =~ /^(.+?\_\_vRhyme.+?)\_\_/;
	if (!exists $Viral_gn2viral_scf{$viral_gn}){
		$Viral_gn2viral_scf{$viral_gn} = $viral_scf;
	}else{
		$Viral_gn2viral_scf{$viral_gn} .= "\t".$viral_scf;
	}
	
	my $seq = $All_phage_species_rep_gn_containing_AMG_seq{$key}; # The sequence of the scaffold
	
	my $length = length($seq);
	$Viral_scf2length{$viral_scf} = $length;
}

## Step 5.2 Filter scaffolds with low depth and breadth
## Only scaffolds that pass filtering will be included in microdiversity analysis
my %Viral_scafflold2IMG2pass_filtering = (); # $viral_scf => $img => 1 (If pass filtering)
open IN, "ls MetaPop/03.Breadth_and_Depth/*.tsv |";
while (<IN>){
	chomp;
	my $file = $_;
	my ($img) = $file =~ /03\.Breadth\_and\_Depth\/(.+?)\.viral.+?\.tsv/;
	open INN, "$file";
	while (<INN>){
		chomp;
		my @tmp = split (/\t/);
		my $viral_scf = $tmp[0];
		my $breath = $tmp[2] / 100;
		my $depth = $tmp[3];
		
		if ($breath >= 0.5 and $depth >= 5){
			$Viral_scafflold2IMG2pass_filtering{$viral_scf}{$img} = 1;
		}
	}
	close INN;
}
close IN;

## Step 5.3 Calculate the nucleotide diversity of each viral genome and write it down
my %Viral_scaffold2IMG2nt_diversity = (); # $viral_scf => $img (the mapping reads source) => $nt_diversity (pi)
my %Viral_gn2IMG2nt_diversity = (); # $viral_gn => $img (the mapping reads source) => $nt_diversity (pi)                                									
									# The $nt_diversity was averaged by scaffold length
open IN, "MetaPop/10.Microdiversity/global_contig_microdiversity.tsv";
while (<IN>){
	chomp;
	if (!/^contig/){
		my $line = $_;
		my @tmp = split (/\t/, $line);
		my $viral_scf = $tmp[0];
		my ($img) = $tmp[1] =~ /^(\d+?)\.viral\_species\_rep\.id90/;
		my $nt_diversity = $tmp[2];
		
		if ($line =~ /TRUE/){
			$Viral_scaffold2IMG2nt_diversity{$viral_scf}{$img} = $nt_diversity;
		}
	}
}
close IN;

foreach my $viral_gn (sort keys %Viral_gn2viral_scf){
	my @Viral_scf = split (/\t/, $Viral_gn2viral_scf{$viral_gn});	
	foreach my $img (sort keys %IMG2date){
		my $logic = 1; # To check if all the scaffolds within the $viral_gn have a valid nt diversity and pass depth and breadth filtering
		
		my $viral_gn_length = 0;
		foreach my $viral_scf (@Viral_scf){
			if (!exists $Viral_scaffold2IMG2nt_diversity{$viral_scf}{$img} or !exists $Viral_scafflold2IMG2pass_filtering{$viral_scf}{$img}){
				$logic = 0;
			}
			$viral_gn_length += $Viral_scf2length{$viral_scf};
		}
		
		if ($logic){ # If all the scaffolds within the $viral_gn have a valid nt diversity
			my $nt_diversity_for_this_gn_this_img = 0;
			foreach my $viral_scf (@Viral_scf){
				my $nt_diversity_for_this_scf_for_this_img = $Viral_scaffold2IMG2nt_diversity{$viral_scf}{$img};
				my $nt_diversity_for_this_scf_for_this_img_normalized_by_length = $nt_diversity_for_this_scf_for_this_img * $Viral_scf2length{$viral_scf} / $viral_gn_length;
				$nt_diversity_for_this_gn_this_img += $nt_diversity_for_this_scf_for_this_img_normalized_by_length;
			}
			$Viral_gn2IMG2nt_diversity{$viral_gn}{$img} = $nt_diversity_for_this_gn_this_img;
		}
	}
}	

## Write down %Viral_gn2IMG2nt_diversity hash
open OUT, ">MetaPop/Viral_gn2IMG2nt_diversity.txt";
my $row=join("\t", sort keys %IMG2date);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %Viral_gn2viral_scf){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %IMG2date) {       
                if (exists $Viral_gn2IMG2nt_diversity{$tmp1}{$tmp2}){
                        push @tmp, $Viral_gn2IMG2nt_diversity{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"NA";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

## Step 5.4 Write down Viral_gn2IMG2nt_diversity.four_AMGs.txt
open OUT, ">MetaPop/Viral_gn2IMG2nt_diversity.four_AMGs.txt";
my $row4=join("\t", sort keys %IMG2date);
print OUT "Head\tViral species\t$row4\n";
foreach my $tmp1 (@Viral_species_containing_four_AMGs){
        print OUT $Viral_species_containing_four_AMGs{$tmp1}."\t".$tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %IMG2date) {       
                if (exists $Viral_gn2IMG2nt_diversity{$tmp1}{$tmp2}){
                        push @tmp, $Viral_gn2IMG2nt_diversity{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"NA";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

# Step 6 Get Viral gn 2 season 2 nucleotide diversity result and 2 year_season 2 nucleotide diversity result
## Step 6.1 Make %Viral_gn2year_season2nt_diversity
my %Viral_gn2year_season2nt_diversity = (); # $viral_gn => $year_season => $nt_diversity
foreach my $viral_gn (sort keys %Viral_gn2distribution_n_cov_mean){
	my $distribution = $Viral_gn2distribution_n_cov_mean{$viral_gn}[0];
	if ($distribution >= 20){
		foreach my $year_season (sort keys %Year_season2img_id){
			my @IMG_ID = split (/\t/,$Year_season2img_id{$year_season}); # Store all metagenomes from this year_season
			my @Nt_diversity_collection = (); # Store all the non-"NA" nt diversity values 
			my $nt_diversity_for_this_year_season = "NA"; # Store the nt diversity for this year_season (the mean value of all non-"NA" values)
			foreach my $img (@IMG_ID){
				my $nt_diversity = "NA";
				if (exists $Viral_gn2IMG2nt_diversity{$viral_gn}{$img}){
					$nt_diversity = $Viral_gn2IMG2nt_diversity{$viral_gn}{$img};
				}
				if ($nt_diversity ne "NA"){
					push @Nt_diversity_collection, $nt_diversity;
				}
			}
			
			if (@Nt_diversity_collection){
				$nt_diversity_for_this_year_season = mean(@Nt_diversity_collection);
			}
			$Viral_gn2year_season2nt_diversity{$viral_gn}{$year_season} = $nt_diversity_for_this_year_season;
		}
	}
}

## Step 6.2 Write down %Viral_gn2year_season2nt_diversity
open OUT, ">MetaPop/Viral_gn2year_season2nt_diversity.txt";
my $row2=join("\t", @Year_season);
print OUT "Head\t$row2\n";
foreach my $tmp1 (sort keys %Viral_gn2year_season2nt_diversity){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (@Year_season) {       
                if (exists $Viral_gn2year_season2nt_diversity{$tmp1}{$tmp2}){
                        push @tmp, $Viral_gn2year_season2nt_diversity{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"NA";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

## Step 6.3 Write down Viral_gn2year_season2nt_diversity.four_AMGs.txt
open OUT, ">MetaPop/Viral_gn2year_season2nt_diversity.four_AMGs.txt";
my $row5=join("\t", @Year_season);
print OUT "Head\tViral species\t$row5\n";
foreach my $tmp1 (@Viral_species_containing_four_AMGs){
        print OUT $Viral_species_containing_four_AMGs{$tmp1}."\t".$tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (@Year_season) {       
                if (exists $Viral_gn2year_season2nt_diversity{$tmp1}{$tmp2}){
                        push @tmp, $Viral_gn2year_season2nt_diversity{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"NA";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

## Step 6.4 Make %Viral_gn2season2nt_diversity
my %Viral_gn2season2nt_diversity = (); # $viral_gn => $season => $nt_diversity
foreach my $viral_gn (sort keys %Viral_gn2distribution_n_cov_mean){
	my $distribution = $Viral_gn2distribution_n_cov_mean{$viral_gn}[0];
	if ($distribution >= 20){
		foreach my $season (sort keys %Season2img_id){
			my @IMG_ID = split (/\t/,$Season2img_id{$season}); # Store all metagenomes from this season
			my @Nt_diversity_collection = (); # Store all the non-"NA" nt diversity values 
			my $nt_diversity_for_this_season = "NA"; # Store the nt diversity for this season (the mean value of all non-"NA" values)
			foreach my $img (@IMG_ID){
				my $nt_diversity = "NA";
				if (exists $Viral_gn2IMG2nt_diversity{$viral_gn}{$img}){
					$nt_diversity = $Viral_gn2IMG2nt_diversity{$viral_gn}{$img};
				}
				if ($nt_diversity ne "NA"){
					push @Nt_diversity_collection, $nt_diversity;
				}
			}
			
			if (@Nt_diversity_collection){
				$nt_diversity_for_this_season = mean(@Nt_diversity_collection);
			}
			$Viral_gn2season2nt_diversity{$viral_gn}{$season} = $nt_diversity_for_this_season;
		}
	}
}

## Step 6.5 Write down %Viral_gn2season2nt_diversity
open OUT, ">MetaPop/Viral_gn2season2nt_diversity.txt";
my $row3=join("\t", sort keys %Season2num);
print OUT "Head\t$row3\n";
foreach my $tmp1 (sort keys %Viral_gn2season2nt_diversity){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %Season2num) {       
                if (exists $Viral_gn2season2nt_diversity{$tmp1}{$tmp2}){
                        push @tmp, $Viral_gn2season2nt_diversity{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"NA";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

## Step 6.6 Write down Viral_gn2season2nt_diversity.four_AMGs.txt
open OUT, ">MetaPop/Viral_gn2season2nt_diversity.four_AMGs.txt";
my $row6=join("\t", sort keys %Season2num);
print OUT "Head\tViral species\t$row6\n";
foreach my $tmp1 (@Viral_species_containing_four_AMGs){
        print OUT $Viral_species_containing_four_AMGs{$tmp1}."\t".$tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %Season2num) {       
                if (exists $Viral_gn2season2nt_diversity{$tmp1}{$tmp2}){
                        push @tmp, $Viral_gn2season2nt_diversity{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"NA";
                }
        }
        print OUT join("\t",@tmp)."\n";
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

