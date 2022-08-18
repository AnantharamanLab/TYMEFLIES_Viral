#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw(sum);

# Aim: Parse the AMG coverage result
# Processing the coverage result with the following criteria:
# 1) screen bin with its any scaffold with < 0.01 coverage 
# 2) normalize the coverage by read numbers 
# 3) maybe set a cutoff to assign a bin is really present or just because of some random read hits to make it have a very low coverage (This was done in the next script)

# Step 1 Store all AMG coverage results
my %Scf2IMG2amg_gene2cov = (); # Store the coverage information of each AMG gene
                               # $scf => $img => $amg_gene => $amg_gene_cov
my %Viral_gn2scfs = (); # $viral_gn => $scf collection separated by "\t" (Only for rep gn with AMG genes)	
my %Scf2IMG2cov = (); # $scf => $img => $partial_scf_cov
my %Scf2amg_genes = (); # $scf => $amg_gene collection separated by "\t"
open IN, "ls MetaPop/AMG_coverage_result/*.viral_species_rep.id90.AMG_cov.txt|";
while (<IN>){
	chomp;
	my $file = $_;
	my ($img) = $file =~ /AMG_coverage_result\/(.+?)\.viral_species_rep\.id90\.AMG_cov\.txt/;
	open INN, "$file";
	while (<INN>){
		chomp;
		if (!/^scaffold/){
			my @tmp = split(/\t/);
			my $scf = $tmp[0];
			my $amg_gene = $tmp[1]; 
			my $whole_scf_cov = $tmp[2];
			my $partial_scf_cov = $tmp[3];
			my $amg_gene_cov = $tmp[4];
			if ($amg_gene_cov ne "NA"){
				$Scf2IMG2amg_gene2cov{$scf}{$img}{$amg_gene} = $amg_gene_cov; # $amg_gene_cov can be "0.0"
			}
			
			# Store %Viral_gn2scfs
			my ($viral_gn) = $scf =~ /^(.+?\_\_.+?)\_\_/;
			if (!exists $Viral_gn2scfs{$viral_gn}){
				$Viral_gn2scfs{$viral_gn} = $scf;
			}else{
				if ($Viral_gn2scfs{$viral_gn} !~ /$scf/){
					$Viral_gn2scfs{$viral_gn} .= "\t".$scf;
				}
			}	
			
			# Store %Scf2IMG2cov
			if ($partial_scf_cov ne "NA"){
				$Scf2IMG2cov{$scf}{$img} = $partial_scf_cov;
			}else{
				$Scf2IMG2cov{$scf}{$img} = $whole_scf_cov;
			}
			
			# Store %Scf2amg_genes
			if (!exists $Scf2amg_genes{$scf}){
				$Scf2amg_genes{$scf} = $amg_gene;
			}else{
				if ($amg_gene ne "NA" and $Scf2amg_genes{$scf} !~ /$amg_gene/){
					$Scf2amg_genes{$scf} .= "\t".$amg_gene;
				}
			}
		}
	}
	close INN;
}
close IN;

# Step 2 Make %Viral_gn2IMG2cov_norm
## Step 2.1 Store %IMG2read_num
my %IMG2read_num = (); # $img => $read_num
open IN, "Read_count_file_for_metapop.txt";
while (<IN>){
	chomp;
	my @tmp = split(/\t/);
	my ($img) = $tmp[0] =~ /^(\d+?)\./;
	my $read_num = $tmp[1];
	$IMG2read_num{$img} = $read_num;
}
close IN;

## Step 2.2 Make %Viral_gn2IMG2cov_norm
my %Viral_gn2IMG2cov_norm = (); # $viral_gn => $img => $cov_norm
foreach my $viral_gn (sort keys %Viral_gn2scfs){
	foreach my $img (sort keys %IMG2read_num){
		my @Scfs = split(/\t/, $Viral_gn2scfs{$viral_gn});
		my @Scf_cov_norm_collection = ();
		foreach my $scf (@Scfs){
			my $scf_cov = $Scf2IMG2cov{$scf}{$img};
			my $scf_cov_norm = $scf_cov * (100000000 / $IMG2read_num{$img}); # Normalized cov, normalized by 100M reads per metagenome
			push @Scf_cov_norm_collection, $scf_cov_norm;
		}
		
		# Filter viral_gn if its any scaffold coverge is < 0.01
		my $logic = 1; # Preset all scaffold coverage is >= 0.01
		foreach my $key (@Scf_cov_norm_collection){
			if ($key < 0.01){
				$logic = 0;
			}
		}
		
		if ($logic){
			my $scf_cov_norm_mean = mean(@Scf_cov_norm_collection);
			$Viral_gn2IMG2cov_norm{$viral_gn}{$img} = $scf_cov_norm_mean;
		}
	}
}

## Step 2.3 Write down %Viral_gn2IMG2cov_norm
open OUT, ">MetaPop/Viral_gn2IMG2cov_norm.txt";
my $row=join("\t", sort keys %IMG2read_num);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %Viral_gn2IMG2cov_norm){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %IMG2read_num) {       
                if (exists $Viral_gn2IMG2cov_norm{$tmp1}{$tmp2}){
                        push @tmp, $Viral_gn2IMG2cov_norm{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"0"; # If a img does not have viral gn, then the coverage is 0
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

# Step 3 Make %Viral_gn2IMG2amg_gene2cov_norm
my %Viral_gn2IMG2amg_gene2cov_norm = (); # $viral_gn => $img => $amg_gene => $cov_norm
my %Viral_gn2amg_genes = (); # $viral_gn => $amg_gene collection separated by "\t"
foreach my $viral_gn (sort keys %Viral_gn2IMG2cov_norm){
	foreach my $img (sort keys %IMG2read_num){
		my @Scfs = split(/\t/, $Viral_gn2scfs{$viral_gn});
		foreach my $scf (@Scfs){
			if (exists $Scf2amg_genes{$scf}){
				my @AMG_genes = split (/\t/, $Scf2amg_genes{$scf});
				foreach my $amg_gene (@AMG_genes){
					if (exists $Scf2IMG2amg_gene2cov{$scf}{$img}{$amg_gene}){
						my $amg_gene_cov = $Scf2IMG2amg_gene2cov{$scf}{$img}{$amg_gene};
						my $amg_gene_cov_norm = $amg_gene_cov * (100000000 / $IMG2read_num{$img});
						$Viral_gn2IMG2amg_gene2cov_norm{$viral_gn}{$img}{$amg_gene} = $amg_gene_cov_norm;
						
						# Store %Viral_gn2amg_genes
						if (!exists $Viral_gn2amg_genes{$viral_gn}){
							$Viral_gn2amg_genes{$viral_gn} = $amg_gene;
						}else{
							if ($Viral_gn2amg_genes{$viral_gn} !~ /$amg_gene/){
								$Viral_gn2amg_genes{$viral_gn} .= "\t".$amg_gene;
							}
						}
					}
				}
			}
		}
	}
}

# Step 4 Make %AMG_gene2IMG2cov_ratio and write it down
my %AMG_gene2IMG2cov_ratio = (); # $amg_gene => $img => $cov_ratio ($amg_gene coverage / $viral_gn coverage)
foreach my $viral_gn (sort keys %Viral_gn2amg_genes){
	my @AMG_genes = split(/\t/, $Viral_gn2amg_genes{$viral_gn});
	foreach my $amg_gene (@AMG_genes){
		foreach my $img (sort keys %IMG2read_num){
			if (exists $Viral_gn2IMG2amg_gene2cov_norm{$viral_gn}{$img}{$amg_gene} and exists $Viral_gn2IMG2cov_norm{$viral_gn}{$img}){
				my $amg_gene_cov = $Viral_gn2IMG2amg_gene2cov_norm{$viral_gn}{$img}{$amg_gene};
				my $viral_gn_cov = $Viral_gn2IMG2cov_norm{$viral_gn}{$img};
				my $cov_ratio = $amg_gene_cov / $viral_gn_cov;
				$AMG_gene2IMG2cov_ratio{$amg_gene}{$img} = $cov_ratio;	
			}elsif(exists $Viral_gn2IMG2amg_gene2cov_norm{$viral_gn}{$img}{$amg_gene} and !exists $Viral_gn2IMG2cov_norm{$viral_gn}{$img}){
				my $cov_ratio = "Viral gn absent, AMG present";
				$AMG_gene2IMG2cov_ratio{$amg_gene}{$img} = $cov_ratio;	
			}elsif(!exists $Viral_gn2IMG2amg_gene2cov_norm{$viral_gn}{$img}{$amg_gene} and exists $Viral_gn2IMG2cov_norm{$viral_gn}{$img}){
				my $cov_ratio = "Viral gn present, AMG absent"; # Does not exit actually
				$AMG_gene2IMG2cov_ratio{$amg_gene}{$img} = $cov_ratio;	
			}elsif(!exists $Viral_gn2IMG2amg_gene2cov_norm{$viral_gn}{$img}{$amg_gene} and !exists $Viral_gn2IMG2cov_norm{$viral_gn}{$img}){
				my $cov_ratio = "Both viral gn and AMG absent";
				$AMG_gene2IMG2cov_ratio{$amg_gene}{$img} = $cov_ratio;					
			}
		}
	}
}

open OUT, ">MetaPop/AMG_gene2IMG2cov_ratio.txt";
my $row2=join("\t", sort keys %IMG2read_num);
print OUT "Head\t$row2\n";
foreach my $tmp1 (sort keys %AMG_gene2IMG2cov_ratio){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %IMG2read_num) {       
                if (exists $AMG_gene2IMG2cov_ratio{$tmp1}{$tmp2}){
                        push @tmp, $AMG_gene2IMG2cov_ratio{$tmp1}{$tmp2};
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
