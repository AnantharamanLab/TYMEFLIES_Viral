#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw(sum);

# Aim: Calculate gene frequency for AMG-containing viral gn (species representatives) for each year

# Note:
# (1) Gene frequency was estimated as the coverage of each gene divided by the mean coverage of all other genes in the genome
# (2) The mean coverage of all other genes in the genome should be >= 5
# (3) The gene number in a genome with a valid coverage (not "NA" and > 0) should be over 50% of the total gene number
# (4) Genes with their length < 450 bp were excluded from analysis
# (5) The positions within the first and last 150 bp of a scaffold are excluded from depth/coverage analysis

# Step 1 Store %Viral_gene2year2depth hash
my %Viral_gene2year2depth = (); # $gene => $year => $gene_depth
my %Year = (); # $year => 1
open IN, "ls MetaPop.for_each_year/MetaPop/04.Depth_per_Pos/Gene_coverage_for_each_year/*.gene_coverage.txt |";
while (<IN>){
	chomp;
	my $file = $_;
	my ($year) = $file =~ /Gene\_coverage\_for\_each\_year\/(\d+?)\./;
	$Year{$year} = 1;
	
	open INN, "$file";
	while (<INN>){
		chomp;
		my @tmp = split (/\t/);
		my $gene = $tmp[0];
		my $gene_depth = $tmp[1];
		$Viral_gene2year2depth{$gene}{$year} = $gene_depth;
	}
	close INN;
}
close IN;

# Step 2 Calculate gene depth
## Step 2.1 Store the hash %Viral_gn2viral_scf and hash %Viral_scf2length
my %All_phage_species_rep_gn_containing_AMG_seq = ();
%All_phage_species_rep_gn_containing_AMG_seq = _store_seq("reference_fasta_for_metapop/All_phage_species_rep_gn_containing_AMG.fasta");

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

## Step 2.2 Store and write down all gene coordinates 
my %All_gene_coordinates = (); # $gene => $gene, $viral_scf, $viral_scf_length, $start, $stop, $gene_size connected by "\t"
my %Viral_scf2genes = (); # $viral_scf => collection of $gene\|$start\|$stop, separated by "\,"
open IN, "All_phage_species_rep_gn_containing_AMG.mdfed.genes"; # Note: Using mdfed gene files here 
while (<IN>){
	chomp;
	if (/^>/){
		my $line = $_;
		my ($gene, $start, $stop) = $line =~ /^>(.+?) \# (.+?) \# (.+?) \#/;
		my ($viral_scf) = $gene =~ /^(.+)\_/;
		my $viral_scf_length = $Viral_scf2length{$viral_scf};
		my $gene_size = $stop - $start;
		$All_gene_coordinates{$gene} = "$gene\t$viral_scf\t$viral_scf_length\t$start\t$stop\t$gene_size";
		
		my $gene_start_stop = "$gene\|$start\|$stop";
		if (!exists $Viral_scf2genes{$viral_scf}){
			$Viral_scf2genes{$viral_scf} = $gene_start_stop;
		}else{
			$Viral_scf2genes{$viral_scf} .= "\,".$gene_start_stop;
		}
	}
}
close IN;

# Step 3 Get viral gene 2 year 2 gene freq hash
my %Viral_gene2year2gene_freq = (); # $gene => $year => $gene_freq
foreach my $gene (sort keys %All_gene_coordinates){
	foreach my $year (sort keys %Year){
		my ($viral_gn) = $gene =~ /^(.+?\_\_.+?)\_\_/;
		my @Viral_scf = split (/\t/, $Viral_gn2viral_scf{$viral_gn});
		
		# The gene freq is the ratio of this gene's depth over the median depth of all other genes in the viral gn
		my $mean_depth_of_all_other_genes = 0;
		my @Depth_of_all_other_genes = ();
		
		my $mean_depth_of_all_genes = 0;
		my @Mean_depth = ();
		
		my @Genes_with_valid_depth = (); # Store the gene with not "NA" and > 0 gene coverage 
		my $num_genes_with_valid_depth = 0;
		my @Genes = (); # Store total genes
		my $num_genes = 0;
		
		foreach my $viral_scf (@Viral_scf){
			my @Gene_start_stop = split (/\,/, $Viral_scf2genes{$viral_scf});
			foreach my $gene_start_stop (@Gene_start_stop){
				my ($gene_all, $start, $stop) = $gene_start_stop =~ /^(.+?)\|(.+?)\|(.+?)$/;
				if ($Viral_gene2year2depth{$gene_all}{$year} ne "NA" and $gene_all ne $gene){
					push @Depth_of_all_other_genes, $Viral_gene2year2depth{$gene_all}{$year};
				}
				if ($Viral_gene2year2depth{$gene_all}{$year} ne "NA"){
					push @Mean_depth, $Viral_gene2year2depth{$gene_all}{$year};
				}
				if ($Viral_gene2year2depth{$gene_all}{$year} ne "NA" and $Viral_gene2year2depth{$gene_all}{$year}){
					push @Genes_with_valid_depth, $gene_all;
				}
				push @Genes, $gene_all;
			}
		}	
		
		if (@Depth_of_all_other_genes){
			$mean_depth_of_all_other_genes = _mean(@Depth_of_all_other_genes);
		}
		
		if (@Mean_depth){
			$mean_depth_of_all_genes = _mean(@Mean_depth);
		}
		
		$num_genes_with_valid_depth = scalar @Genes_with_valid_depth;
		$num_genes = scalar @Genes;
		
		my $ratio_genes_with_valid_depth = $num_genes_with_valid_depth / $num_genes; 
		
		my $gene_freq = "NA";
		# We set the following three requirements to gene freq that are not assigned as "NA":
		# (1) The gene coverage should not be "NA"
		# (2) The mean depth of all other genes should not be 0
		# (3) The mean depth of all genes within this genome should be >= 5
		# (4) The gene number in a genome with a valid coverage (not "NA" and > 0) should be over 50% of the total gene number
		
		if ($Viral_gene2year2depth{$gene}{$year} ne "NA" and $mean_depth_of_all_other_genes and $mean_depth_of_all_genes >= 5 and $ratio_genes_with_valid_depth >= 0.5){
			$gene_freq = $Viral_gene2year2depth{$gene}{$year} / $mean_depth_of_all_other_genes;
		}
		
		$Viral_gene2year2gene_freq{$gene}{$year} = $gene_freq;
	}
}

# Step 4 Write down Viral_gene2year2gene_freq.txt
open OUT, ">MetaPop.for_each_year/MetaPop/Viral_gene2year2gene_freq.txt";
my $row=join("\t", sort keys %Year);
print OUT "Viral gene\t$row\n";
foreach my $tmp1 (sort keys %Viral_gene2year2gene_freq){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %Year) {       
                if (exists $Viral_gene2year2gene_freq{$tmp1}{$tmp2}){
                        push @tmp, $Viral_gene2year2gene_freq{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"NA";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

# Step 5 Write down Viral_gene2year2gene_freq.for_four_AMGs.txt
## Step 5.1 Store the viral species containing four AMGs
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

# Step 5.2 Write down Viral_gene2year2gene_freq.for_four_AMGs.txt
open OUT, ">MetaPop.for_each_year/MetaPop/Viral_gene2year2gene_freq.for_four_AMGs.txt";
my $row2=join("\t", sort keys %Year);
print OUT "Viral gene\t$row2\n";
foreach my $tmp1 (sort keys %Viral_gene2year2gene_freq){
        my ($viral_gn) = $tmp1 =~ /^(.+?)\_\_Ga/;
		if (exists $Viral_species_containing_four_AMGs{$viral_gn}){
			print OUT $tmp1."\t";
			my @tmp = ();
			foreach my $tmp2 (sort keys %Year) {       
					if (exists $Viral_gene2year2gene_freq{$tmp1}{$tmp2}){
							push @tmp, $Viral_gene2year2gene_freq{$tmp1}{$tmp2};
					}else{              
							push @tmp,"NA";
					}
			}
			print OUT join("\t",@tmp)."\n";
		}
}
close OUT;



## Subroutine

sub _mean { # Get mean value within an array
    return sum(@_)/@_;
}

sub _median { # Get median value within an array
    my (@data) = sort { $a <=> $b } @_;
    if ( scalar(@data) % 2 ) {
        return ( $data[ @data / 2 ] );
    } else {
        my ( $upper, $lower );
        $lower = $data[ @data / 2 ];
        $upper = $data[ @data / 2 - 1 ];
        return ( _mean( $lower, $upper ) );
    }
}

sub _store_seq{ # Store seq; the head will be truncated from the first space if space(s) are present
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