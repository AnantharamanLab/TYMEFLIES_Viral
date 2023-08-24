#!/usr/bin/perl

use strict;
use warnings;

# Aim: Get the viral species (four AMG containing viral species) pNpS results for 2000-2003 and 2016-2019, and compare them to see if there are genes with elevated pNpS values

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

# Step 2 Get all viral species rep gene annotation
## Step 2.1 Store all the gene ID
my %All_gene_seq = _store_seq("All_phage_species_rep_gn_containing_AMG.mdfed.genes");
my %All_gene_seq_ID = (); # Store only the gene ID
foreach my $key (sort keys %All_gene_seq){
	my ($gene) = $key =~ /^>(.+?)$/;
	$All_gene_seq_ID{$gene} = 1;
}

my %All_gene_seq_ID_map = (); # $gene_short_wo_fragment_x => $gene
foreach my $gene (sort keys %All_gene_seq_ID){
	my $gene_new = $gene;
	
	my $gene_old = $gene;
	# Convert gene to gene_old 
	if (exists $New_gene2old_gene_map{$gene}){
		$gene_old = $New_gene2old_gene_map{$gene};
	}
	# Get the $gene_short
	my ($gene_short) = $gene_old =~ /^.+?\_\_.+?\_\_(.+?)$/;
	# Get the $gene_short_wo_fragment_x
	my $gene_short_wo_fragment_x = $gene_short;
	if ($gene_short_wo_fragment_x =~ /fragment\_\d\_/){
		$gene_short_wo_fragment_x =~ s/fragment\_\d\_//g;
	}
	
	$All_gene_seq_ID_map{$gene_short_wo_fragment_x} = $gene_new;
}

## Step 2.2 Store all gene annotation results
my $annotation_header = "";
my %All_gene_annotation = (); # $gene => $annotation
open IN, "ls /storage1/data11/TYMEFLIES_phage/*/VIBRANT_*.a/VIBRANT_results_*.a/VIBRANT_annotations_*.a.tsv |";
while (<IN>){
	chomp;
	my $file = $_;
	open INN, $file;
	while (<INN>){
		chomp;
		if (/^protein/){
			$annotation_header = $_;
		}else{
			my $line = $_;
			my @tmp = split (/\t/, $line);
			my $gene_short = $tmp[0];
			
			# Get the $gene_short_wo_fragment_x
			my $gene_short_wo_fragment_x = $gene_short;
			if ($gene_short_wo_fragment_x =~ /fragment\_\d\_/){
				$gene_short_wo_fragment_x =~ s/fragment\_\d\_//g;
			}			
			
			my $annotation = $line;
			if (exists $All_gene_seq_ID_map{$gene_short_wo_fragment_x}){
				my $gene = $All_gene_seq_ID_map{$gene_short_wo_fragment_x};
				$All_gene_annotation{$gene} = $annotation;
			}
		}
	}
	close INN;
}
close IN;

# Step 3 Store all pNpS result for viral species containing four AMGs
## Step 3.1 Store the viral species containing four AMGs
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

## Step 3.2 Store all pNpS result for all genes
my %Gene2year2pNpS = (); # $gene => $year => $pNpS
my %Year = (); # $year => 1
open IN, "MetaPop.for_each_year/MetaPop/10.Microdiversity/global_gene_microdiversity.tsv";
while (<IN>){
	chomp;
	if (!/^contig/){
		my @tmp = split (/\t/);
		my $gene = $tmp[0];
		my $snp_represent = $tmp[-1];
		my $pNpS = $tmp[-2];
		if ($pNpS eq ""){
			$pNpS = "NA";
		}
		my ($year) = $tmp[1] =~ /^(.+?)\./;
		$Year{$year} = 1;
		
		if ($snp_represent eq "TRUE"){
			$Gene2year2pNpS{$gene}{$year} = $pNpS;
		}
	}
}
close IN;

foreach my $gene (sort keys %All_gene_seq_ID){
	foreach my $year (sort keys %Year){
		if (!exists $Gene2year2pNpS{$gene}{$year}){
			$Gene2year2pNpS{$gene}{$year} = "NA";
		}
	}
}

# Step 4 Write down results
## Get the Viral gn to gene number hash
my %Viral_gn2gene_num = (); # $viral_gn => $gene_num
my %Viral_gn2genes = (); # $viral_gn => $genes (collection of $gene, separated by "\t")
foreach my $gene (sort keys %All_gene_seq_ID){
	my ($viral_gn) = $gene =~ /^(.+?\_\_.+?)\_\_.+?$/;
	$Viral_gn2gene_num{$viral_gn}++;
	if (!exists $Viral_gn2genes{$viral_gn}){
		$Viral_gn2genes{$viral_gn} = $gene;
	}else{
		$Viral_gn2genes{$viral_gn} .= "\t".$gene;
	}
}

## Step 4.1 Write down gene pN/pS results for 2000-2003 and 2016-2019 in viral species containing four AMGs
open OUT, ">MetaPop.for_each_year/MetaPop/pNpS_result.for_2000-2003_and_2016-2019.for_four_AMGs.txt";
print OUT "Gene\tViral gn\tViral gn characteristics\t2000\t2001\t2002\t2003\t2016\t2017\t2018\t2019\t$annotation_header\n";
foreach my $viral_gn (sort keys %Viral_species_containing_four_AMGs){
	my @Genes = split (/\t/, $Viral_gn2genes{$viral_gn});
	foreach my $gene (@Genes){
		my @Year_2000_2003_and_2016_2019 = qw/2000 2001 2002 2003 2016 2017 2018 2019/;
		my @pNpS_values_for_year_2000_2003_and_2016_2019 = ();
		
		foreach my $year (@Year_2000_2003_and_2016_2019){
			push @pNpS_values_for_year_2000_2003_and_2016_2019, $Gene2year2pNpS{$gene}{$year};
		}
		
		my $pNpS_values_for_year_2000_2003_and_2016_2019 = join("\t", @pNpS_values_for_year_2000_2003_and_2016_2019);
		
		my $gene_annotation = $All_gene_annotation{$gene};
		if (! $gene_annotation){
			print "$gene annotation is not present\n";
			$gene_annotation = "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
		}
		print OUT "$gene\t$viral_gn\t$Viral_species_containing_four_AMGs{$viral_gn}\t$pNpS_values_for_year_2000_2003_and_2016_2019\t$gene_annotation\n";
	}
}
close OUT;

## Step 4.2 Write down gene pN/pS results for 2000-2003 and 2016-2019 in viral species
open OUT, ">MetaPop.for_each_year/MetaPop/pNpS_result.for_2000-2003_and_2016-2019.txt";
print OUT "Gene\tViral gn\t2000\t2001\t2002\t2003\t2016\t2017\t2018\t2019\t$annotation_header\n";
foreach my $viral_gn (sort keys %Viral_gn2genes){
	my @Genes = split (/\t/, $Viral_gn2genes{$viral_gn});
	foreach my $gene (@Genes){
		my @Year_2000_2003_and_2016_2019 = qw/2000 2001 2002 2003 2016 2017 2018 2019/;
		my @pNpS_values_for_year_2000_2003_and_2016_2019 = ();
		
		foreach my $year (@Year_2000_2003_and_2016_2019){
			push @pNpS_values_for_year_2000_2003_and_2016_2019, $Gene2year2pNpS{$gene}{$year};
		}
		
		my $pNpS_values_for_year_2000_2003_and_2016_2019 = join("\t", @pNpS_values_for_year_2000_2003_and_2016_2019);
		
		my $gene_annotation = $All_gene_annotation{$gene};
		if (! $gene_annotation){
			print "$gene annotation is not present\n";
			$gene_annotation = "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
		}
		print OUT "$gene\t$viral_gn\t$pNpS_values_for_year_2000_2003_and_2016_2019\t$gene_annotation\n";
	}
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
