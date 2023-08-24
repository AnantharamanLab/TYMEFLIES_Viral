#!/usr/bin/perl
  
use strict;
use warnings;

# Aim: Get gene frequency distribution for selected viral genomes

# Step 1 Get gene gain and loss information for selected viral genomes
my @Viral_gn_selected = qw/3300043571__vRhyme_51 3300042549__vRhyme_1492/;
my %Viral_gn2gene_gain_collection = (); # $viral_gn => $collection of gene gain
my %Viral_gn2gene_loss_collection = (); # $viral_gn => $collection of gene loss
open IN, "/storage1/data11/TYMEFLIES_phage/MetaPop.for_each_year/MetaPop/Viral_gene2gain_or_loss.for_four_AMGs.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $viral_gn = $tmp[0];
	my $gene_gain = $tmp[3];
	my $gene_loss = $tmp[5];
	foreach my $key (@Viral_gn_selected){
		if ($key eq $viral_gn){
			$Viral_gn2gene_gain_collection{$viral_gn} = $gene_gain;
			$Viral_gn2gene_loss_collection{$viral_gn} = $gene_loss;
		}			
	}
}
close IN;

my %Gene_loss_or_gain = (); # $gene => "gain" or "loss"
foreach my $viral_gn (sort keys %Viral_gn2gene_gain_collection){
	my @Gene_gain_collection = split(/\,/,$Viral_gn2gene_gain_collection{$viral_gn});
	foreach my $key (@Gene_gain_collection){
		$Gene_loss_or_gain{$key} = "gain";
	}
	my @Gene_loss_collection = split(/\,/,$Viral_gn2gene_loss_collection{$viral_gn});
	foreach my $key (@Gene_loss_collection){
		$Gene_loss_or_gain{$key} = "loss";
	}	
}

# Step 2 Get gene frequency distribution for selected viral genomes
## Step 2.1 Store Viral_gene2year2gene_freq.for_four_AMGs.txt hash
my %Viral_gene2year2gene_freq_for_four_AMGs = (); # $viral_gene => $year => $gene_freq
my @Head = ();
my %Year = (); # $year => 1
open IN, "/storage1/data11/TYMEFLIES_phage/MetaPop.for_each_year/MetaPop/Viral_gene2year2gene_freq.for_four_AMGs.txt";
while (<IN>){
	chomp;
	if (/^Viral/){
		my @tmp = split(/\t/);
		@Head = @tmp;
	}else{
		my @tmp = split(/\t/);
		my $viral_gene = $tmp[0];
		for(my $i=1; $i<=$#tmp; $i++){
			my $year = $Head[$i];
			$Year{$year} = 1;
			my $gene_freq = $tmp[$i];
			$Viral_gene2year2gene_freq_for_four_AMGs{$viral_gene}{$year} = $gene_freq;
		}
	}
}
close IN;

## Step 2.2 Write down gene frequency for selected viral genomes
open OUT, ">MetaPop.for_each_year/MetaPop/Viral_gene2year2gene_freq.for_four_AMGs.for_selected_genomes.txt";
my $row=join("\t", sort keys %Year);
print OUT "Head\t$row\tGene gain or loss\n";
foreach my $tmp1 (sort keys %Gene_loss_or_gain){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %Year) {       
                if (exists $Viral_gene2year2gene_freq_for_four_AMGs{$tmp1}{$tmp2}){
                        push @tmp, $Viral_gene2year2gene_freq_for_four_AMGs{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"NA";
                }
        }
		push @tmp, $Gene_loss_or_gain{$tmp1};
        print OUT join("\t",@tmp)."\n";
}
close OUT;

