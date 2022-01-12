#!/usr/bin/perl

use strict;
use warnings;

###################################################################
# Aim: Add metagenome info and KEGG pathway info to the AMG table #
###################################################################

# Store all metagenomes info
my %Meta_info = (); # $img_id => $date_and_season (for example, "2000-03-15	| 75 | Spring")
open IN, "TYMEFLIES_metagenome_info.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		my @tmp = split (/\t/);
		my $img_id = $tmp[0];
		my $date_and_season = "$tmp[8] \| $tmp[9] \| $tmp[10]";
		$Meta_info{$img_id} = $date_and_season;
	}
}
close IN;

# Store KEGG pathway info
my %KEGG_pathway_info = (); # $map => [0] $metabolism [1] $pathway [2] $kos
open IN, "VIBRANT_KEGG_pathways_summary.tsv";
while (<IN>){
	chomp;
	if (!/^Entry/){
		my @tmp = split (/\t/);
		my $map = $tmp[0];
		my $metabolism = $tmp[1];
		my $pathway = $tmp[2];
		my $kos = $tmp[3];
		$KEGG_pathway_info{$map}[0] = $metabolism;
		$KEGG_pathway_info{$map}[1] = $pathway;
		$KEGG_pathway_info{$map}[2] = $kos;
	}
}
close IN;

# Store AMG summary table
my %AMG_summary = (); # $pro => [0] $amg_ko [1] $amg_ko_name [2] $pfam [3] $pfam_name
open IN, "AMG_analysis/AMG_summary.txt";
while (<IN>){
	chomp;
	if (/^protein/){
		#$header = $_;
	}else{
		my @tmp = split(/\t/);
		my $pro = $tmp[0];
		my $amg_ko = $tmp[1];
		my $amg_ko_name = $tmp[2];
		if ($amg_ko_name =~ /^\".+\"$/){ # Delete the quote marks in the front and in the end
			$amg_ko_name =~ s/^\"//g;
			$amg_ko_name =~ s/\"$//g;
		}
		my $pfam = $tmp[3];
		my $pfam_name = $tmp[4];
		
		$AMG_summary{$pro}[0] = $amg_ko;
		$AMG_summary{$pro}[1] = $amg_ko_name;
		$AMG_summary{$pro}[2] = $pfam;
		$AMG_summary{$pro}[3] = $pfam_name;
	}	
}
close IN;

# Add metagenome info and KEGG pathway info to the AMG summary hash
foreach my $pro (sort keys %AMG_summary){
	my $amg_ko = $AMG_summary{$pro}[0];
	my $amg_ko_name = $AMG_summary{$pro}[1];
	my $pfam = $AMG_summary{$pro}[2];
	my $pfam_name = $AMG_summary{$pro}[3];
	
	my ($img_id) = $pro =~ /^(.+?)\_\_/;
	my $date_and_season = $Meta_info{$img_id};
	
	my $metabolisms = "";
	my $pathways = "";
	
	foreach my $map (sort keys %KEGG_pathway_info){
		my $kos = $KEGG_pathway_info{$map}[2];
		if ($kos =~ /$amg_ko/){
			my $metabolism = $KEGG_pathway_info{$map}[0];
			my $pathway = $KEGG_pathway_info{$map}[1];
			
			if (!$metabolisms){
				$metabolisms = $metabolism;
			}else{
				$metabolisms .= " \| ".$metabolism; 
			}
			
			if (!$pathways){
				$pathways = $pathway;
			}else{
				$pathways .= " \| ".$pathway; 
			}			
		}
	}
	
	$AMG_summary{$pro}[4] = $date_and_season;
	$AMG_summary{$pro}[5] = $metabolisms;
	$AMG_summary{$pro}[6] = $pathways;
}

# Write down new AMG_summary.txt
open OUT, ">AMG_analysis/AMG_summary_new.txt";
print OUT "Protein\tMetagenome date and season ('date | date in the year | season')\tAMG KO\tAMG KO name\tPfam\tPfam name\tMetabolisms (mutiple metabolisms separated by '|')\tPathways (mutiple pathways separated by '|')\n";
foreach my $pro (sort keys %AMG_summary){
	my $date_and_season = $AMG_summary{$pro}[4];
	my $amg_ko = $AMG_summary{$pro}[0];
	my $amg_ko_name = $AMG_summary{$pro}[1];
	my $pfam = $AMG_summary{$pro}[2];
	my $pfam_name = $AMG_summary{$pro}[3];
	my $metabolisms = $AMG_summary{$pro}[5];
	my $pathways = $AMG_summary{$pro}[6];
	print OUT "$pro\t$date_and_season\t$amg_ko\t$amg_ko_name\t$pfam\t$pfam_name\t$metabolisms\t$pathways\n";
}
close OUT;

# Sort the result by date
`cat AMG_analysis/AMG_summary_new.txt | sort -k 2 -n > tmp`;
`mv tmp AMG_analysis/AMG_summary_new.txt`;

# Replace the old summary file
`rm AMG_analysis/AMG_summary.txt`;
`mv AMG_analysis/AMG_summary_new.txt AMG_analysis/AMG_summary.txt`;


