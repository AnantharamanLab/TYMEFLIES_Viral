#!/usr/bin/perl

use strict;
use warnings;

################################################
# Aim: generate AMG table for all metagenomes  #
################################################

# Store all metagenomes
my %IMGID = (); # $img_id => 1
open IN, "ls -l | sed -r s'/ +/\t/g' | grep ^d | cut -f 9 | grep '^33' |";
while (<IN>){
	chomp;
	my $img_id = $_;
	$IMGID{$img_id} = 1;
}
close IN;

my %AMG_summary = (); # $pro_full => [0] $amg_ko, [1] $amg_ko_name, [2] $pfam, [3] $pfam_name
foreach my $img_id (sort keys %IMGID){
	
	# Store AMG annotations from each metagenome's VIBRANT result
	my %AMG_summary_for_each_metagenome = (); # $pro (the original protein name) => [0] $amg_ko, [1] $amg_ko_name, [2] $pfam, [3] $pfam_name
	open IN, "$img_id/VIBRANT_$img_id.a/VIBRANT_results_$img_id.a/VIBRANT_AMG_individuals_$img_id.a.tsv";
	while (<IN>){
		chomp;
		if (!/^protein/){
			my @tmp = split (/\t/);
			my $pro = $tmp[0];
			my $amg_ko = "";
			my $amg_ko_name = "";
			my $pfam = "n/a";
			my $pfam_name = "n/a";
			$amg_ko = $tmp[2];
			$amg_ko_name = $tmp[3];
			if ($tmp[4]){
				$pfam = $tmp[4];
			}
			if ($tmp[5]){			
				$pfam_name = $tmp[5];	
			}			
			$AMG_summary_for_each_metagenome{$pro}[0] = $amg_ko;
			$AMG_summary_for_each_metagenome{$pro}[1] = $amg_ko_name;
			$AMG_summary_for_each_metagenome{$pro}[2] = $pfam;
			$AMG_summary_for_each_metagenome{$pro}[3] = $pfam_name;
		}
	}
	close IN;
	
	# Store all sequence head from all phage genomes	
	my %Seq_head_for_each_metagenome = ();	# $pro_full (for example here, 3300020498__vRhyme_unbinned237__Ga0208050_1001025_1) => 1
	open IN, "ls $img_id/vRhyme_best_bins_fasta_parsed/*vRhyme_*.faa | ";
	while (<IN>){
		chomp;
		my $file = $_;
		my %Hash_tmp = _store_seq_head($file); 
		%Seq_head_for_each_metagenome = (%Seq_head_for_each_metagenome, %Hash_tmp);
	}
	close IN;
	
	foreach my $pro_full (sort keys %Seq_head_for_each_metagenome){
		my ($pro_original) = $pro_full =~ /^33.+?\_\_vRhyme.+?\_\_(.+?)$/;
		if (exists $AMG_summary_for_each_metagenome{$pro_original}){
			$AMG_summary{$pro_full}[0] = $AMG_summary_for_each_metagenome{$pro_original}[0];
			$AMG_summary{$pro_full}[1] = $AMG_summary_for_each_metagenome{$pro_original}[1];
			$AMG_summary{$pro_full}[2] = $AMG_summary_for_each_metagenome{$pro_original}[2];
			$AMG_summary{$pro_full}[3] = $AMG_summary_for_each_metagenome{$pro_original}[3];	
		}else{
			foreach my $pro (sort keys %AMG_summary_for_each_metagenome){
				if ($pro =~ /fragment/){ # If this protein is from prophage, contains "fragment" inside; after prophage combining step, this protein will not contain "fragment" in its name
					my ($pro_wo_fragment_1) = $pro =~ /^(.+?)\_fragment/;
					my ($pro_wo_fragment_2) = $pro =~ /fragment\_.+?(\_\d+?)$/;
					my $pro_wo_fragment = $pro_wo_fragment_1.$pro_wo_fragment_2; # The protein name without fragment_X inside
					if ($pro_wo_fragment eq $pro_original){
						$AMG_summary{$pro_full}[0] = $AMG_summary_for_each_metagenome{$pro}[0];
						$AMG_summary{$pro_full}[1] = $AMG_summary_for_each_metagenome{$pro}[1];
						$AMG_summary{$pro_full}[2] = $AMG_summary_for_each_metagenome{$pro}[2];
						$AMG_summary{$pro_full}[3] = $AMG_summary_for_each_metagenome{$pro}[3];						
					}
				}
			}
		}
	}
}

# Store the my %AMG_summary = (); # $pro_full => [0] $amg_ko, [1] $amg_ko_name, [2] $pfam, [3] $pfam_name
`mkdir AMG_analysis`;
open OUT, ">AMG_analysis/AMG_summary.txt";
print OUT "protein\tAMG KO\tAMG KO name\tPfam\tPfam name\n";
foreach my $pro_full (sort keys %AMG_summary){
	print OUT "$pro_full\t$AMG_summary{$pro_full}[0]\t$AMG_summary{$pro_full}[1]\t$AMG_summary{$pro_full}[2]\t$AMG_summary{$pro_full}[3]\n";
}
close OUT;

# Subroutine

# Store sequence head (clean) of input 
sub _store_seq_head{ 
	my $file = $_;
	my %Seq_head = ();
	open _IN, "$file";
	while (<_IN>){
		chomp;
		if (/^>/){
			my $head = $_;
			$head =~ s/^>//g;
			$Seq_head{$head} = 1;
		}
	}
	close _IN;
	return %Seq_head;
}	


