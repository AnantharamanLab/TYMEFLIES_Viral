#!/usr/bin/perl

use strict;
use warnings;
use File::Slurp;

my %IMGID = (); # $img_id => 1
open IN, "ls -l | sed -r s'/ +/\t/g' | grep ^d | cut -f 9 | grep '33' |";
while (<IN>){
	chomp;
	my $img_id = $_;
	$IMGID{$img_id} = 1;
}
close IN;

my %Result = (); # $img_id => [0] total scaffold num [1] scaffold over min2000 [2] phage num in total [3] lytic phage num [4] lysogenic phage num [5] prophage num [6] complete phage num
foreach my $img_id (sort keys %IMGID){
	if (-d "$img_id/VIBRANT_$img_id.a.v2.min2000"){
		my $total_scaf_num = 0; my $scaf_num_over_min2000 = 0; my $phage_num_in_total = 0; 
		my $lytic_phage_num = 0; my $lysogenic_phage_num = 0; my $prophage_num = 0; my $complete_phage_num = 0; # lysogenic phage here refers to lysogenic but not prophage
		
		my %Seq_total = _store_seq("$img_id/$img_id.a.fna"); 
		
		foreach my $key (sort keys %Seq_total){
			$total_scaf_num ++;
			if (length($Seq_total{$key}) >= 2000){
				$scaf_num_over_min2000++;
			}
		}
		
		open IN, "$img_id/VIBRANT_$img_id.a.v2.min2000/phage_results.min2000.txt";
		while (<IN>){
			chomp;
			my @tmp = split (/\t/);
			my $scf = $tmp[0];
			my $state = $tmp[1];
			my $scf_length = $tmp[2];
			
			$phage_num_in_total++;
			if ($state eq "lytic"){
				$lytic_phage_num++;
			}elsif($state eq "lysogenic_but_non-prophage"){
				$lysogenic_phage_num++;
			}elsif($state eq "prophage"){
				$prophage_num++;
			}
		}
		close IN;
		
		open IN, "$img_id/CheckV_phage_scaffold/quality_summary.tsv";
		while (<IN>){
			chomp;
			if (!/^contig/){
				my @tmp = split (/\t/);
				my $checkv_quality = $tmp[7];
				if ($checkv_quality eq "Complete"){
					$complete_phage_num++;
				}
			}
		}
		close IN;
		
		$Result{$img_id}[0] = $total_scaf_num;
		$Result{$img_id}[1] = $scaf_num_over_min2000;
		$Result{$img_id}[2] = $phage_num_in_total;
		$Result{$img_id}[3] = $lytic_phage_num;
		$Result{$img_id}[4] = $lysogenic_phage_num;
		$Result{$img_id}[5] = $prophage_num;
		$Result{$img_id}[6] = $complete_phage_num;
		
	}	
}

open OUT, ">VIBRANT_result_summary.txt";
print OUT "IMG ID\tTotal scaffold num\tScaffold num over min2000\tPhage num in total\tLytic phage num\tLysogenic phage (excluding prophage) num\tProphage num\tComplete phage num\n";
foreach my $img_id (sort keys %Result){
	print OUT "$img_id\t$Result{$img_id}[0]\t$Result{$img_id}[1]\t$Result{$img_id}[2]\t$Result{$img_id}[3]\t$Result{$img_id}[4]\t$Result{$img_id}[5]\t$Result{$img_id}[6]\n";
}
close OUT;



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
