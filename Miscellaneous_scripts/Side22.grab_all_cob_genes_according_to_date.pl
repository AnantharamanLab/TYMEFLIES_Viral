#!/usr/bin/perl

use strict;
use warnings;

# Aim: Grab all the Cob proteins and put them into individual files according to date

# Step 1 Grab all cob proteins
my %Cob2KO = ('cobS'=>'K09882', 'cobT'=>'K09883', 'cobL'=>'K00595', 'cobA'=>'K19221', 'cobC'=>'K22316');
my %CobS_seq = (); # Store all the proteins belonging to CobS
my %CobT_seq = (); # Store all the proteins belonging to CobT
my %CobL_seq = (); # Store all the proteins belonging to CobL
my %CobA_seq = (); # Store all the proteins belonging to CobA
my %CobC_seq = (); # Store all the proteins belonging to CobC

## Step 1.1 Store AMG summary table 
my %AMG_summary = (); # $pro => $ko
my %KOs= (); # $ko => 1;
my %IMG2date = (); # $img_id => $date_n_season
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
	}
}
close IN;

## Step 1.2 Store all the viral protein sequences 
my %All_viral_seq = (); 
%All_viral_seq = _store_seq("/storage1/data11/TYMEFLIES_phage/Cluster_phage_genomes/All_phage_genome.faa");

## Step 1.3 Make a hash of Cob protein to Cob
my %Cob_pro2cob = ();
foreach my $pro (sort keys %AMG_summary){
	my $ko = $AMG_summary{$pro};
	foreach my $cob (sort keys %Cob2KO){
		my $ko_from_cob = $Cob2KO{$cob};
		if ($ko eq $ko_from_cob){
			$Cob_pro2cob{">$pro"} = $cob;
		}
	}
}

## Step 1.4 Store only Cob proteins
foreach my $cob_pro (sort keys %Cob_pro2cob){
	my $cob = $Cob_pro2cob{$cob_pro};
	if ($cob eq "cobS"){
		$CobS_seq{$cob_pro} = $All_viral_seq{$cob_pro};
	}elsif ($cob eq "cobT"){
		$CobT_seq{$cob_pro} = $All_viral_seq{$cob_pro};
	}elsif ($cob eq "cobL"){
		$CobL_seq{$cob_pro} = $All_viral_seq{$cob_pro};
	}elsif ($cob eq "cobA"){
		$CobA_seq{$cob_pro} = $All_viral_seq{$cob_pro};
	}elsif ($cob eq "cobC"){
		$CobC_seq{$cob_pro} = $All_viral_seq{$cob_pro};
	}
}

# Step 2 Separate and write down Cob proteins according to the date
`mkdir Cob_proteins`;
my %CobS_IMG = (); # Store all the IMG ID in CobS proteins; $img_id => 1
foreach my $pro (sort keys %CobS_seq){
	my ($img_id) = $pro =~ /^>(33.+?)\_\_/;
	$CobS_IMG{$img_id} = 1;
}

foreach my $img_id (sort keys %CobS_IMG){
	my $date_n_season = $IMG2date{$img_id};
	my ($date) = $date_n_season =~ /^(.+?)\s/;
	open OUT, ">Cob_proteins/CobS.$img_id.$date.faa";
	foreach my $pro (sort keys %CobS_seq){
		if ($pro =~ />$img_id\_\_/){
			print OUT "$pro\n$CobS_seq{$pro}\n";
		}
	}
	close OUT;
}

my %CobT_IMG = (); # Store all the IMG ID in CobT proteins; $img_id => 1
foreach my $pro (sort keys %CobT_seq){
	my ($img_id) = $pro =~ /^>(33.+?)\_\_/;
	$CobT_IMG{$img_id} = 1;
}

foreach my $img_id (sort keys %CobT_IMG){
	my $date_n_season = $IMG2date{$img_id};
	my ($date) = $date_n_season =~ /^(.+?)\s/;
	open OUT, ">Cob_proteins/CobT.$img_id.$date.faa";
	foreach my $pro (sort keys %CobT_seq){
		if ($pro =~ />$img_id\_\_/){
			print OUT "$pro\n$CobT_seq{$pro}\n";
		}
	}
	close OUT;
}

my %CobL_IMG = (); # Store all the IMG ID in CobL proteins; $img_id => 1
foreach my $pro (sort keys %CobL_seq){
	my ($img_id) = $pro =~ /^>(33.+?)\_\_/;
	$CobL_IMG{$img_id} = 1;
}

foreach my $img_id (sort keys %CobL_IMG){
	my $date_n_season = $IMG2date{$img_id};
	my ($date) = $date_n_season =~ /^(.+?)\s/;
	open OUT, ">Cob_proteins/CobL.$img_id.$date.faa";
	foreach my $pro (sort keys %CobL_seq){
		if ($pro =~ />$img_id\_\_/){
			print OUT "$pro\n$CobL_seq{$pro}\n";
		}
	}
	close OUT;
}

my %CobA_IMG = (); # Store all the IMG ID in CobA proteins; $img_id => 1
foreach my $pro (sort keys %CobA_seq){
	my ($img_id) = $pro =~ /^>(33.+?)\_\_/;
	$CobA_IMG{$img_id} = 1;
}

foreach my $img_id (sort keys %CobA_IMG){
	my $date_n_season = $IMG2date{$img_id};
	my ($date) = $date_n_season =~ /^(.+?)\s/;
	open OUT, ">Cob_proteins/CobA.$img_id.$date.faa";
	foreach my $pro (sort keys %CobA_seq){
		if ($pro =~ />$img_id\_\_/){
			print OUT "$pro\n$CobA_seq{$pro}\n";
		}
	}
	close OUT;
}

my %CobC_IMG = (); # Store all the IMG ID in CobC proteins; $img_id => 1
foreach my $pro (sort keys %CobC_seq){
	my ($img_id) = $pro =~ /^>(33.+?)\_\_/;
	$CobC_IMG{$img_id} = 1;
}

foreach my $img_id (sort keys %CobC_IMG){
	my $date_n_season = $IMG2date{$img_id};
	my ($date) = $date_n_season =~ /^(.+?)\s/;
	open OUT, ">Cob_proteins/CobC.$img_id.$date.faa";
	foreach my $pro (sort keys %CobC_seq){
		if ($pro =~ />$img_id\_\_/){
			print OUT "$pro\n$CobC_seq{$pro}\n";
		}
	}
	close OUT;
}



## Subroutine

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

