#!/usr/bin/perl

use strict;
use warnings;

# AIM: Get phage scaffold and vMAGs length and completeness info for plotting
# Including the following groups:
# 1) Phage scaffolds: for length plotting n = 1,820,639; for completeness plotting: 1,618,826
# 2) vMAGs + unbinned scaffolds: for length plotting n = 1,307,400; for completeness plotting: 1,182,522
# 3) vMAGs: for length plotting n = 222,390; for completeness plotting: 218,547
# 4) Binned scaffolds (within vMAGs): for length plotting n = 735,629 ; for completeness plotting: 654,851
# 5) Unbinned scaffolds: for length plotting n = 1,085,010; for completeness plotting: 963,975

# Hash for Phage scaffolds: %Scf2length_n_completeness; [0] for length; [1] for completeness; Should be 1,820,639 lines
# Hash for vMAGs + unbinned scaffolds: %vMAG_and_unbinned_scaffold2length_n_completeness; [0] for length; [1] for completeness
# Hash for vMAGs: %vMAG2length_n_completeness; [0] for length; [1] for completeness
# Hash for Binned scaffolds (within vMAGs): %Binned_scaffold2length_n_completeness; [0] for length; [1] for completeness
# Hash for Unbinned scaffolds: %Unbinned_scaffold2length_n_completeness; [0] for length; [1] for completeness

# Step 1. Get all scaffolds length and completeness information
my %Scf2length_n_completeness = (); # $scf => [0] $length [1] $completeness
open IN, "ls /storage1/data11/TYMEFLIES_phage/33*/CheckV_phage_scaffold/quality_summary.tsv |";
while (<IN>){
	chomp;
	my $file = $_;
	open INN, "$file";
	while (<INN>){
		chomp;
		if (!/^contig_id/){
			my @tmp = split (/\t/);
			my $scf = $tmp[0];
			my $length = $tmp[1];
			my $completeness = $tmp[9];
			$Scf2length_n_completeness{$scf}[0] = $length;
			$Scf2length_n_completeness{$scf}[1] = $completeness;
		}
	}
	close INN;
}
close IN;

# Step 2. Get scf 2 bin information
my %Scf_full2bin = (); # $scf_full => $bin
my %Bin = (); # $bin => $scf_full collection, separated by "\t"
my %Scf_full2scf = (); # $scf_full => $scf
my %Scf2scf_full = (); # $scf => $scf_full

open IN, "find /storage1/data11/TYMEFLIES_phage/33* -maxdepth 1 -name 'vRhyme_best_bins_fasta_parsed' -type d |";
while (<IN>){
	chomp;
	my $folder_adr = $_;
	open INN, "grep '>' $folder_adr/*.fasta |";
	while (<INN>){
		chomp; # for example, a line looks like: 
	# /storage1/data11/TYMEFLIES_phage/3300044855/vRhyme_best_bins_fasta_parsed/3300044855__vRhyme_100.fasta:>3300044855__vRhyme_100__Ga0453711_005625
		my ($bin) = $_ =~ /vRhyme_best_bins_fasta_parsed\/(.+?)\.fasta\:\>/;
		my ($scf_full) = $_ =~ /\:\>(.+?)$/;
		my ($scf) = $scf_full =~ /^.+?\_\_.+?\_\_(.+?)$/;
		
		$Scf_full2bin{$scf_full} = $bin;
		$Scf_full2scf{$scf_full} = $scf;
		
		if (!exists $Bin{$bin}){
			$Bin{$bin} = $scf_full;
		}else{
			$Bin{$bin} .= "\t".$scf_full;
		}
		
		$Scf2scf_full{$scf} = $scf_full;
	}
	close INN;
}
close IN;

# Step 3. Test the code to see if scaffold number and bin number are both correct
my $num_Bin = 0;
my $num_Scf_full2scf = 0;
my $num_Scf = 0;

$num_Bin = scalar (keys %Bin);
$num_Scf_full2scf = scalar (keys %Scf_full2scf);
$num_Scf = scalar (keys %Scf2length_n_completeness);

print "The number of bins is: $num_Bin\n"; # Should be 1,307,400
print "The number of scaffolds (full) is: $num_Scf_full2scf\n"; # Should be 1,820,639
print "The number of scaffolds is: $num_Scf\n"; # Should be 1,820,639


# Step 4. Get bin completeness hash
my %Bin2completeness = (); # $bin => $completeness
open IN, "/storage1/data11/TYMEFLIES_phage/CheckV_phage_bin_all/quality_summary.tsv";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $bin = $tmp[0];
	my $completeness = $tmp[9];
	$Bin2completeness{$bin} = $completeness;
}
close IN;

# Step 5. Get vMAG length and completeness hash and unbinned scaffolds length and completeness
my %vMAG2length_n_completeness = (); # $vMAG (vRhyme bin) => [0] $length [1] $completeness
my %Unbinned_scaffold2length_n_completeness = (); # $unbinned_scaffold (vRhyme unbinned) => [0] $length [1] $completeness
foreach my $bin (sort keys %Bin){
	if ($bin !~ /unbinned/){ # If the bin name does not contain "unbinned" (it is a vMAG)
		my $vMAG = $bin;
		my $length = 0;
		my $completeness = "";
		
		my @Scf_full_collection = split (/\t/, $Bin{$bin});
		foreach my $scf_full (@Scf_full_collection){
			my $scf = $Scf_full2scf{$scf_full};
			my $scf_length = $Scf2length_n_completeness{$scf}[0];
			$length += $scf_length;
		}
		
		$completeness = $Bin2completeness{$bin};
		
		$vMAG2length_n_completeness{$vMAG}[0] = $length;
		$vMAG2length_n_completeness{$vMAG}[1] = $completeness;
	}else{ # The unbinned scaffold
		my $unbinned_scaffold = $bin;
		my $length = 0;
		my $completeness = "";
		
		my $scf_full = $Bin{$bin}; # Just one scaffold in the bin
		my $scf = $Scf_full2scf{$scf_full};
		$length = $Scf2length_n_completeness{$scf}[0];
		
		$completeness = $Bin2completeness{$bin};
		
		$Unbinned_scaffold2length_n_completeness{$unbinned_scaffold}[0] = $length;
		$Unbinned_scaffold2length_n_completeness{$unbinned_scaffold}[1] = $completeness;	
	}
}

my %vMAG_and_unbinned_scaffold2length_n_completeness = (); # $bin (contains both vMAG bin and unbinned scaffold) => [0] $length [1] $completeness
foreach my $vMAG (sort keys %vMAG2length_n_completeness){
	my $length = $vMAG2length_n_completeness{$vMAG}[0];
	my $completeness = $vMAG2length_n_completeness{$vMAG}[1];
	
	$vMAG_and_unbinned_scaffold2length_n_completeness{$vMAG}[0] = $length;
	$vMAG_and_unbinned_scaffold2length_n_completeness{$vMAG}[1] = $completeness;
}

foreach my $unbinned_scaffold (sort keys %Unbinned_scaffold2length_n_completeness){
	my $length = $Unbinned_scaffold2length_n_completeness{$unbinned_scaffold}[0];
	my $completeness = $Unbinned_scaffold2length_n_completeness{$unbinned_scaffold}[1];
	
	$vMAG_and_unbinned_scaffold2length_n_completeness{$unbinned_scaffold}[0] = $length;
	$vMAG_and_unbinned_scaffold2length_n_completeness{$unbinned_scaffold}[1] = $completeness;
}

# Step 6. Get the Binned scaffolds length and completeness hash
my %Binned_scaffold2length_n_completeness = (); # $binned_scaffold => [0] $length [1] $completeness
foreach my $vMAG (sort keys %vMAG2length_n_completeness){
	my @Scf_full_collection = split (/\t/, $Bin{$vMAG});
	
	foreach my $scf_full (@Scf_full_collection){
		my $scf = $Scf_full2scf{$scf_full};
		
		my $length = $Scf2length_n_completeness{$scf}[0];
		my $completeness = $Scf2length_n_completeness{$scf}[1];
		
		$Binned_scaffold2length_n_completeness{$scf_full}[0] = $length;
		$Binned_scaffold2length_n_completeness{$scf_full}[1] = $completeness;
	}
}

# Step 7. Write down all groups
open OUT, ">phage_scaffold_and_vMAGs_length_and_completeness_info.txt";
# Write down the hash for Phage scaffolds: %Scf2length_n_completeness; [0] for length; [1] for completeness
foreach my $scf (sort keys %Scf2length_n_completeness){
	print OUT "$scf\tPhage scaffolds\t$Scf2length_n_completeness{$scf}[0]\t$Scf2length_n_completeness{$scf}[1]\n";
}
# Write down the hash for vMAGs + unbinned scaffolds: %vMAG_and_unbinned_scaffold2length_n_completeness; [0] for length; [1] for completeness
foreach my $key (sort keys %vMAG_and_unbinned_scaffold2length_n_completeness){
	print OUT "$key\tvMAGs + unbinned scaffolds\t$vMAG_and_unbinned_scaffold2length_n_completeness{$key}[0]\t$vMAG_and_unbinned_scaffold2length_n_completeness{$key}[1]\n";
}
# Write down the hash for vMAGs: %vMAG2length_n_completeness; [0] for length; [1] for completeness
foreach my $key (sort keys %vMAG2length_n_completeness){
	print OUT "$key\tvMAGs\t$vMAG2length_n_completeness{$key}[0]\t$vMAG2length_n_completeness{$key}[1]\n";
}
# Write down the hash for Binned scaffolds (within vMAGs): %Binned_scaffold2length_n_completeness; [0] for length; [1] for completeness
foreach my $key (sort keys %Binned_scaffold2length_n_completeness){
	print OUT "$key\tBinned scaffolds (within vMAGs)\t$Binned_scaffold2length_n_completeness{$key}[0]\t$Binned_scaffold2length_n_completeness{$key}[1]\n";
}
# Write down the hash for Unbinned scaffolds: %Unbinned_scaffold2length_n_completeness; [0] for length; [1] for completeness
foreach my $key (sort keys %Unbinned_scaffold2length_n_completeness){
	print OUT "$key\tUnbinned scaffolds\t$Unbinned_scaffold2length_n_completeness{$key}[0]\t$Unbinned_scaffold2length_n_completeness{$key}[1]\n";
}
close OUT;





