#!/usr/bin/perl

use strict;
use warnings;

# Aim: Get the monthly viral tax (family level) and host tax (family level) distribution pattern

# Step 1 Store all viral abundance info
## Step 1.1 Store scaffold coverage for each phage genome
my %Scf2cov = (); # $scf => $cov (within individual metagenome); a typical $scf is like this: "Ga0453730_0000001_fragment_3"
open IN, "ls 33**/vRhyme_result/vRhyme_coverage_files/vRhyme_coverage_values.tsv | ";
while (<IN>){
	chomp;
	my $file = $_;
	my ($img_id) = $file =~ /(33\d+?)\//;
	open INN, "$file";
	while (<INN>){
		chomp;
		if (!/^scaffold/){
			my @tmp = split (/\t/);
			my $scf = $tmp[0];
			my $cov = $tmp[1];
			$Scf2cov{$scf} = $cov;
		}
	}
	close INN;
}
close IN;

## Step 1.2 Store the viral gn to viral scf hash
my %Viral_gn2viral_seq = (); # $viral_gn => $viral_seq collection separated by "\t"; a typical $viral_seq is like this: "3300033816__vRhyme_unbinned70__Ga0334980_0000585_fragment_1"
open IN, "/storage1/data11/TYMEFLIES_phage/Host_prediction/All_phage_genomes.fasta";
while (<IN>){
	chomp;
	if (/^>/){
		my $line = $_;
		my ($viral_seq) = $line =~ /^>(.+?)$/;
		my ($viral_gn) = $viral_seq =~ /^(.+?\_\_.+?)\_\_.+?$/;
		if (! exists $Viral_gn2viral_seq{$viral_gn}){
			$Viral_gn2viral_seq{$viral_gn} = $viral_seq;
		}else{
			$Viral_gn2viral_seq{$viral_gn} .= "\t".$viral_seq;
		}
	}
}
close IN;

## Step 1.3 Get AMG abundances (normalzied by 100M reads per metagenome)
my %IMG_ID2read_num = ();
open IN, "TYMEFLIES_metagenome_info.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		my @tmp = split (/\t/);
		my $img_id = $tmp[0];
		my $read_num = $tmp[13];
		$IMG_ID2read_num{$img_id} = $read_num;
	}
}	
close IN;

## Step 1.4 Get the Viral_gn2IMG2abun hash; abundance here is normalized by 100M reads/metagenome
my %Viral_gn2IMG2abun = (); # $viral_gn => $img_id => $abun
my %IMG2viral_gn = (); # $img_id => $viral_gns (collection of $viral_gn, separated by "\t")
foreach my $viral_gn (sort keys %Viral_gn2viral_seq){
	my @Viral_seqs = split (/\t/, $Viral_gn2viral_seq{$viral_gn}); 
	my @Scf_cov = (); # Store all the coverages of scaffolds within this genome
	foreach my $viral_seq (@Viral_seqs){
		my ($scf) = $viral_seq =~ /^.+?\_\_.+?\_\_(.+?)$/;
		my $cov = $Scf2cov{$scf};
		push @Scf_cov, $cov;
	}
	my $viral_gn_abun = _avg(@Scf_cov);
	
	my ($img_id) = $viral_gn =~ /^(.+?)\_\_/;
	my $read_num = $IMG_ID2read_num{$img_id};
	
	$viral_gn_abun = $viral_gn_abun / ($read_num / 100000000); # Normalize the abundance
	
	$Viral_gn2IMG2abun{$viral_gn}{$img_id} = $viral_gn_abun;
	
	if (!exists $IMG2viral_gn{$img_id}){
		$IMG2viral_gn{$img_id} = $viral_gn;
	}else{
		$IMG2viral_gn{$img_id} .= "\t".$viral_gn;
	}
}

# Step 2 Store the viral gn to tax hash
my %Viral_gn2tax = (); # $gn => $tax;
open IN, "/storage1/data11/TYMEFLIES_phage/Taxonomic_classification/Each_bin_tax_combined_result.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	$Viral_gn2tax{$tmp[0]} = $tmp[1];
}
close IN;

# Step 3 Store month to img_id hash
my %Month2num = (); # $month => $num_metagenome; Store how many samples (metagenomes) are there in each month for the whole datasets
my %Month2img_id = (); # $month => collection of $img_id separeted by "\t"
my %Year_month2num = (); # $year_month => $num_metagenome
my %Year_month2img_id = (); # $year_month (for example, "2010-06") => collection of $img_id separeted by "\t"
open IN, "TYMEFLIES_metagenome_info.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		my @tmp = split (/\t/);
		my $img_id = $tmp[0];
		my $date = $tmp[8];
		my ($month) = $date =~ /\d\d\d\d-(\d\d)-\d\d/; 
		$Month2num{$month}++;
		if (!exists $Month2img_id{$month}){
			$Month2img_id{$month} = $img_id;
		}else{
			$Month2img_id{$month} .= "\t".$img_id;
		}
		
		my ($year_month) = $date =~ /(\d\d\d\d-\d\d)-\d\d/;
		$Year_month2num{$year_month}++;
		if (!exists $Year_month2img_id{$year_month}){
			$Year_month2img_id{$year_month} = $img_id;
		}else{
			$Year_month2img_id{$year_month} .= "\t".$img_id;
		}		
	}
}
close IN;

# Step 4 Get the Month2family2abun hash
my %Month2family2abun = (); # $month => $family => $abun
my %Family = (); # $family => 1

foreach my $month (sort keys %Month2img_id){
	my @IMG_ID = split (/\t/, $Month2img_id{$month});
	my @Viral_gns = (); # All the viral genomes belonging to this $month
	foreach my $img_id (@IMG_ID){
		my @tmp = split (/\t/,$IMG2viral_gn{$img_id});
		foreach my $key (@tmp){
			push @Viral_gns, $key;
		}
	}
	
	foreach my $viral_gn (@Viral_gns){
		my $family = "Unclassified";
		if (exists $Viral_gn2tax{$viral_gn}){
			my $tax = $Viral_gn2tax{$viral_gn};
			my @Tax = split (/\;/, $tax);
			$family = "$Tax[4]\;$Tax[5]";
		}
		$Family{$family} = 1;
		
		my ($img_id) = $viral_gn =~ /^(.+?)\_/;
		$Month2family2abun{$month}{$family} += $Viral_gn2IMG2abun{$viral_gn}{$img_id};
	}
}

# Step 5 Write down Month2family2abun.txt
open OUT, ">Month2family2abun.txt";
my $row=join("\t", sort keys %Month2img_id);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %Family){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %Month2family2abun) {       
                if (exists $Month2family2abun{$tmp2}{$tmp1}){
                        push @tmp, $Month2family2abun{$tmp2}{$tmp1};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

## Step 6 Store the viral gn to host tax hash
my %Viral_gn2host_tax = (); # $gn => $tax;
open IN, "/storage1/data11/TYMEFLIES_phage/Host_prediction/Viral_gn2host_tax_final.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	$Viral_gn2host_tax{$tmp[0]} = $tmp[1];
}
close IN;

# Step 7 Get the Month2host_family2abun hash
my %Month2host_family2abun = (); # $month => $host_family => $abun
my %Host_family = (); # $host_family => 1

foreach my $month (sort keys %Month2img_id){
	my @IMG_ID = split (/\t/, $Month2img_id{$month});
	my @Viral_gns = (); # All the viral genomes belonging to this $month
	foreach my $img_id (@IMG_ID){
		my @tmp = split (/\t/,$IMG2viral_gn{$img_id});
		foreach my $key (@tmp){
			push @Viral_gns, $key;
		}
	}
	
	foreach my $viral_gn (@Viral_gns){
		my $host_family = "Unclassified";
		if (exists $Viral_gn2host_tax{$viral_gn}){
			my $host_tax = $Viral_gn2host_tax{$viral_gn};
			my @Host_tax = split (/\;/, $host_tax);
			$host_family = "$Host_tax[3]\;$Host_tax[4]";
		}
		$Host_family{$host_family} = 1;
		
		my ($img_id) = $viral_gn =~ /^(.+?)\_/;
		$Month2host_family2abun{$month}{$host_family} += $Viral_gn2IMG2abun{$viral_gn}{$img_id};
	}
}

# Step 8 Write down Month2host_family2abun.txt
open OUT, ">Month2host_family2abun.txt";
my $row2=join("\t", sort keys %Month2img_id);
print OUT "Head\t$row2\n";
foreach my $tmp1 (sort keys %Host_family){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %Month2host_family2abun) {       
                if (exists $Month2host_family2abun{$tmp2}{$tmp1}){
                        push @tmp, $Month2host_family2abun{$tmp2}{$tmp1};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;



# Subroutine 
sub _avg {
    my $total;
    $total += $_ foreach @_;
    # sum divided by number of components.
    return $total / @_;
}

