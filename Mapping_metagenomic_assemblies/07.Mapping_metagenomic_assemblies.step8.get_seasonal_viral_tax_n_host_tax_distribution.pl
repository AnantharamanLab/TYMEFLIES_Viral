#!/usr/bin/perl

use strict;
use warnings;

# Aim: Get the seasonal viral tax (family level) and host tax (family level) distribution pattern

# Step 1 Store all viral abundance info
my %Viral_gn2IMG2abun = (); # $viral_gn => $img_id => $abun
my %IMG2viral_gn = (); # $img_id => $viral_gns (collection of $viral_gn, separated by "\t")
open IN, "viral_gn2depth_normalized.txt";
while (<IN>){
	chomp;
	my @tmp = split(/\t/);
	my $viral_gn = $tmp[0];
	my $depth_normalized = $tmp[1];
	my ($img_id) = $viral_gn =~ /^(.+?)\_\_/;
	
	$Viral_gn2IMG2abun{$viral_gn}{$img_id} = $depth_normalized;
	
	if (!exists $IMG2viral_gn{$img_id}){
		$IMG2viral_gn{$img_id} = $viral_gn;
	}else{
		$IMG2viral_gn{$img_id} .= "\t".$viral_gn;
	}	
}
close IN;

# Step 2 Store the viral gn to tax hash
my %Viral_gn2tax = (); # $gn => $tax;
open IN, "/storage1/data11/TYMEFLIES_phage/Taxonomic_classification/Each_bin_tax_combined_result.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	$Viral_gn2tax{$tmp[0]} = $tmp[1];
}
close IN;

# Step 3 Store season to img_id hash
my %Season2num = (); # $season => $num_metagenome; Store how many samples (metagenomes) are there in each season for the whole datasets
my %Season2img_id = (); # $season => collection of $img_id separeted by "\t"
open IN, "TYMEFLIES_metagenome_info.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		my @tmp = split (/\t/);
		my $img_id = $tmp[0];
		my $date = $tmp[8];
		my $season = $tmp[10];
		$Season2num{$season}++;
		if (!exists $Season2img_id{$season}){
			$Season2img_id{$season} = $img_id;
		}else{
			$Season2img_id{$season} .= "\t".$img_id;
		}		
	}
}
close IN;

# Step 4 Get the Season2family2abun hash
my %Season2family2abun = (); # $season => $family => $abun
my %Family = (); # $family => 1

foreach my $season (sort keys %Season2img_id){
	my @IMG_ID = split (/\t/, $Season2img_id{$season});
	my @Viral_gns = (); # All the viral genomes belonging to this $season
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
			$family = "$Tax[3]\;$Tax[4]\;$Tax[5]";
		}
		$Family{$family} = 1;
		
		my ($img_id) = $viral_gn =~ /^(.+?)\_/;
		$Season2family2abun{$season}{$family} += $Viral_gn2IMG2abun{$viral_gn}{$img_id};
	}
}

my @Season = ("Spring", "Clearwater", "Early Summer", "Late Summer", "Fall", "Ice-on");
# Step 5 Write down Season2family2abun.txt
open OUT, ">Season2family2abun.txt";
my $row=join("\t", @Season);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %Family){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %Season2family2abun) {       
                if (exists $Season2family2abun{$tmp2}{$tmp1}){
                        push @tmp, $Season2family2abun{$tmp2}{$tmp1};
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

# Step 7 Get the Season2host_family2abun hash
my %Season2host_family2abun = (); # $season => $host_family => $abun
my %Host_family = (); # $host_family => 1

foreach my $season (sort keys %Season2img_id){
	my @IMG_ID = split (/\t/, $Season2img_id{$season});
	my @Viral_gns = (); # All the viral genomes belonging to this $season
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
		$Season2host_family2abun{$season}{$host_family} += $Viral_gn2IMG2abun{$viral_gn}{$img_id};
	}
}

# Step 8 Write down Season2host_family2abun.txt
open OUT, ">Season2host_family2abun.txt";
my $row2=join("\t", @Season);
print OUT "Head\t$row2\n";
foreach my $tmp1 (sort keys %Host_family){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %Season2host_family2abun) {       
                if (exists $Season2host_family2abun{$tmp2}{$tmp1}){
                        push @tmp, $Season2host_family2abun{$tmp2}{$tmp1};
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

