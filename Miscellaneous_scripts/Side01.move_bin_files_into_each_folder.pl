#!/usr/bin/perl

use strict;
use warnings;

my %MAG2IMGID = (); #$mag => $img_id
my %IMGID = (); #$img_id => 1
open IN, "ls *.tar.gz |";
while (<IN>){
	chomp;
	my $file = $_;
	my ($mag) = $file =~ /^(.+?)\.tar\.gz/;
	my ($img_id) = $mag =~ /^(.+?)\_/;
	$MAG2IMGID{$mag} = $img_id;
	$IMGID{$img_id} = 1;
}
close IN;

my $num_of_folders = 0;
$num_of_folders = scalar (keys %IMGID);
print "$num_of_folders\n";

foreach my $img_id (sort keys %IMGID){
	if (-d "$img_id"){
		open IN, "ls ${img_id}\_* |";
		while (<IN>){
			chomp;
			my $file = $_;
			if ($file !~ /\.tar\.gz/){
				`mv $file $img_id`;
			}
		}
		close IN;
=pod		
		open IN, "ls mbin_datafile_*.txt |";
		while (<IN>){
			chomp;
			my $file = $_;
			if ($file !~ /\.tar\.gz/){
				`mv $file $img_id`;
			}
		}
		close IN;	
=cut		
	}else{
		`mkdir $img_id`;
		open IN, "ls ${img_id}\_* |";
		while (<IN>){
			chomp;
			my $file = $_;
			if ($file !~ /\.tar\.gz/){
				`mv $file $img_id`;
			}
		}
		close IN;	
=pod	
		open IN, "ls mbin_datafile_*.txt |";
		while (<IN>){
			chomp;
			my $file = $_;
			if ($file !~ /\.tar\.gz/){
				`mv $file $img_id`;
			}
		}
		close IN;	
=cut		
	}
}

open IN, "ls ../ -l | grep ^d| grep '33' | sed -r s'/ +/\t/g' | cut -f 9 |";
while (<IN>){
	chomp;
	my $line = $_;
	if (!exists $IMGID{$line}){
		print "$line\n";
	}
}
close IN;



