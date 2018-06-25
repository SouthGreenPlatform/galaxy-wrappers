#!/usr/bin/perl

use warnings;
use strict;

#Command line processing.
use Getopt::Long;

my $input;
my $output;
Getopt::Long::Configure ('bundling');
GetOptions ('i|input=s' => \$input,
	    'o|output=s' => \$output)
or die("Error in command line arguments\n");

open(F1,"$input") or die ("Erreur d'ouverture de F1\n");

my @chro;
my $chr;
my $data;
while (my $li = <F1>){
	chomp($li);
	if($li !~ m/^#(.+)/){
		my@l = split(/\t/, $li);
		if(!grep { $_ eq $l[0] } @chro){
			push (@chro, $l[0]);
		}
	}
}
close F1;

foreach my$i(@chro){
	$data = "";
	open(F2,">$i.vcf") or die ("Erreur d'ouverture de F2\n");
	open(F1,"$input") or die ("Erreur d'ouverture de F1\n");
	while (my $li = <F1>){
		chomp($li);
		if($li =~ m/^#(.+)/){
			if($1 =~ m/#contig=<ID=(\S+),.+/){
				if ($1 eq $i){
					print F2 "$li\n";
				}
			}
			else{
				print F2 "$li\n";
			}
		}
		else{
			if($li =~ m/^$i/){
				print F2 "$li\n";
			}
		}
	}
	close F1;
	close F2;
}




