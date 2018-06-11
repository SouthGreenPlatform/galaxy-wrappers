#!/usr/bin/perl

# Version 1 : 01/06/18


use warnings;
use strict;

#Command line processing.
use Getopt::Long;

#initialization of files and variables that can be chosen by the user to launch the script.
my $input;
my $vcf;
Getopt::Long::Configure ('bundling');
GetOptions ('t|input=s' => \$input, ## reference matrix
	    'v|vcf=s' => \$vcf) ##
or die("Error in command line arguments\n");


my%positions; #hash table containing
open(F1,"$input") or die ("Error: DSNP matrix wont open\n");  # opening the DSNP matrix 
while (my $li = <F1>){ # pour chaque dsnps
	chomp($li);
	if ($li !~ m/^ancestor/){
		my@ligne = split("\t", $li);
		#Par marqueur
		my$ancetre = $ligne[0];
		my$chromosome = $ligne[1];
		my$position = $ligne[2];
		my$allele = $ligne[3];
		
		$positions{$chromosome}{$position} = 1;
	}
}
close F1;

open(F2,"$vcf") or die ("Error: vcf wont open\n");  # opening the DSNP matrix and read the first line.
open(FW, ">newVCF");
my$chrVCF;
my$posVCF;
while (my $li = <F2>){ # pour chaque dsnps
	chomp($li);
	if ($li =~ m/^#/){
		print FW "$li\n";
	}
	elsif ($li =~ m/^(\S+)\t(\d+)\t.+/){
		$chrVCF = $1;
		$posVCF = $2;
		
		my@chr = split("chr0",$chrVCF);
		my$chromosome = $chr[1];
		print "$chromosome\n";
		#if ($positions{$chrVCF}{$posVCF} == 1 or $positions{"chr".$chrVCF}{$posVCF} or $positions{"chr0".$chrVCF}{$posVCF}){
			#print FW "$li\n";
		#} 
		if (exists($positions{$chromosome}{$posVCF}) && $positions{$chromosome}{$posVCF} == 1){
			print FW "$li\n";
		}
	}	
}
close F2;
close FW;

