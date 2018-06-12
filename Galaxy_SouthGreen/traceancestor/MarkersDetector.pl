#!/usr/bin/perl

# Version 1 : 01/06/18 
# This script is used to define a matrix of ancestry informative markers from GST informations.

use warnings;
use strict;

# Command line processing.
use Getopt::Long;

# Initialization of files and variables that can be chosen by the user to launch the script.
my $input;
my $gst = 0.9;
my $missing = 0.3;
my $help = "";
Getopt::Long::Configure ('bundling');
GetOptions ('t|input=s' => \$input, ## Reference matrix with GST and Alternative allele frequence (F) informations. The first column names schould be under the format: #CHROM POS ID REF ALT %Nref GST1 GST2 GST3 F1 F2 F3
	    'g|gst=f' => \$gst, ## value of GST (inter-population differentiation parameter) necessary for an allele to be an ancestry informative marker.
	    'm|missing=f' => \$missing, ## maximal frequency of missing data for an allele to be defined as an ancestry informative marker
	    'h|help!'=> \$help, ## display the help
		)
or die("Error in command line arguments\n");
if($help) { ## display this help when -help is called
	print "\nusage: script.pl [-t matrix file] [-g threshold for gst] [-m threshold for missing data] [-h help]\n\n";
	print '-t | --input : reference matrice. The first column names schould be under the format: #CHROM POS ID REF ALT %Nref GST1 GST2 GST3 F1 F2 F3 (1, 2, 3 ancestry names)'."\n";
	print "\t#CHROM = chromosome name\n";
    	print "\tPOS = position of DSNP\n";
	print "\tREF = Base of the reference allele of this DSNP\n";
	print "\tALT = Base of the alternative allele of this DSNP\n";
	print "\t".'%Nref = Percentage of maximal missing data for this DSNP'."\n";
	print "\tGST = Allele frequency for each ancestor (With 1,2,3 the ancestors names)\n";
	print "\tF = Alternative allele frequency for each ancestor (With 1,2,3 the ancestors names)\n";
	print "-g | --gst : threshold for gst\n";
	print "-m | --missing : threshold for missing data\n";
	print "-h | --help : display this help\n\n";
	exit;
}   

#creation of 2 lists: @ref and @header:
open(F1,"$input") or die ("Error: DSNP matrix wont open\n");  # open the DSNP matrix and read the first line.
my $first_line = <F1>;
chomp($first_line);
my @header = split("\t", $first_line); # header of the DSNP matrix 
close F1;

my@ref; # list of the name of GST for each ancestral genomes.
foreach my$head(@header){
	if ($head =~ m/^GST.+$/){
		push (@ref, $head);
	}
}

#------------------------------------------------------------------------------------------------------------------------------------------------------

# Selection of all ancestry informative markers for each ancestry and associated allele (alt or ref)

my%allele_refalt; # key 1 = parent / key 2 = chromosome / key 3 = position => value = allele
my$chro = ""; # chromosome of the read line 
my$pos = ""; # postiton of the read line 
my%isAMarker; # key 1 = parent / key 2 = chromosome / key 3 = position => value = 1 if the position is an marker for this parent
my$Bref; # value of reference allele for the position of the read line 
my$Balt; # value of alternative allele for the position of the read line 
open(F1,"$input") or die ("Erreur d'ouverture de F1\n"); #ouverture du fichier matrice
while (my $li = <F1>){ # pour chaque dsnps
	chomp($li);
	if ($li !~ m/^#CHROM/){
		my@ligne = split("\t", $li);
		#Par marqueur
		foreach my$r(@ref){
			my$ancetre;
			if ($r =~ m/^GST(.+)/){
				$ancetre = $1;
			}
			my$booMiss = 0;
			for (my$h = 0; $h < scalar(@header); $h++){
				if($header[$h] eq '#CHROM'){
					$chro = $ligne[$h];
				}
				elsif($header[$h] eq 'POS'){
					$pos = $ligne[$h];
				}
				elsif($header[$h] eq 'REF'){
					$Bref = $ligne[$h];
				}
				elsif($header[$h] eq 'ALT'){
					$Balt = $ligne[$h];
				}
				elsif($header[$h] eq '%Nref'){
					if($ligne[$h] <= $missing){
						$booMiss = 1;
					}
				}
				elsif ($header[$h] eq $r){
					# Selection of marker with the threshold of the missing data and the threshold of GST
					if($ligne[$h] >= $gst && $booMiss == 1){
						$isAMarker{$ancetre}{$chro}{$pos} = 1;							
					}
				}
				# We are wondering if the ancestral allele is the reference or the alternative allele. To see that, we are checking F values (alternative allele frequency). If > 0.8, ancestral allele is alternative allele. If < 0.2, ancestral allele is reference allele.
				elsif ($header[$h] =~ m/^F$ancetre$/){
					if (exists($isAMarker{$ancetre}{$chro}{$pos})){
						if($ligne[$h] > 0.8){ # alternative 
							$allele_refalt{$ancetre}{$chro}{$pos} = $Balt;
						}
						elsif ($ligne[$h] < 0.2){  # reference
							$allele_refalt{$ancetre}{$chro}{$pos} = $Bref;
						}
					}
				}

			 }

		}
	}
}
close F1;

# Fill the Diagnosis matrix
open(Fmat, ">Diagnosis_matrix");
print Fmat "ancestor\tchromosome\tposition\tallele\n"; # header
for my$ancetre (sort keys %isAMarker){
	for my$chro (sort keys %{$isAMarker{$ancetre}}){
		for my$pos (sort {$a<=>$b} keys %{$isAMarker{$ancetre}{$chro}}){
			if ($isAMarker{$ancetre}{$chro}{$pos} == 1){
				print Fmat join("\t", $ancetre, $chro ,$pos, $allele_refalt{$ancetre}{$chro}{$pos}),"\n";
			}
		}
	}
}
close Fmat;

