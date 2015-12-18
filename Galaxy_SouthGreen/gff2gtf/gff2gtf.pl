#!/usr/local/bin/perl
#-*-coding:UTF-8-*-

use strict;
use Carp;
use Getopt::Long;

#############
# Variables #
#############

my @f;
#- array   : GFF fields
my $gene;
#- string  : gene identifier
my $transcript;
#- string  :
my $coding_start;
#- integer : marker for finding start codon of CDS
my $coding_end;
#- integer : marker for finding end codon of CDS
my $parsed_CDS_exon;
#-         :
my %exon_no;
#- hash    :
my $ATG_start;
#- string  :
my $ATG_end;
#- string  :
my $STOP_start;
#- string  :
my $STOP_end;
#- string  :
my %exon_limit;
#- hash    :
my $pseudo;
#- integer : marker for pseudogene tag

########################
# Command-line options #
########################

my $infile;
#- string  : input GFF file
my $outfile;
#- string  : output GTF file
my $verbose = 1;
#- scalar  : Set for verbose output

#&GetOptions('file:s'=> \$infile,'output:s'=> \$outfile,'verbose'=> \$verbose);
$infile = shift;
$outfile = $infile.".gtf";


####################
# trap some errors #
####################

if ( $infile eq "" ) {
	print "// gff2gtf.pl : No input file name supplied. Nothing to parse, Using ".$ARGV[0]."\n";
	#$infile = shift;
	#$outfile = shift;
	print "// gff2gtf.pl : $infile\n";
	#    exit;
}

if ( $outfile eq "" ) {
	print "// gff2gtf_.pl : No output file name supplied. Using default - output.gff\n" if ( $verbose);
	exit;
}

#################################################################
# calculate the number of exons and coding exons per transcript #
#################################################################


my ($no_exons,$no_cds_exons) = &get_data($infile);
my %no_exons = %$no_exons;
my %no_cds_exons = %$no_cds_exons;



##------------------- MAIN LOOP ---------------------------##

############################
# Open outfile to write to #
############################
open (OUT, ">$outfile") or die "Failed to open output file: $outfile\n";

##########################################
# Open the GFF file and start the rename #
##########################################
open (GFF, "$infile") or die "cannot open $infile\n";
while (<GFF>) {

	chomp;
	@f = split (/\t/);

	# parse some fields from the free text information
	($transcript) = $f[8] =~ (/Parent=(\S+)\;/);
	
	
	
	# ignore polypeptide lines
	if ( $f[2] eq "polypeptide" ) {
		($pseudo) = $f[8] =~ (/pseudogene\=1/);
	}	
	
	
	########
	# gene #
	########
	if ( $f[2] eq "gene" ) {
	
		# get some data from the free text field
		($gene) = $f[8] =~ (/ID=(\S+)\;/);
		
		# Verbose/debug line
		print "\n// New gene $gene\n"if ( $verbose );
		
	}

	########
	# mRNA #
	########
	if ( $f[2] eq "mRNA" ) {
		
		# get some data from the free text field
		($transcript) = $f[8] =~ (/ID=(\S+)\;/);

		# Verbose/debug line
		print "// New transcript $transcript [$no_cds_exons{$transcript} coding exons of $no_exons{$transcript} total exons]\n"if ( $verbose );

		# Housekeeping
		undef $coding_start;
		undef $coding_end;
		undef $pseudo;
		$parsed_CDS_exon = 0;

	} # end of mRNA parsing

	########
	# Exon #
	########
	if ( $f[2] eq "exon" ) {

		$parsed_CDS_exon++;

		# get some data from the free text field
		($transcript) = $f[8] =~ (/Parent=(\S+)/);


		# verbose/debug line
		print "// $f[0] $f[2] $f[3] - $f[4] ($f[6]) [ $parsed_CDS_exon of $no_exons{$transcript} ] $f[8]\n" if ( $verbose );

		##############
		# +ve strand #
		##############
		if ( $f[6] eq "+" ) {
			#print "// Exon on +ve strand\n";

			$f[8] = $f[8] . " exon_number=$parsed_CDS_exon;";

			# First exon
			unless ( $coding_start ) {
				$ATG_start = $f[3];
				$ATG_end = $f[3] + 2;
				print OUT "$f[0]\t$f[1]\tstart_codon\t$ATG_start\t$ATG_end\t$f[5]\t$f[6]\t$f[7]\t$f[8]\n";
				$coding_start = 1;
			} #_ end of first exon block

			# Final exon
			unless ( $coding_end ) {
				if ( $parsed_CDS_exon == $no_cds_exons{$transcript} ) {
					$STOP_start = $f[4] - 2;
					$STOP_end = $f[4];
					print OUT "$f[0]\t$f[1]\tstop_codon\t$STOP_start\t$STOP_end\t$f[5]\t$f[6]\t$f[7]\t$f[8]\n";
					$coding_end = 1;
				}
			} #_ end of final exon block

		} #_ end of +ve strand parsing

		################
		## -ve strand ##
		################
		if ( $f[6] eq "-" ) {

			$f[8] = $f[8]." exon_number=".(($no_cds_exons{$transcript} - $parsed_CDS_exon) +1) . ";";

			# First exon
			unless ( $coding_end ) {

				$STOP_start = $f[3];
				$STOP_end = $f[3] + 2;
				print OUT "$f[0]\t$f[1]\tstop_codon\t$STOP_start\t$STOP_end\t$f[5]\t$f[6]\t$f[7]\t$f[8]\n";
				$coding_end = 1;

			} #_ end of first exon block

			# Final exon
			unless ( $coding_start ) {

				if ( $parsed_CDS_exon == $no_cds_exons{$transcript} ) {
					$ATG_start = $f[4] - 2;
					$ATG_end = $f[4];
					print OUT "$f[0]\t$f[1]\tstart_codon\t$ATG_start\t$ATG_end\t$f[5]\t$f[6]\t$f[7]\t$f[8]\n";
					$coding_start = 1;
				}
			} #_ end of final exon block

		} #_ end of -ve strand parsing

	} #_ end of CDS parsing

	#  polypeptide lines (get pseudogene information)

	if ( $f[2] eq "polypeptide" ) {
		($pseudo) = $f[8] =~ (/pseudogene\=1/);

		# this gets printed to the GTF file already (ignore for now)
	}

	# print GTF format line
	print OUT "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\t$f[8]\n";
	if ( $f[2] eq "exon" ) {
		print OUT "$f[0]\t$f[1]\tCDS\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\t$f[8]\n";
	}
}
close GFF;

close OUT;

# hasta luego
exit(0);


############################################
# get exon & CDS information from GFF file #
############################################
sub get_data {

my $infile = shift;
#- string  : path of the input GFF file
my @f;
#- array   : GFF fields
my %no_exons;
#- hash    : number of exons keyed by transcript name
my %no_cds_exons;
#- hash    : number of coding exons keyed by transcript name

open (GFF, "<$infile");
while (<GFF>) {

	chomp;
	@f = split (/\t/);

	# parse exon information
	if ( $f[2] eq "exon") {
		if ( $f[8] =~ /Parent=(\S+)/ ) {
			$no_exons{$1}++;
		}
	}

	# parse CDS information
	if ( $f[2] eq "exon") {
		if ( $f[8] =~ /Parent=(\S+)/ ) {
			$no_cds_exons{$1}++;
		}
	}
}
close GFF;

return (\%no_exons,\%no_cds_exons);

} #_ end of get_data sub
