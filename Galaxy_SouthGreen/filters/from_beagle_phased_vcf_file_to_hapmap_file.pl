#!/usr/bin/perl i-w

#Take an input VCF file and print a hapmap file (for one site in the vcf file, make two in the output hapmap format to keep phaasing informations)
#Check IUPAC hash to see how genotypes are coded

use warnings;
use strict;
use Getopt::Long;

#variables
my ($input_phasedvcf, $out_hapmap, $out_del_memory);
my @individuals;


#Modified IUPAC Code
my %IUPAC = ("AA" => "A",
	     "TT" => "T",
	     "CC" => "C",
	     "GG" => "G",
	     "AG" => "R",
	     "GA" => "R",
	     "CT" => "Y",
	     "TC" => "Y",
	     "AC" => "M",
	     "CA" => "M",
	     "TG" => "K",
	     "GT" => "K",
	     "AT" => "W",
	     "TA" => "W",
	     "CG" => "S",
	     "GC" => "S",
	     "A-" => "W",
	     "-A" => "W",
	     "T-" => "W",
	     "-T" => "W",
	     "G-" => "S",
	     "-G" => "S",
	     "C-" => "S",
	     "-C" => "S");

my %compl = ("A" => "T",
	     "T" => "A",
	     "C" => "G",
	     "G" => "C");


if (scalar(@ARGV) == 0) {
	print "\n\nUsage : perl $0 -i vcf_in -o hapmap_output -d how_del_was_coded_output\n\n\n";
	exit;
}

GetOptions ("i|vcf_in=s" => \$input_phasedvcf,	#string
	    "o|out=s" => \$out_hapmap,		#string
	    "d|del_out=s" => \$out_del_memory)	#string
or die ("Error in command line arguments\n");

die "Error: phased vcf not specified (use '-i [FILENAME]')\n" unless defined $input_phasedvcf;
die "Error: hapamp filename not specified (use '-o [FILENAME]'\n" unless defined $out_hapmap;
die "Error: txt filename (how deletions are coded not specified (use '-d [FILENAME]'\n" unless defined $out_del_memory;

#Open files
open(VCF_IN,"<$input_phasedvcf") or die("Cannot open file $input_phasedvcf : $!");
open(HAP,">$out_hapmap") or die("Cannot open file $out_hapmap : $!");
open(DELETION,">$out_del_memory") or die("Cannot open file $out_del_memory : $!");

#head line hapmap and txt files
print HAP "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode";
print DELETION "CHROM\tPOS\tREF\tALT\n";

#Running through vcf file
while (defined (my $vcf_line = <VCF_IN>)) {
	next if ($vcf_line =~ m/^##/);
	chomp($vcf_line);
	my @infos_line = split(m/\t/,$vcf_line);
	my @individuals_infos = @infos_line[9..$#infos_line];

	#Print individuals name in hapmap file
	my $individuals_string="";
	if ($vcf_line =~ m/^#CHROM/) {
		@individuals = @individuals_infos;
		foreach my $indiv (@individuals) {
			$individuals_string .= "\t" . $indiv . "A" . "\t" . $indiv . "B";
		}
		print HAP $individuals_string . "\n";
		next;
	}
	
	#For one site
	my $chr = $infos_line[0];
	my $pos = $infos_line[1];
	my $ref = $infos_line[3];
	my $alt = $infos_line[4];

	#Change deletion and save changes in a file
	if (length($ref) > 1) {
		print DELETION "$chr\t$pos\t$ref\t$alt\n";
		#Change ref (compl of alternate allele)! This is arbitrary
		$ref = $compl{$alt};
	}
	my @ref_alt_list = ($ref,$alt);

	#line_begins
	my $strand = "+";
	my $line_to_print = "S$chr\_$pos\t$ref/$alt\t$chr\t$pos\t$strand\tNA\tNA\tNA\tNA\tNA\tNA";

	foreach my $indiv_infos (@individuals_infos) { #a loop for each indiv
		my @indiv_infos_split = split(m/:/,$indiv_infos);
		my $genotype = $indiv_infos_split[0];
		my @genotype_split = split(m/\|/,$genotype);
		$line_to_print .=  "\t" . $ref_alt_list[$genotype_split[0]] . "\t" . $ref_alt_list[$genotype_split[1]];
	}
	#print line 
	print HAP $line_to_print . "\n";
}

close(VCF_IN);
close(HAP);
close(DELETION);

exit 1;
