#!/usr/bin/perl -w

#This script takes an input VCF file and discard individuals with missingness percentage superior to INT (in parameter)
#Yann Hueber, july 2014

use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
use Carp;


# Forward declarations
sub parse_command_line();
sub open_vcf_out_file($);
sub load_vcf_in_file($);
sub indv_with_high_missingness();
sub get_ac_af_an_dp_indvstring($$);
sub print_results(@);
sub usage();

#Global arguments set by command line
my $vcf_in;
my $vcf_out = "output.vcf";
my $max_missingness_perc = 100;
my $ploidy;

#Gobal variables
my @individuals;
my @indv_to_del; #remove those indv from the vcf file
my $filehandle;
my %snps;
my $count_snps = 0;
my %count_missing;
my $headline ="";


################### Start of program #####################

parse_command_line();

open_vcf_out_file( $vcf_out );

load_vcf_in_file ( $vcf_in );

@indv_to_del = indv_with_high_missingness();

if (scalar (@indv_to_del) == 0) {
	#print ("No individuals with missingness superior to $max_missingness_perc\n");
	my $cp_vcf_in = "cp $vcf_in $vcf_out";
	system("$cp_vcf_in");
}
else {
	#print "HEY THERE we will discard some individuals\n";
	foreach ( @indv_to_del ) { print "Indiv to delete $individuals[$_] \n"; }

	print_results( @indv_to_del );
}


#################### End of program ######################




sub parse_command_line() {
	my $help;

	usage() if (scalar @ARGV==0);

	my $result = GetOptions ( "i|vcfin=s" => \$vcf_in,
				  "o|vcfout=s" => \$vcf_out,
				  "m|max_missingness_percentage=i" => \$max_missingness_perc,
				  "p|ploidy=i" => \$ploidy,
				  "h|help" => \$help);

	usage() if ($help);

	die "Error: vcf file not specified (use '-i or --vcfin [FILENAME]')\n" unless defined $vcf_in;
	
	die "Error invalid maximum missingness percentage (valid values are from 0 to 100)\n" if ($max_missingness_perc < 0 or $max_missingness_perc > 100);

	die "Error invalid ploidy (valid values are 2 or 3)\n" if (! ($ploidy == 2 or $ploidy == 3));

	die "Error in command line aguments\n" unless $result;
}





sub open_vcf_out_file($) {
	my $filename_out = shift or croak("Missing vcf_out file name");
	open($filehandle,">$filename_out") or die("Cannot open file $filename_out : $!");
}





sub load_vcf_in_file($) {
	my $filename_in = shift or croak("Missing vcf_in file name");
	
	open(VCF_IN,"<$filename_in") or die("Cannot open vcf file $filename_in : $!\n");
	
	croak("Bad file handle") unless(defined($filehandle));

	while (my $vcf_line = <VCF_IN>) {
		

		chomp($vcf_line);
		my @infos_line = split(m/\t/,$vcf_line);
		
		#File infos
		if ($vcf_line =~ m/^#/) {
			if ($vcf_line !~ m/^#CHROM/) {
				$headline.=  $vcf_line . "\n";
			} else {
				@individuals = @infos_line[9..$#infos_line];
			}
			next;
		}
		
		#Snps
		$count_snps++;
		#my $id = $infos_line[2];
		my $chr = $infos_line[0];
		my $chr_number = $chr;
		$chr_number =~ s/chr//;
		my $position = $infos_line[1];
		my $id = "S" . $chr_number . "_" . $position;
		
		my @snp_infos = @infos_line[0..8];
		my @individuals_infos = @infos_line[9..$#infos_line];
		$snps{$chr}{$position}{"snp_infos"} = \@snp_infos;
		$snps{$chr}{$position}{"indvs_infos"} = \@individuals_infos;
		
		
		#Count missing per indv
		my $i=0;
		foreach my $indv_info (@individuals_infos) {
			my @info_split = split(m/:/,$indv_info);
			my $genotype = $info_split[0];
			if ($genotype =~ m/\.\/\./) {
				$count_missing{$individuals[$i]}++;
			}
			$i++;
		}
	}
	close(VCF_IN);
}



sub indv_with_high_missingness() {
	
	my @indv_to_del;

	my $num_indiv=0;
	foreach my $indv (@individuals) {
		my $perc_miss = ($count_missing{$indv} / $count_snps) * 100;
		#print "HELLO indiv $num_indiv with percent miss : $perc_miss\n";
		if ($perc_miss > $max_missingness_perc) {

			#print "HERE : indiv to del : indiv number $num_indiv with percent miss : $perc_miss\n";
			push(@indv_to_del,$num_indiv);
		}
		$num_indiv++;
	}
	return @indv_to_del;
}



sub print_results(@) {
	
	my @indv_to_delete = @_;
	#foreach (@indv_to_delete) {
		#print "\nCOUCOU INdiv to delete : $_ \n";
	#}



	my %bad_indvs = map { $_ => 1 } @indv_to_delete;

	#print headlines with ##
	print $filehandle $headline;

	#print #CHROM headline
	my $chrom_header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	my $incr=0;
	foreach (@individuals) {
		if (! $bad_indvs{$incr}) {
			$chrom_header_line .= "\t$_";
		}
		$incr++;
	}
	print $filehandle $chrom_header_line . "\n";

	#SNP by SNP
	foreach my $chr (sort {$a <=> $b} keys (%snps)) {
		
		my $pos_ref = $snps{$chr};
		my %positions = %$pos_ref;
	
		foreach my $pos (sort {$a <=> $b} keys (%positions)) {

			my $snp_infos_ref = $snps{$chr}{$pos}{"snp_infos"};
			my @snp_infos = @$snp_infos_ref;
			my $indvs_infos_ref = $snps{$chr}{$pos}{"indvs_infos"};
			my @indvs_infos = @$indvs_infos_ref;

			my ($AC,$AF,$AN,$DP,$genotypes_string,$variant) = get_ac_af_an_dp_indvstring(\@indvs_infos,\%bad_indvs);
			next if (! $variant);

			#Line to print
			my $line_to_print="";
			my $chr = $snp_infos[0];
			my $position = $snp_infos[1];
			my $id_snp = $snp_infos[2];
			my $ref = $snp_infos[3];
			my $alt = $snp_infos[4];
			my $qual = $snp_infos[5];
			my $filter = $snp_infos[6];
			my $infos_updated = "AC=" . $AC . ";AF=" . $AF . ";AN=" . $AN . ";DP=" . $DP;
			my $previous_infos = $snp_infos[7];
			my @previous_infos_split = split(m/;/,$previous_infos);
			foreach (@previous_infos_split) {
				next if ($_ =~ m/A[CFN]|DP/);
				$infos_updated .= $_;
			}
			my $format = $snp_infos[8];
			
			$line_to_print.= "$chr\t$position\t$id_snp\t$ref\t$alt\t$qual\t$filter\t$infos_updated\t$format\t$genotypes_string\n";
			print $filehandle $line_to_print;
		}
	}
	close($filehandle);
}



sub get_ac_af_an_dp_indvstring($$) {
	my $indiv_infos_ref = shift;
	my @indivs_infos = @$indiv_infos_ref;
	my $bad_indiv_ref = shift;
	my %bad_indiv = %$bad_indiv_ref;

	my %allele_count;
	my $an=0; #total allele number
	my $total_DP=0;
	my $indvs_genotype_string="";


	for (my $i=0; $i< scalar(@indivs_infos); $i++) {
		next if (exists($bad_indiv{$i}));

		my $indiv_infos = $indivs_infos[$i];
		$indvs_genotype_string.= $indiv_infos . "\t";

		my @indiv_infos_split = split(m/:/,$indiv_infos);
		my $genotype = $indiv_infos_split[0];
		next if ($genotype =~ m/\.\/\./);

		#allele total number and allele count
		$an+=2;
		my $first_allele = substr($genotype,0,1);
		my $second_allele = substr($genotype,2,1);
		$allele_count{$first_allele}++;
		$allele_count{$second_allele}++;
		if ($ploidy == 3) {
			my $third_allele = substr($genotype,4,1);
			$allele_count{$third_allele}++;
			$an+=1;
		}

		
		#read depth
		my $DP = $indiv_infos_split[2];
		$total_DP+=$DP;

	}
	chop($indvs_genotype_string);


	#SNP still a variant after removing individuals with high missingness?
	my $is_variant=0;
	if (scalar(keys(%allele_count)) >= 2) {
		$is_variant=1;
	}

	#allele count and allele frequency
	my $ac="";
	my $af="";
	for (my $i = 1; $i<= 3; $i++) {
		if (exists($allele_count{$i})) {
			$ac.= $allele_count{$i} . ",";
			$af.= sprintf("%.2f", $allele_count{$i} / $an) . ",";
		}
	}
	chop($ac);
	chop($af);

	return($ac,$af,$an,$total_DP,$indvs_genotype_string,$is_variant);
}


sub usage() {
print<<EOF;

This program reads a VCF file and remove individuals based 
on missingness provided in the parameters

usage: perl $0 -i VCF_FILENAME_IN -m MISSINGNESS -p PLOIDY [-o VCF_FILENAME_OUT]

Arguments:

-i|--vcfin FILE	                      - VCF_IN filename

-m|--max_missingness_percentage INT   - Set the maximum missingness percentage allowed to keep individuals

-o|--vcfout FILENAME		      - VCF_OUT filename

-p|--ploidy INT			      - Set the ploidy (2 or 3)

--help 				      - This helpful help screen

EOF

exit 1;
}



