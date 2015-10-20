#!/usr/bin/perl -w


#This script takes an input VCF file from Cornell TASSEL pipeline (plugin tbt2vcfPlugin) and makes it standard, that is to say :
#change major/minor SNPs format to reference/alternate
#correct deletion (ex : ref A alt -) to standard format (ref TA alt T)
#Yann Hueber, June 2014

use warnings;
use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;


#Forward declarations
sub parse_command_line();
sub open_vcf_out_file($);
sub load_ref_file($);
sub runthrough_vcf_and_print_out_file($);
sub getPrecedentalleles($$);


#Global arguments set by command line
my $input_vcf;
my $output_vcf = "output.vcf";
my $ref_file;


#Global variables
my %seq_hash;
my $seqio_obj;
my $filehandle;

################### Start of program #####################

parse_command_line();

open_vcf_out_file( $output_vcf );

load_ref_file( $ref_file );

runthrough_vcf_and_print_out_file( $input_vcf );


#################### End of program ######################



sub parse_command_line() {
	my $help;

	usage() if (scalar @ARGV==0);
	
	my $result = GetOptions ("i|vcfin=s" => \$input_vcf,       #string
				 "r|refin=s" => \$ref_file,        # string
				 "o|output=s" => \$output_vcf,      #string
				 "h|help" => \$help);

	usage() if ($help);

	die "Error vcf file not specified (use '-i or --vcfin [FILENAME]')\n" unless defined $input_vcf;
		
	die "Error reference fasta file not specified (use '-r or --ref_in [FILENAME]')\n" unless defined $ref_file;

	die "Error in command line aguments\n" unless $result;
}




sub open_vcf_out_file($) {
	my $filename_out = shift or croak("Missing output_vcf file name");
	open($filehandle,">$filename_out") or die("Cannot oepn file $filename_out : $!");
}



#retrieve seqIO object with reference file and charge it in memory
sub load_ref_file ($) {
	$seqio_obj = Bio::SeqIO->new(-file => "<$ref_file", -format => "fasta");
	while (my $seq_obj = $seqio_obj->next_seq) {
		my $seq_name = $seq_obj->display_id;
		$seq_hash{$seq_name} = $seq_obj;
	}
}





#Running through Cornell VCF file and print a new standard VCF
sub runthrough_vcf_and_print_out_file($) {
	my $filename_in = shift or croak("Missing input_vcf filename");
	
	open(VCF_IN,"<$filename_in") or die("Cannot open vcf file $filename_in : $!\n");

	croak("Bad filehandle") unless(defined($filehandle));


	while (my $vcf_line = <VCF_IN>) {
	
		#header lines
		if ($vcf_line =~ m/^#/) {
			print $filehandle $vcf_line;
			next;
		}
	
		#sites (SNPs)
		chomp($vcf_line);
		my @infos_line = split(m/\t/,$vcf_line);
		my $chr = $infos_line[0];
		my $pos = $infos_line[1];
		my $major_allele = $infos_line[3];
		my $minor_allele = $infos_line[4];
		my $ref_allele = $seq_hash{$chr}->subseq($pos,$pos);
		my $new_line="";

		if ($major_allele eq $ref_allele and $minor_allele !~ m/-/) {
			$new_line .= $vcf_line . "\n";
		} else {
			#variables
			my ($new_ref,$new_id,$new_ref_string,$precedent_ref_allele,$new_position);
			my $new_alt_string="";
			my @minor_allele_split = split(m/,/,$minor_allele);

			#Major allele equal to ref allele --> No problem
			if ($major_allele eq $ref_allele) {
				($new_position,$precedent_ref_allele) = getPrecedentalleles($chr,$pos);
				$new_ref_string = $precedent_ref_allele . $ref_allele;
				$new_id = "S" . $chr . "_$new_position";
				foreach(@minor_allele_split) {
					if ($_ eq "-") {
						$new_alt_string .= $precedent_ref_allele . ",";
					} else {
						$new_alt_string .= $precedent_ref_allele . $_ . ",";
					}
				}
				chop($new_alt_string);
				my $end_line;
				foreach(@infos_line[5..7]) {
					$end_line .= $_ . "\t";
				}
				foreach(@infos_line[8..$#infos_line]) {
					if ($minor_allele =~ m/,/ and $_ =~ m/^GT/) {
						$end_line .= "GT:AD:DP:GQ\t";
					} elsif ($minor_allele =~ m/,/) {
						my @individual_genotype_split = split(m/:/,$_);
						$end_line .= $individual_genotype_split[0] . ":" . $individual_genotype_split[1] . ":" .  $individual_genotype_split[2] . ":" . $individual_genotype_split[3] . "\t";
					} else {
						$end_line .= $_ . "\t";
					}
				}
				chop($end_line);
				$new_line .= "$chr\t$new_position\t$new_id\t$new_ref_string\t$new_alt_string\t$end_line\n";
			} else {
	
				my %genotype_link;
	
				#new position and new ref string
				if ($major_allele eq "-" or $minor_allele =~ m/-/) {
					($new_position, $precedent_ref_allele) = getPrecedentalleles($chr,$pos);
					$new_ref_string = $precedent_ref_allele . $ref_allele;
				} else {
					$new_position = $pos;
					$new_ref_string = $ref_allele;
				}
				
				#print "ICI : precendent : $precedent_ref_allele et ref allele $ref_allele\n";
	
	
				#new alternate string
				my $i =1;
				foreach(@minor_allele_split) {
					#print "LA $_ equivaut a $ref_allele\n";
					next if ($minor_allele eq ".");
					if ($_ eq "-") {
						$new_alt_string .= $precedent_ref_allele . ",";
						$genotype_link{$i} = $i;
					} elsif ($_ eq $ref_allele) {
						$genotype_link{$i} = "0";
						$genotype_link{"0"} = $i;
						if ($major_allele eq "-" or $minor_allele =~ m/-/) {
							if ($major_allele eq "-") {
								$new_alt_string .= $precedent_ref_allele . ",";
							} else {
								$new_alt_string .= $precedent_ref_allele . $major_allele .",";
							}
						} else {
							$new_alt_string .= $major_allele . ",";
						}
					}  else {
						$genotype_link{$i} = $i;
						if ($major_allele eq "-" or $minor_allele =~ m/-/) {
							$new_alt_string .= $precedent_ref_allele . $_ . ",";
						}  else {
							$new_alt_string .= $_ . ",";
						}
					}
					$i++;
				}
				if ($minor_allele !~ m/$ref_allele/) { #in case of the reference is neither in major allele nor in minor allele 
					$genotype_link{$i} = "0";
					$genotype_link{"0"} = $i;
					if ($major_allele eq "-") {
						$new_alt_string .= $precedent_ref_allele . ",";
					}
					else {
						if ($major_allele ne $ref_allele and $minor_allele =~ m/-/) {
							$new_alt_string .= $precedent_ref_allele . $major_allele . ",";
						} else {
							$new_alt_string .= $major_allele . ",";
						}
					}
				}
				chop($new_alt_string);
	
	
				#New line beginning
				my $new_format;
				if ($new_alt_string =~ m/,/) {
					$new_format = "GT:AD:DP:GQ";
				} else {
					$new_format = $infos_line[8];
				}
				$new_line .= $chr . "\t" .$new_position . "\tS" . $chr . "_$new_position\t$new_ref_string\t$new_alt_string\t$infos_line[5]\t$infos_line[6]\t$infos_line[7]\t$new_format\t";
				
	
	
				#New line end
				my $end_line="";
				foreach my $individual_genotype(@infos_line[9..$#infos_line]) {
					
					my $new_individual_genotype;
					my @individual_genotype_split = split(m/:/,$individual_genotype);
	
	
					#dealing with genotypes
					my $old_genotype = $individual_genotype_split[0];
					if ($old_genotype eq "./.") {
						if ($new_alt_string =~ m/,/) {
							$new_individual_genotype = $individual_genotype_split[0] . ":";
							my @split_new_alt_string = split(m/,/,$new_alt_string);
							my $new_indiv_DP = "0,"; #for the reference
							foreach (@split_new_alt_string) { #for alt 
								$new_indiv_DP .= "0,";
							}
							chop($new_indiv_DP);
							
							$new_individual_genotype .= $new_indiv_DP . ":" . $individual_genotype_split[2] . ":" . $individual_genotype_split[3];
						} else {
							$new_individual_genotype = $individual_genotype;
						}
					} else {
						my @old_genotype_split = split(m/\//,$old_genotype);
	
						my @new_genotype = ($genotype_link{$old_genotype_split[0]},$genotype_link{$old_genotype_split[1]});
						@new_genotype = sort ( {$a <=> $b } @new_genotype );
						my $new_genotype = $new_genotype[0] . "/" . $new_genotype[1];
						
						#dealing with individual DP
						my $old_indiv_DP = $individual_genotype_split[1];
						my @old_indiv_DP_split = split(m/,/,$old_indiv_DP);
						my @new_indiv_DP;
						my $j=0;
						while ($j < scalar(@old_indiv_DP_split)) {
							if (! defined ($old_indiv_DP_split[$genotype_link{$j}])) { #in case of referece has been added
								@new_indiv_DP = ();
								$old_indiv_DP_split[$genotype_link{$j}] = "0";
								foreach my $key(keys(%genotype_link)) {
									$new_indiv_DP[$key] = $old_indiv_DP_split[$genotype_link{$key}];
								}
								$j = scalar(@old_indiv_DP_split);
							} else {
								$new_indiv_DP[$j] = $old_indiv_DP_split[$genotype_link{$j}];
								$j++;
							}
						}
						my $new_indiv_DP_string = join(",",@new_indiv_DP);
		
						#dealing with individual PL // Not applicable if site is not biallelic
						my $new_indiv_PL_string = "";
						if ($new_alt_string !~ m/,/) {
							my $old_indiv_PL = $individual_genotype_split[4];
							my @old_indiv_PL_split = split(m/,/,$old_indiv_PL);
							if ($old_genotype eq "0/0" or $old_genotype eq "1/1") {
								$new_indiv_PL_string .= $old_indiv_PL_split[2] . "," . $old_indiv_PL_split[1] . "," . $old_indiv_PL_split[0];
							} else {
								$new_indiv_PL_string = $old_indiv_PL;
							}
						}
					
						$new_individual_genotype = $new_genotype . ":" . $new_indiv_DP_string . ":" . $individual_genotype_split[2] . ":" . $individual_genotype_split[3];
						if ($new_indiv_PL_string ne "") {
							$new_individual_genotype .= ":" . $new_indiv_PL_string;
						}
					}
					$end_line .= $new_individual_genotype . "\t";
				}
				chop($end_line);
				$new_line .= $end_line . "\n";
			}
		}
		#print new line in output
		print $filehandle $new_line;
	}
	
	close(VCF_IN);
	close($filehandle);
}	
	





sub getPrecedentalleles($$) {
	my ($chrom,$current_position) = @_;
	my @precedents =();
	my $precedent_position = $current_position -1;
	my $precedent_ref_allele = $seq_hash{$chrom}->subseq($precedent_position,$precedent_position);
	push(@precedents,$precedent_position,$precedent_ref_allele);
	return @precedents;
}






sub usage() {
print<<EOF;
This script takes an input VCF file from Cornell TASSEL pipeline (plugin tbt2vcfPlugin) and makes it standard, that is to say :
#change major/minor SNPs format to reference/alternate
#correct deletion (ex : ref A alt -) to standard format (ref TA alt T)

usage: perl $0 -i CORNELL_VCF_IN -r REFERENCE.FA [-o VCF_FILENAME_OUT]

Arguments:

-i|vcf_in FILE			- VCF_IN filename

-r|ref_in FILE			- Reference filename(in FASTA)

-o|vcf_out			- VCF_OUT filename

--help				- This helpful help screen

EOF

exit 1;
}




