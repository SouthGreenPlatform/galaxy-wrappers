#!/usr/bin/perl

#blat_results_analyzer.pl - Analyse and return in text format the BLAT's results.

## command line : perl blat_results_analyzer.pl --blat "file_BLAST/BLAT" --identity "value_of _percentage_identity" --outfile "outfile_name"

## file_BLAST/BLAT => name of your BLAST/BLAT out file (in pairwise format).
## value_of _percentage_identity => value between 0 and 100 (ex : 80).
## outfile_name => the name desired for program outfile, just the name without extension. The outfile name will be named like this : BRA_"outfile name".txt

use strict;
use Carp qw(confess);
use Bio::SearchIO;
use Getopt::Long;
use Pod::Usage;
use warnings;

my $nombre;
my $help;
my $blat_result;
my $identity_percent;
# my $outfile;
my $blat_result_in;
my $result;
my %result_query_length;
my $hit;
my %result_hit_length;
my $hsp;
my %hsp_start_end;
my %hsp_strand;
my $verbose = 0;
GetOptions(
	   "help" => \$help,
	   "blat=s" => \$blat_result, 
	   "identity=i" => \$identity_percent,
	   # "outfile=s" => \$outfile,
	   "nombre=i" => \$nombre,
	   "verbose" => \$verbose);
if (defined $help) { pod2usage("\n#blat_results_analyzer.pl - Analyse and return in text format the BLAT's results.\n
## command line : perl blat_results_analyzer.pl --blat \"file_BLAST/BLAT\" --identity \"value_of _percentage_identity\" --nombre \"max number different hits\"\n
## file_BLAST/BLAT => name of your BLAST/BLAT out file (in pairwise format).
## value_of _percentage_identity => value between 0 and 100 (ex : 80).") && exit; }
unless(defined $blat_result) { pod2usage("\n!!! You must enter a blat result file with \"--blat\" option - for more informations, enter \"perl blat_results_analyzer.pl --help\" !!!\n") && exit; }
unless(defined $identity_percent) { pod2usage("\n!!! You must enter an identity percent with \"--identity\" option - for more informations, enter \"perl blat_results_analyzer.pl --help\" !!!\n") && exit; }
# unless(defined $outfile) { pod2usage("\n!!! You must enter an outfile name (without extension) with \"--outfile\" option - for more informations, enter \"perl blat_results_analyzer.pl --help\" !!!\n") && exit; }

$blat_result_in = new Bio::SearchIO(-format => 'blast', -file => $blat_result);

#Store four data structure : 
# - First : query with association to query's length.
# - Second : subject with association to subject's length.
# - Third : all data for each hsp.
# - Four : relative data for hsp strand orientation.
print "   1. Recovery of all BLAT/BLAST data.\n\n" if $verbose;
print "   Markers with no hits found : \n" if $verbose;
while ($result = $blat_result_in -> next_result){
    if ($result -> num_hits == 0){
	print " - ".$result -> query_name."\n" if $verbose;
    }
    $result_query_length{$result -> query_name} = $result -> query_length;
    while ($hit = $result -> next_hit){
	if ($hit -> length != 0){
	    $result_hit_length{$hit -> name} = $hit -> length;
	}
	else {
	    next;
	}
	while ($hsp = $hit -> next_hsp){
	    $hsp_start_end{$result -> query_name}{$hit -> name}{$hsp -> start('hit')}{$hsp -> end('hit')}{$hsp -> start('query')}{$hsp -> end('query')}{$hsp -> num_identical} = $hsp;# -> length('total');
	    $hsp_strand{$result -> query_name}{$hit -> name}{$hsp -> start('hit')} = $hsp -> strand('hit');
	}
    }
}

my $ref_HSE = \%hsp_start_end;
my $ref_HS = \%hsp_strand;

# unlink ("BRA_".$outfile."_marker_no_results.txt") if -e "BRA_".$outfile."_marker_no_results.txt";
# open (BMNR, ">>BRA_".$outfile."_marker_no_results.txt") or confess ("File BRA_".$outfile."_marker_no_results can't be created.");
# print (BMNR "Markers rejected at data presort : \n");

print "\n   2. Presorting data ...\n\n" if $verbose;
my @data_of_presort = &data_presort($ref_HSE, $ref_HS);
my $query_number_presort = $data_of_presort[2];

print "    * $query_number_presort markers rejected at data presort.\n\n" if $verbose;

# unlink ("BRA_".$outfile."_more_five_scaffolds.txt") if -e "BRA_".$outfile."_more_five_scaffolds.txt";
# open (MFS, ">>BRA_".$outfile."_more_five_scaffolds.txt") or confess ("File BRA_".$outfile."_more_five_scaffolds can't be created.");
# print (MFS "Markers who match on five or more scaffolds : \n");

# unlink ("BRA_".$outfile.".txt") if -e "BRA_".$outfile.".txt";
# open (FBR, ">>BRA_".$outfile.".txt") or confess ("File BRA_".$outfile.".txt can't be created.");
# print (FBR "Query\tHit_name\tHit_length\tHSP_number\tIdentity percent\tStart_position_subject\tEnd_position_subject\n");

# print (BMNR "\nMarkers without significant hit on the BLAST/BLAT : \n");

print "   3. Final sorting ...\n\n" if $verbose;
my @final_sorting_result = &final_sorting($data_of_presort[0], $data_of_presort[1], $ref_HS);
# close (BMNR);
# close (MFS);

my $number_query_without_significant_hit = $final_sorting_result[2];
my $query_with_many_scaffold = $final_sorting_result[1];
my $ref_MI = $final_sorting_result[0];
my %markers_infos = %$ref_MI;
my $number_of_good_markers = 0;

print "    * $number_query_without_significant_hit markers without significant hit on the BLAST/BLAT.\n" if $verbose;
print "    * $query_with_many_scaffold markers who match on five or more scaffold.\n\n" if $verbose;

print "   4. Results : \n\n" if $verbose;
foreach my $query_name (keys %markers_infos){
    $number_of_good_markers ++;
    foreach my $scaffold (keys %{$markers_infos{$query_name}}){
	foreach my $scaffold_length (keys %{$markers_infos{$query_name}{$scaffold}}){
	    foreach my $hsp_number_per_marker_localisation (keys %{$markers_infos{$query_name}{$scaffold}{$scaffold_length}}){
		foreach my $global_hsp_query_identity (keys %{$markers_infos{$query_name}{$scaffold}{$scaffold_length}{$hsp_number_per_marker_localisation}}){
		    foreach my $start (sort {$a <=> $b} keys %{$markers_infos{$query_name}{$scaffold}{$scaffold_length}{$hsp_number_per_marker_localisation}{$global_hsp_query_identity}}){
			foreach my $end (keys %{$markers_infos{$query_name}{$scaffold}{$scaffold_length}{$hsp_number_per_marker_localisation}{$global_hsp_query_identity}{$start}}){
			    my $strand = $markers_infos{$query_name}{$scaffold}{$scaffold_length}{$hsp_number_per_marker_localisation}{$global_hsp_query_identity}{$start}{$end}->strand('hit');
			   
			    # print (FBR "$query_name\t$scaffold\t$scaffold_length\t$hsp_number_per_marker_localisation\t$global_hsp_query_identity\t$start\t$end\t$strand\n");
			    print "$query_name\t$scaffold\t$scaffold_length\t$hsp_number_per_marker_localisation\t$global_hsp_query_identity\t$start\t$end\t$strand\n" unless $verbose;
			    
			}
		    }
		}
	    }
	}
    }
}
# close (FBR);

print "    * $number_of_good_markers markers who match correctly on $blat_result\n\n" if $verbose;
print "   ------------------------------------------------------\n\n" if $verbose;
print "    -> Analysis results are available in : BLAT_results_analyzer.txt\n    -> Markers who match on more of five scaffolds are available in : BLAT_results_analyzer_more_five_scaffolds.txt\n    -> Markers without hits are available in : BLAT_results_analyzer_marker_no_results.txt\n\n" if $verbose;
print "   Done!\n\n" if $verbose;

#Data presort, conservation if : 
# - hsp has a coverage between 80 an 120% compared with query length.
# - also if the distance with an other hsp is less than one thousand base pairs and the distance between two extremities of query is included between 0 and 150.
sub data_presort(){
    my ($ref_HSE, $ref_HS) = @_;
    my %hsp_start_end = %$ref_HSE;
    my %hsp_strand = %$ref_HS;
    my $old_start_hit;
    my $end_hit_interest;
    my $start_query_check;
    my $end_query_check;
    my $old_hsp_identical_length;
    my %pre_interest_hsp;
    my %scaffold_hsp_number;
    my $query_number_presort = 0;
    
    foreach my $query_name (keys %hsp_start_end){
	my $rejected_hit = 0;
	foreach my $scaffold (keys %{$hsp_start_end{$query_name}}){
	    my $i = 0;
	    my $no_old = 0;
	    my $hsp_number = 0;
	    for my $hsp (keys %{$hsp_start_end{$query_name}{$scaffold}}){
		$hsp_number ++;
	    }
	    my $new_hsp_number = 0;
	    foreach my $start_hit (sort {$a <=> $b} keys %{$hsp_start_end{$query_name}{$scaffold}}){
		foreach my $end_hit (keys %{$hsp_start_end{$query_name}{$scaffold}{$start_hit}}){
		    foreach my $start_query (keys %{$hsp_start_end{$query_name}{$scaffold}{$start_hit}{$end_hit}}){
			foreach my $end_query (keys %{$hsp_start_end{$query_name}{$scaffold}{$start_hit}{$end_hit}{$start_query}}){
			    foreach my $hsp_identical_length (keys %{$hsp_start_end{$query_name}{$scaffold}{$start_hit}{$end_hit}{$start_query}{$end_query}}){
				my $hsp_length  = $hsp_start_end{$query_name}{$scaffold}{$start_hit}{$end_hit}{$start_query}{$end_query}{$hsp_identical_length}->length('total');
				my $hsp = $hsp_start_end{$query_name}{$scaffold}{$start_hit}{$end_hit}{$start_query}{$end_query}{$hsp_identical_length};
				if ($i == 0){
				    my $query_length = $result_query_length{$query_name};
				    my $hsp_query_covering = $hsp_length / $query_length;
				    if ($hsp_number == 1){
					if (($hsp_query_covering > 0.6) && ($hsp_query_covering < 1.4)){
					    $pre_interest_hsp{$query_name}{$scaffold}{$start_hit}{$end_hit} = {
						hsp_identical_length => $hsp_identical_length,
						hsp                  => $hsp
						};
					    $new_hsp_number ++;
					    $rejected_hit ++;
					}
				    }
				    else{
					if (($hsp_query_covering > 0.6) && ($hsp_query_covering < 1.4)){
					    $pre_interest_hsp{$query_name}{$scaffold}{$start_hit}{$end_hit} = {
						hsp_identical_length => $hsp_identical_length,
						hsp                  => $hsp
						};#= $hsp_identical_length;
					    $new_hsp_number ++;
					    $rejected_hit ++;
					    $no_old ++;
					    $end_hit_interest = $end_hit + 1000;
					    $start_query_check = $start_query;
					    $end_query_check = $end_query;
					    $i ++;
					}
					else{
					    $old_start_hit = $start_hit;
					    $end_hit_interest = $end_hit + 1000;
					    $start_query_check = $start_query;
					    $end_query_check = $end_query;
					    $old_hsp_identical_length = $hsp_identical_length;
					    $i ++;
					}
				    }
				}
				else{
				    my $old_end = $end_hit_interest - 1000;
				    my $end_start_query_difference;
				    my $hsp_check;
				    if ($hsp_strand{$query_name}{$scaffold}{$start_hit} == 1){
					$end_start_query_difference = $start_query - $end_query_check;
					$hsp_check = $start_hit - $old_end;
				    }
				    else{
					$end_start_query_difference = $start_query_check - $end_query;
					$hsp_check = $start_hit - $old_end;
				    }
				    my $query_length = $result_query_length{$query_name};
				    my $hsp_query_covering = $hsp_length / $query_length;
				    if (($start_hit < $end_hit_interest) && ($end_start_query_difference < 500)){
					if (($hsp_check > 0) && ($end_start_query_difference > -40)){
					    if ($no_old == 0){
						$pre_interest_hsp{$query_name}{$scaffold}{$old_start_hit}{$old_end} = {
						hsp_identical_length => $old_hsp_identical_length,
						hsp                  => $hsp
						};#= $old_hsp_identical_length;
						$pre_interest_hsp{$query_name}{$scaffold}{$start_hit}{$end_hit} = {
						hsp_identical_length => $hsp_identical_length,
						hsp                  => $hsp
						};#= $hsp_identical_length;
						$new_hsp_number += 2;
						$rejected_hit ++;
					    }
					    else{
						$pre_interest_hsp{$query_name}{$scaffold}{$start_hit}{$end_hit}= {
						hsp_identical_length => $hsp_identical_length,
						hsp                  => $hsp
						};# = $hsp_identical_length;
						$new_hsp_number ++;
						$rejected_hit ++;
					    }
					    $end_hit_interest = $end_hit + 1000;
					    $start_query_check = $start_query;
					    $end_query_check = $end_query;
					    $no_old ++;
					}
									}
				    elsif (($hsp_query_covering > 0.6) && ($hsp_query_covering < 1.4)){
					$pre_interest_hsp{$query_name}{$scaffold}{$start_hit}{$end_hit} = {
						hsp_identical_length => $hsp_identical_length,
						hsp                  => $hsp
						};#= $hsp_identical_length;
					$new_hsp_number ++;
					$rejected_hit ++;
					$end_hit_interest = $end_hit + 1000;
					$start_query_check = $start_query;
					$end_query_check = $end_query;
					$no_old ++;
				    }
				    else{
					$old_start_hit = $start_hit;
					$end_hit_interest = $end_hit + 1000;
					$start_query_check = $start_query;
					$end_query_check = $end_query;
					$old_hsp_identical_length = $hsp_identical_length;
					$no_old = 0;
				    }
				}
			    }
			}	
		    }
		}
	    }
	    $scaffold_hsp_number{$query_name}{$scaffold} = $new_hsp_number;
	}
	if ($rejected_hit == 0){
	    # print (BMNR "$query_name\n");
	    $query_number_presort ++;
	}
    }
    my $ref_PIH = \%pre_interest_hsp;
    my $ref_SHN = \%scaffold_hsp_number;
    my @data_returned = ($ref_PIH, $ref_SHN, $query_number_presort);
    return @data_returned;
}

#Final sorting, conservation and writing in BLAT_results_analyzer.txt if :
# - strand orientation is the same with the various hsp of the same marker.
# - hsp has a coverage between 80 and 120% and an identity percentage higher than the desired percentage, all compared with query length.
#Else writing in BLAT_markers_no_results.txt for indicating markers who don't have a significant hit.
sub final_sorting(){
	my ($ref_PIH, $ref_SHN, $ref_HS) = @_;
	my %pre_interest_hsp = %$ref_PIH;
	my %scaffold_hsp_number = %$ref_SHN;
	my %hsp_strand = %$ref_HS;
	my $start_conserved;
	my $end_conserved;
	my $hsp_strand;
	my $number_query_without_significant_hit = 0;
	my $query_with_many_scaffold = 0;
	my %no_duplicate_marker_many_scaffold;
	my %markers_infos;

	foreach my $query_name (keys %pre_interest_hsp){
		my $no_significant_hit = 0;
		my $scaffold_number_per_marker = 0;
		my $previous_scaffold = "";
		foreach my $scaffold (keys %{$pre_interest_hsp{$query_name}}){
			my $i = 0;
			my $j = 0;
			my $hsp_number_per_marker_localisation = 0;
			my $scaffold_length = $result_hit_length{$scaffold};
			my $global_hsp_identical_length;
			foreach my $start (sort {$a <=> $b} keys %{$pre_interest_hsp{$query_name}{$scaffold}}){
				foreach my $end (keys %{$pre_interest_hsp{$query_name}{$scaffold}{$start}}){
					$j ++;
					$hsp_number_per_marker_localisation ++;
					my $keep_hsp = $pre_interest_hsp{$query_name}{$scaffold}{$start}{$end}->{hsp};
					my $new_hsp_number = $scaffold_hsp_number{$query_name}{$scaffold};
					if ($i == 0){
						$start_conserved = $start;
						$end_conserved = $end;
						if ($new_hsp_number == 1){
							my $position_hsp_length = $end_conserved - $start_conserved + 1;
							$global_hsp_identical_length = $pre_interest_hsp{$query_name}{$scaffold}{$start}{$end}->{hsp_identical_length};
							my $global_hsp_query_identity = $global_hsp_identical_length / $position_hsp_length * 100;
							if ($global_hsp_query_identity > $identity_percent){
								if ($previous_scaffold ne $scaffold){
									$scaffold_number_per_marker ++;
									$previous_scaffold = $scaffold;
									if ($scaffold_number_per_marker < $nombre){
										$global_hsp_query_identity = sprintf ("%.3f", $global_hsp_query_identity);
										$markers_infos{$query_name}{$scaffold}{$scaffold_length}{$hsp_number_per_marker_localisation}{$global_hsp_query_identity}{$start_conserved}{$end_conserved} = $keep_hsp;
									}
									elsif (!(exists $no_duplicate_marker_many_scaffold{$query_name})){
										delete ($markers_infos{$query_name});
										$no_duplicate_marker_many_scaffold{$query_name} = 1;
										$query_with_many_scaffold ++;
										# print (MFS "$query_name\n");
									}
								}
								elsif ($scaffold_number_per_marker < $nombre){
									$global_hsp_query_identity = sprintf ("%.3f", $global_hsp_query_identity);
									$markers_infos{$query_name}{$scaffold}{$scaffold_length}{$hsp_number_per_marker_localisation}{$global_hsp_query_identity}{$start_conserved}{$end_conserved} = $keep_hsp;
								}
								$no_significant_hit ++;
							}
						}
						else{
							$global_hsp_identical_length = $pre_interest_hsp{$query_name}{$scaffold}{$start}{$end}->{hsp_identical_length};
							$hsp_strand = $hsp_strand{$query_name}{$scaffold}{$start};
							$i ++;
						}
					}
					else{
						my $hsp_position_comparator = $start - $end_conserved;
						if (($hsp_position_comparator < 1000) && ($j == $new_hsp_number)){
							if ($hsp_strand == $hsp_strand{$query_name}{$scaffold}{$start}){
								$end_conserved = $end;
								my $query_length = $result_query_length{$query_name};
								my $position_hsp_length = $end_conserved - $start_conserved + 1;
								my $position_hsp_query_covering = $position_hsp_length / $query_length;
								if (($position_hsp_query_covering > 0.6) && ($position_hsp_query_covering < 1.4)){
									$global_hsp_identical_length += $pre_interest_hsp{$query_name}{$scaffold}{$start}{$end}->{hsp_identical_length};
									my $global_hsp_query_identity = $global_hsp_identical_length / $position_hsp_length * 100;
									if ($global_hsp_query_identity > $identity_percent){
										if ($previous_scaffold ne $scaffold){
											$scaffold_number_per_marker ++;
											$previous_scaffold = $scaffold;
											if ($scaffold_number_per_marker < $nombre){
												$global_hsp_query_identity = sprintf ("%.3f", $global_hsp_query_identity);
												$markers_infos{$query_name}{$scaffold}{$scaffold_length}{$hsp_number_per_marker_localisation}{$global_hsp_query_identity}{$start_conserved}{$end_conserved} = $keep_hsp;
											}
											elsif (!(exists $no_duplicate_marker_many_scaffold{$query_name})){
												delete ($markers_infos{$query_name});
												$no_duplicate_marker_many_scaffold{$query_name} = 1;
												$query_with_many_scaffold ++;
												# print (MFS "$query_name\n");
											}
										}
										elsif ($scaffold_number_per_marker < $nombre){
											$global_hsp_query_identity = sprintf ("%.3f", $global_hsp_query_identity);
											$markers_infos{$query_name}{$scaffold}{$scaffold_length}{$hsp_number_per_marker_localisation}{$global_hsp_query_identity}{$start_conserved}{$end_conserved} = $keep_hsp;
										}
										$no_significant_hit ++;
									}
								}
							}
						}
						elsif ($hsp_position_comparator < 1000){
							if ($hsp_strand == $hsp_strand{$query_name}{$scaffold}{$start}){
								$end_conserved = $end;
								$global_hsp_identical_length += $pre_interest_hsp{$query_name}{$scaffold}{$start}{$end}->{hsp_identical_length};
							}
						}
						elsif (($hsp_position_comparator > 1000) && ($j == $new_hsp_number)){
							my $query_length = $result_query_length{$query_name};
							my $position_hsp_length = $end_conserved - $start_conserved + 1;
							my $position_hsp_query_covering = $position_hsp_length / $query_length;
							$hsp_number_per_marker_localisation = $hsp_number_per_marker_localisation - 1;
							if (($position_hsp_query_covering > 0.6) && ($position_hsp_query_covering < 1.4)){
								my $global_hsp_query_identity = $global_hsp_identical_length / $position_hsp_length * 100;
								if ($global_hsp_query_identity > $identity_percent){
									if ($previous_scaffold ne $scaffold){
										$scaffold_number_per_marker ++;
										$previous_scaffold = $scaffold;
										if ($scaffold_number_per_marker < $nombre){
											$global_hsp_query_identity = sprintf ("%.3f", $global_hsp_query_identity);
											$markers_infos{$query_name}{$scaffold}{$scaffold_length}{$hsp_number_per_marker_localisation}{$global_hsp_query_identity}{$start_conserved}{$end_conserved} = $keep_hsp;
										}
										elsif (!(exists $no_duplicate_marker_many_scaffold{$query_name})){
											delete ($markers_infos{$query_name});
											$no_duplicate_marker_many_scaffold{$query_name} = 1;
											$query_with_many_scaffold ++;
											# print (MFS "$query_name\n");
										}
									}
									elsif ($scaffold_number_per_marker < $nombre){
										$global_hsp_query_identity = sprintf ("%.3f", $global_hsp_query_identity);
										$markers_infos{$query_name}{$scaffold}{$scaffold_length}{$hsp_number_per_marker_localisation}{$global_hsp_query_identity}{$start_conserved}{$end_conserved} = $keep_hsp;
									}
									$no_significant_hit ++;
								}
							}
							$hsp_number_per_marker_localisation = 1;
							$start_conserved = $start;
							$end_conserved = $end;
							$position_hsp_length = $end_conserved - $start_conserved + 1;
							$position_hsp_query_covering = $position_hsp_length / $query_length;
							if (($position_hsp_query_covering > 0.6) && ($position_hsp_query_covering < 1.4)){
								$global_hsp_identical_length = $pre_interest_hsp{$query_name}{$scaffold}{$start}{$end}->{hsp_identical_length};
								my $global_hsp_query_identity = $global_hsp_identical_length / $position_hsp_length * 100;
								if ($global_hsp_query_identity > $identity_percent){
									if ($previous_scaffold ne $scaffold){
										$scaffold_number_per_marker ++;
										$previous_scaffold = $scaffold;
										if ($scaffold_number_per_marker < $nombre){
											$global_hsp_query_identity = sprintf ("%.3f", $global_hsp_query_identity);
											$markers_infos{$query_name}{$scaffold}{$scaffold_length}{$hsp_number_per_marker_localisation}{$global_hsp_query_identity}{$start_conserved}{$end_conserved} = $keep_hsp;
										}
										elsif (!(exists $no_duplicate_marker_many_scaffold{$query_name})){
											delete ($markers_infos{$query_name});
											$no_duplicate_marker_many_scaffold{$query_name} = 1;
											$query_with_many_scaffold ++;
											# print (MFS "$query_name\n");
										}
									}
									elsif ($scaffold_number_per_marker < $nombre){
										$global_hsp_query_identity = sprintf ("%.3f", $global_hsp_query_identity);
										$markers_infos{$query_name}{$scaffold}{$scaffold_length}{$hsp_number_per_marker_localisation}{$global_hsp_query_identity}{$start_conserved}{$end_conserved} = $keep_hsp;
									}
									$no_significant_hit ++;
								}
							}
						}
						elsif ($hsp_position_comparator > 1000){
							my $query_length = $result_query_length{$query_name};
							my $position_hsp_length = $end_conserved - $start_conserved + 1;
							my $position_hsp_query_covering = $position_hsp_length / $query_length;
							$hsp_number_per_marker_localisation = $hsp_number_per_marker_localisation - 1;
							if (($position_hsp_query_covering > 0.6) && ($position_hsp_query_covering < 1.4)){
								my $global_hsp_query_identity = $global_hsp_identical_length / $position_hsp_length * 100;
								if ($global_hsp_query_identity > $identity_percent){
									if ($previous_scaffold ne $scaffold){
										$scaffold_number_per_marker ++;
										$previous_scaffold = $scaffold;
										if ($scaffold_number_per_marker < $nombre){
											$global_hsp_query_identity = sprintf ("%.3f", $global_hsp_query_identity);
											$markers_infos{$query_name}{$scaffold}{$scaffold_length}{$hsp_number_per_marker_localisation}{$global_hsp_query_identity}{$start_conserved}{$end_conserved} = $keep_hsp;
										}
										elsif (!(exists $no_duplicate_marker_many_scaffold{$query_name})){
											delete ($markers_infos{$query_name});
											$no_duplicate_marker_many_scaffold{$query_name} = 1;
											$query_with_many_scaffold ++;
											# print (MFS "$query_name\n");
										}
									}
									elsif ($scaffold_number_per_marker < $nombre){
										$global_hsp_query_identity = sprintf ("%.3f", $global_hsp_query_identity);
										$markers_infos{$query_name}{$scaffold}{$scaffold_length}{$hsp_number_per_marker_localisation}{$global_hsp_query_identity}{$start_conserved}{$end_conserved} = $keep_hsp;
									}
									$no_significant_hit ++;
								}
							}
							$hsp_number_per_marker_localisation = 1;
							$start_conserved = $start;
							$end_conserved = $end;
							$hsp_strand = $hsp_strand{$query_name}{$scaffold}{$start};
							$global_hsp_identical_length = $pre_interest_hsp{$query_name}{$scaffold}{$start}{$end}->{hsp_identical_length};
						}
					}
				}
			}
		}
		if ($no_significant_hit == 0){
			# print (BMNR "$query_name\n");
			$number_query_without_significant_hit ++;
		}
	}
	my $ref_MI = \%markers_infos;
	my @final_sorting_result = ($ref_MI, $query_with_many_scaffold, $number_query_without_significant_hit);
	return @final_sorting_result;
}

#Philippe FRANCOIS (CIRAD) - philippe.francois@cirad.fr
#Version [3.0.0]
#Ver. 1.0.1 : Fixed a bug in scaffold's size recovery.
#Ver. 1.1.0 : Rewriting analysis core of BLAST results for more specific results.
#Ver. 1.1.1 : Fixed a bug who count 2x global_hsp_identical_length when new_hsp_number equals 1 and added display of markers who don't pass the data presort phase.
#Ver. 1.1.2 : Fixed a bug that consider the hsp who overlap with the good hsp, and a bug with display of hsp number.
#Ver. 1.1.3 : Markers who matches on five or more scaffolds are written in a new file, markers who are eliminated at presort data are written in a different new file, display of markers with "No results found" and modification of identity calculation methodology for better results.
#Ver. 3.0.0 : Removing the creation of the 3 BRA* file, addition of a new parametre --nombre | Modif Guillaume MARTIN - guillaume.martin@cirad.fr
#Date [16/04/14]
