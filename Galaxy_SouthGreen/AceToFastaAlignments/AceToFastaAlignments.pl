#!/usr/bin/perl

use strict;
use Switch;
use Getopt::Long;

my $usage = qq~Usage:$0 <args> [<opts>]
where <args> are:
    -i, --input_ACE         <ACE input file>
    -m, --minimum_depth     <minimum depth to call a position>
    -h, --hetero_depth      <minimum depth to estimate heterozygosity>
    -o, --occ_min           <minimum allele occurence to call heterozygous>
    -f, --freq_min          <minimum allele frequency to call heterozygous>
    -p, --project_name      <project_name>
    
    <opts> are:
    -r, --remove_gaps       <remove positions with gaps: yes/no. Default: no>
~;
$usage .= "\n";

my ($ace_file,$minimum_depth,$minimum_depth_for_heterozygosity,$project,$minimal_allele_occurence,$minimum_allele_frequency,$remove_position_with_gaps);


GetOptions(
	"input_ACE=s"      => \$ace_file,
	"minimum_depth=s"  => \$minimum_depth,
	"hetero_depth=s"   => \$minimum_depth_for_heterozygosity,
	"occ_min=s"        => \$minimal_allele_occurence,
	"freq_min=s"       => \$minimum_allele_frequency,
	"project_name=s"   => \$project,
	"remove_gaps=s"    => \$remove_position_with_gaps
);


die $usage
  if ( !$ace_file || !$minimum_depth || !$project  || !$minimal_allele_occurence  || !$minimum_allele_frequency);
  


my %contigs;

mkdir($project);
	
##########################################################
# instanciation of hash table containing reads
##########################################################
open(ACE,$ace_file);
my $contig_name;
my $read_name;
my $go_for_contig = 0;
my $go_for_read_sequence = 0;
my %consensi;
my %positions;
while(<ACE>)
{
	if ($go_for_read_sequence && /^QA\s*\d+\s+\d+/)
	{
		$go_for_read_sequence = 0;
	}
	if ($go_for_contig && $contig_name && /[ATGCN\*]+/)
	{
		my $line = $_;
		$line =~s/\n//g;
		$line =~s/\r//g;
		if ($consensi{$contig_name})
		{
			$consensi{$contig_name} .= $line;
		}
		else
		{
			$consensi{$contig_name} = $line;
		}
	}
	if ($go_for_read_sequence && $contig_name && /[AaTtGgCcNn\*]+/)
	{
		my $line = $_;
		$line =~s/\n//g;
		$line =~s/\r//g;
		if ($contigs{$contig_name}{$read_name})
		{
			$contigs{$contig_name}{$read_name} .= $line;
		}
		else
		{
			$contigs{$contig_name}{$read_name} = $line;
		}

	}
	if ($go_for_contig && /^BQ\s*$/)
	{
		$go_for_contig = 0;
	}
		
	if (/^CO\s+([\w-\|\.]+)\s+\d+\s+\d+\s+\d+\s+/)
	{
		$contig_name = $1;
		$contig_name =~s/_//g;
		$go_for_contig = 1;
	}
	if (/^AF\s+([\w-\|\.]+)\s+[UC]+\s+(-*\d+)\s+$/)
	{
		my $read = $1;
		my $pos = $2;
		$positions{$contig_name}{$read} = $pos;
	}

	
	if (/^RD\s+([\w-\|\.]+)\s+(\d+)\s+\d+\s+\d+\s*$/)
	{
		$go_for_read_sequence = 1;
		$read_name = $1;
		my $read_length = $2;
	}
	if (/^QA\s+(\d+)\s+(\d+)/)
	{
		my $start = $1;
		my $end = $2;
			
		# change position according to quality limits
		my $position = $positions{$contig_name}{$read_name};
		my $new_position = $position + $start - 1;
		$positions{$contig_name}{$read_name} = $new_position;
			
		# change sequence according to quality limits
		my $sequence = $contigs{$contig_name}{$read_name};
		
		my $length = $end - $start + 1;
		my $new_sequence = substr($sequence,($start - 1),$length);
		$new_sequence =~s/\*/-/g;
		$contigs{$contig_name}{$read_name} = $new_sequence;
	}
}
close(ACE);

	
############################################################
# estimate heterozygosity and creation of a fasta alignment
############################################################
my $depth_for_heterozygosity_estimation = 6;
my $minimum_allele_frequency = 0.2;
my $remove_position_with_gaps = "yes";
my %sum_positions;
my %sum_heterozygous_positions;
open(HET,">$project.heterozygous_sites");
	
	
if (-e "$project.SNP")
{
	unlink("$project.SNP");
	unlink("$project.contigs_fasta_alignments.zip");
}

my $nb_contig = 0;
foreach my $contig(keys(%contigs))
{
	#############################################
	# creation of a fasta alignment
	#############################################
	my $ref = $contigs{$contig};
	my %reads = %$ref;
	my $consensus = $consensi{$contig};
	my $consensus_length = length($consensus);
	
	$nb_contig++;
	open(FAS,">$project/$contig.reads.alignment.fas");
	foreach my $read(keys(%reads))
	{
		my $seq = $contigs{$contig}{$read};

		my $position = $positions{$contig}{$read};
		
		my $sequence = "";
		for (my $j=0; $j < ($position - 1); $j++)
		{
			$sequence .= "?";
		}
		$sequence .= $seq;
			
		my $seq_length = length($sequence);
		if ($seq_length < $consensus_length)
		{
			my $ajout = "?" x ($consensus_length - $seq_length);
			$sequence = $sequence . $ajout;
		}
		$sequence =~s/\n//g;
		$sequence =~s/\r//g;
		print FAS ">$read\n$sequence\n";
	}
	close(FAS);
	
	
	
	###############################
	# estimate heterozygosity 
	###############################
		
	my $result = estimateHeterozygosity($contig);
	my @lines_result = split(/\n/,$result);
	

	my $nb_heterozygous_positions = 0;
	my $nb_positions = 0;
	foreach my $line(@lines_result)
	{
		my @infos = split(/\t/,$line);
		my $ind = $infos[0];
		my $nb_heterozygous_positions = $infos[1];
		my $nb_positions = $infos[2];
		if ($sum_positions{$ind})
		{
			$sum_positions{$ind} += $nb_positions;
		}
		else
		{
			$sum_positions{$ind} = $nb_positions;
		}
		if ($sum_heterozygous_positions{$ind})
		{
			$sum_heterozygous_positions{$ind} += $nb_heterozygous_positions;
		}
		else
		{
			$sum_heterozygous_positions{$ind} = $nb_heterozygous_positions;
		}
		if ($nb_heterozygous_positions > 0)
		{
			print HET "$ind	$contig	$nb_heterozygous_positions	$nb_positions\n";
		}
	}
	

	#################
	# cleaning
	#################
	if ($contig)
	{
		system("rm -rf $project/$contig.reads*");
	}
}


print HET "\n\nFINAL RESULTS:\n\n";
foreach my $ind(keys(%sum_heterozygous_positions))
{
	print HET "$ind:\n";
	my $portion = 0;
	if ($sum_heterozygous_positions{$ind})
	{
		$portion = $sum_positions{$ind} / $sum_heterozygous_positions{$ind};
	}
	$portion = sprintf("%.1f", $portion);
	
	print HET $sum_heterozygous_positions{$ind} . " / " . $sum_positions{$ind} . "\n";
	print HET "==> ~ 1 heterozygous position every $portion pb\n\n";
}
close(HET);

#chdir($project);

#########################################################
# reorganize alignments in different zip files
#########################################################

use File::Copy;
opendir (DIR, "$project" ) || die "Error in opening dir $project\n";
my $count = 0;
my $count_dir = 0;
my $set_dir = "set0";
mkdir("$project/$set_dir");
my $nb_contig2 = 0;
while( (my $fic = readdir(DIR)))
{
	if ($fic =~/align.fas/)
	{
		$count++;
		copy("$project/$fic","$project/$set_dir");
	}
	if ($fic =~/depth/)
	{
		copy("$project/$fic","$project/$set_dir");
	}
	if ($count == 8000)
	{
		$count_dir++;
		$set_dir = "set" . $count_dir;
		mkdir("$project/$set_dir");
		$count = 0;
	}
	$nb_contig2++;
}
close(DIR);

opendir (DIR2, "$project" ) || die "Error in opening dir $project\n";
while( (my $subdir = readdir(DIR2)))
{
	if ($subdir =~/set/)
	{
		chdir("$project/$subdir");
		system("zip ../../$project.$subdir.fasta_alignments.zip ./*.align.fas >$project.transfer.log");
		system("zip ../../$project.$subdir.depths.zip ./*.depth >$project.transfer.log");
		chdir("../..");
	}
}
close(DIR2);


#system("zip ../$project.contigs_fasta_alignments.zip ./*.align.fas >$project.transfer.log");
#chdir("..");

if ($project && -e $project)
{
	system("rm -rf ./$project");
}


sub estimateHeterozygosity($)
{
	my $contig = $_[0];
	my $result;
	my %IUPAC = 
	(
		'[A/G]' => "R",
		'[G/A]' => "R",
		'[C/T]' => "Y",
		'[T/C]' => "Y",
		'[G/T]' => "K",
		'[T/G]' => "K",
		'[C/G]' => "S",
		'[G/C]' => "S",
		'[A/T]' => "W",
		'[T/A]' => "W",
		'[A/C]' => "M",
		'[C/A]' => "M",
		'[A/-]' => "L",
		'[C/-]' => "L",
		'[G/-]' => "L",
		'[T/-]' => "L",
		'[-/A]' => "L",
		'[-/C]' => "L",
		'[-/G]' => "L",
		'[-/T]' => "L"
	);


	my $seq;
	my $seq_name;
	my $ind;
	my %sequences;
	open(FASTA,"$project/$contig.reads.alignment.fas");
	while(<FASTA>)
	{
		if (/>(\w+)/)
		{
			if ($seq_name && $seq && $ind)
			{
				$sequences{$ind}{$seq_name} = $seq;
			}
			$seq_name = $1;
			if ($seq_name =~/([^_]+)_(.*)/)
			{
				$ind = $1;
				$seq_name = $2;
			}
			#my @infos = split("_",$seq_name);
			#$ind = $infos[0];
			#$seq_name = $infos[1];
			$seq = "";
		}
		else
		{
			my $line = $_;
			chomp($line);
			$seq .= $line;
		}
	}
	close(FASTA);
	if ($seq_name && $seq && $ind)
	{
		$sequences{$ind}{$seq_name} = $seq;
	}
	
	open(SNP,">>$project.SNP");
	open(NEW_ALIGN,">$project/$contig.align.fas");
	open(DEPTH,">$project/$contig.depth");
	foreach my $ind(keys(%sequences))
	{
		my $ref = $sequences{$ind};
		my %seqs = %$ref;
		print NEW_ALIGN "\n>$ind" . "_" . $contig . "\n";
		print DEPTH "\n>$ind" . "_" . $contig . "\n";
		
		
		###############################################################
		# get nucleotids for each position
		###############################################################
		my %positions;
		my $homopolymere_position;
		
		foreach my $key(keys(%seqs))
		{
			my $sequence = $seqs{$key};

			################################
			# mask homopolymeres
			################################
			## TO DO	
	
			#$sequence =~s/N/-/g;
			my @nts = split(//,$sequence);
			for (my $i = 0; $i <= $#nts; $i++)
			{
				my $nucleotid = $nts[$i];
				if ($positions{$i}{$nucleotid})
				{
					$positions{$i}{$nucleotid}++;
				}
				else
				{
					$positions{$i}{$nucleotid} = 1;
				}	
			}
		}
		
		
		my $nb_snp = 0;
		my $nb_positions_to_consider_for_heterozygosity = 0;
		my $nb_positions_to_keep = 0;
		my @sorted_positions = sort{$a<=>$b}keys(%positions);
		foreach my $pos(@sorted_positions)
		{
			my $ref = $positions{$pos};
			$pos++;
			my %nts = %$ref;
	
			my $nb_total = 0;
			my $nb_missing = 0;
			foreach my $nt(keys(%nts))
			{
				$nb_total += $nts{$nt};
				if ($nt eq '?')
				{
					$nb_missing = $nts{$nt};
				}
			}
			my $depth = $nb_total - $nb_missing;
			
			if ($depth >= $minimum_depth_for_heterozygosity)
			{
				$nb_positions_to_consider_for_heterozygosity++;
			}
			
			

			# if insufficient depth, add missing data for the individual
			if ($depth < $minimum_depth)
			{
				%nts = ();
				$nts{"?"} = 1;
			}

			
			print DEPTH "$depth ";
			
			###################################
			# different letters for this column
			###################################
			if (scalar keys(%nts) > 1)
			{
				my @alleles_for_snp;
				my @occurences_for_snp;
				my @whole_nt;
				my @occurences_for_whole;
				my %alleles2;
				my $unique_letter;
						
				#########################
				# get alleles
				#########################
				foreach my $nt(keys(%nts))
				{
					if ($nt ne '?' && $nt ne 'N' && $nt ne '-')
					{
						push(@alleles_for_snp,$nt);
						push(@occurences_for_snp,$nts{$nt});
						$alleles2{$nts{$nt}} = $nt;
					}
					if ($nt ne '?')
					{
						push(@whole_nt,$nt);
						push(@occurences_for_whole,$nts{$nt});
					}
				}
				
				###############################
				# manage occurence and frequency
				###############################
				my @sorted_occurences_for_snp = sort{$a<=>$b}@occurences_for_snp;
				my @sorted_occurences_for_whole = sort{$a<=>$b}@occurences_for_whole;
				my $snp_frequency = 0;
				my $occurence_snp = $sorted_occurences_for_snp[$#sorted_occurences_for_snp - 1];
				my $occurence_majoritaire = $sorted_occurences_for_snp[$#sorted_occurences_for_snp];
				if ($depth && $occurence_snp)
				{
					$snp_frequency = $occurence_snp / $depth;
				}
				
				##############################
				# heterozygous position found 
				##############################
				if (scalar @alleles_for_snp > 1)
				{
					####################################################################
					# sufficient depth, allele frequency, and minority allele occurence
					####################################################################
					if ($occurence_snp >= $minimal_allele_occurence && $snp_frequency >= $minimum_allele_frequency && $depth >= $minimum_depth_for_heterozygosity)
					{
						my $snp;
						if (scalar @alleles_for_snp == 2)
						{
							@alleles_for_snp = sort(@alleles_for_snp);
							$snp = "[" . join("/",@alleles_for_snp). "]";
						}
						else
						{
							$snp = "[" . $alleles2{$occurence_majoritaire} . "/" . $alleles2{$occurence_snp} . "]";

							# if the two majoritary alleles show the same abundance (example: A:3, G:6, C:6)
							if ($occurence_majoritaire == $occurence_snp)
							{
								my @new_alleles_for_snp;
								foreach my $nt(keys(%nts))
                                				{
				                                        if ($nt ne '?' && $nt ne 'N' && $nt ne '-' && $nts{$nt} == $occurence_majoritaire)
                                				        {
										push(@new_alleles_for_snp,$nt);
									}
								}
								$snp = "[" . join("/",@new_alleles_for_snp). "]";
							}
						}
						 
						#print "$pos $depth $snp $frequency\n";
						if ($pos && $snp)
						{
							print SNP "$contig	$ind	$snp	$pos\n";
						}
						
						$nb_snp++;
						$snp=uc($snp);
						if ($IUPAC{$snp})
						{
							print NEW_ALIGN $IUPAC{$snp};
						}
						else
						{
							print NEW_ALIGN '?';
						}
					}
					#####################################
					# add majority allele in consensus
					#####################################
					else
					{
						my $nt_for_consensus;
							
						# same number of A and G for example
						if ($sorted_occurences_for_whole[$#sorted_occurences_for_whole] == $sorted_occurences_for_whole[$#sorted_occurences_for_whole - 1])
						{
							$nt_for_consensus = "?";
						}
						# add majority allele
						else
						{
							$nt_for_consensus = $alleles2{$sorted_occurences_for_whole[$#sorted_occurences_for_whole]};
						}
						if ($nt_for_consensus)
						{
							print NEW_ALIGN $nt_for_consensus;
						}
						else
						{
							print NEW_ALIGN '?';
						}
					}
				}
				################################
				# C/- A/N 
				################################
				elsif (scalar @whole_nt > 1)
				{
					@whole_nt = sort(@whole_nt);
					my $variation = "[" . join("/",@whole_nt). "]";
					
					##########################
					# if indel, add majority
					##########################
					if ($variation =~/\-/)
					{
						my $nt_for_consensus = $alleles2{$sorted_occurences_for_whole[$#sorted_occurences_for_whole]};
						if (!$nt_for_consensus)
						{
							$nt_for_consensus = "-";
						}
						print NEW_ALIGN $nt_for_consensus;
					}
					##########################
					# if A/N, add A
					##########################
					else
					{
						my $ok = 0;
						foreach my $nt(keys(%nts))
						{
							if ($nt =~/[ATGC]/)
							{
								print NEW_ALIGN $nt;
								$ok = 1;
							}	
						}
						if (!$ok)
						{
							print NEW_ALIGN '?';
						}
					}
				}
				##############################
				# only missing data and letter
				##############################
				else
				{
					my $ok = 0;
					foreach my $nt(keys(%nts))
					{
						if ($nt ne '?')
						{
							print NEW_ALIGN $nt;
							$ok = 1;
						}	
					}
					if ($ok == 0)
					{
						print NEW_ALIGN '?';
					}
				}
			}
			###################################
			# no difference (only one letter)
			###################################
			else
			{
				foreach my $nt(keys(%nts))
				{
					print NEW_ALIGN $nt;
				}
			}
		}
		
		$result .= "$ind	$nb_snp	$nb_positions_to_consider_for_heterozygosity\n";
	}
	close(NEW_ALIGN);
	close(SNP);
	close(DEPTH);
	
	
	################################################################
	# remove positions with gaps
	################################################################
	
	if ($remove_position_with_gaps eq "yes")
	{
		my %sequences2;
		my $ind2;
		my $seq_name2;
		my $seq2;
		
		open(O,">$project/$contig.align.fas.2");
		open(I,"$project/$contig.align.fas");
		while(<I>)
		{
			if (/>(.+)\s*$/)
			{
				if ($seq_name2 && $seq2)
				{
					$sequences2{$seq_name2} = $seq2;
				}
				$seq_name2 = $1;
				$seq2 = "";
			}
			else
			{
				my $line = $_;
				chomp($line);
				$seq2 .= $line;
			}			
		}
		close(I);
		if ($seq_name2 && $seq2)
		{
			$sequences2{$seq_name2} = $seq2;
		}
		
		my %depths2;
		$ind2 = "";
		$seq_name2 = "";
		$seq2 = "";
		open(D,"$project/$contig.depth");
		open(D2,">$project/$contig.depth.2");
                while(<D>)
                {
                        if (/>(.+)\s*$/)
                        {
                                if ($seq_name2 && $seq2)
                                {
                                        $depths2{$seq_name2} = $seq2;
                                }
                                $seq_name2 = $1;
                                $seq2 = "";
                        }
                        else
                        {
                                my $line = $_;
                                chomp($line);
                                $seq2 .= $line;
                        }
                }
                close(D);
                if ($seq_name2 && $seq2)
                {
                        $depths2{$seq_name2} = $seq2;
                }
		
		my %positions2;
		my %positions2_for_depth;
		my %positions_to_remove;
		foreach my $ind(keys(%sequences2))
		{
			my $sequence = $sequences2{$ind};
			my $depth = $depths2{$ind};
			my @nts = split(//,$sequence);
			for (my $i = 0; $i <= $#nts; $i++)
			{
				my $nucleotid = $nts[$i];
				if ($nucleotid eq '-')
				{
					$positions_to_remove{$i} = 1;
				}
				$positions2{$i}{$ind} = $nucleotid;
			}
			
			my @depth_values = split(/ /,$depth);
			for (my $i = 0; $i <= $#depth_values; $i++)
                        {
				my $depth_value = $depth_values[$i];
				$positions2_for_depth{$i}{$ind} = $depth_value;
			}
		}
		
		
		foreach my $ind(keys(%sequences2))
		{
			my $sequence = $sequences2{$ind};
			print O ">$ind\n";
			print D2 ">$ind\n";
			my @sorted_positions = sort{$a<=>$b}keys(%positions2);
			foreach my $pos(@sorted_positions)
			{
				if (!$positions_to_remove{$pos})
				{
					print O $positions2{$pos}{$ind};
					print D2 $positions2_for_depth{$pos}{$ind} . " ";
				}
			}
			print O "\n";
			print D2 "\n";
		}
		close(O);
		close(D2);
		
		rename("$project/$contig.align.fas.2","$project/$contig.align.fas");
		rename("$project/$contig.depth.2","$project/$contig.depth");
	}
	
	return $result;
}
