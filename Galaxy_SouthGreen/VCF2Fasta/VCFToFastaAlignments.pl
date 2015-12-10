#!/usr/bin/perl

use strict;
use Switch;
use Getopt::Long;

my $usage = qq~Usage:$0 <args> [<opts>]
where <args> are:
    -i, --input_VCF         <VCF input file>
    -r, --reference         <reference FASTA file>
    -d, --depth_file        <depth file>
    -p, --project_name      <project_name>
    -n, --name_of_reference <name of reference>
    
    <opts> are:
    -m, --max_nb_positions  <maximum number of positions. Default: 10000000>
    -c, --coverage_min      <minimal depth of coverage to consider position>
    -f, --frequency_min     <minimal frequency of variation to consider as heterozygote>
~;
$usage .= "\n";

my ($vcf_file,$reference_file,$depth_file,$tmpdir,$name_ref);
my $max_nb_pos = 10000000;
my $min_coverage = 1;
my $min_frequency_for_heterozygote = 30;

GetOptions(
	"input_VCF=s"         => \$vcf_file,
	"reference=s"         => \$reference_file,
	"depth_file=s"        => \$depth_file,
	"project_name=s"      => \$tmpdir,
	"name_of_reference=s" => \$name_ref,
	"max_nb_positions=s"  => \$max_nb_pos,
	"coverage_min=s"      => \$min_coverage,
	"frequency_min=s"     => \$min_frequency_for_heterozygote
);


die $usage
  if ( !$vcf_file || !$reference_file || !$tmpdir  || !$name_ref);
  

my %individus;
$individus{$name_ref} = 1;

#my $tmpdir = int(rand(10000000000000));
mkdir $tmpdir;

my %sequences;
my $seq_name = "";
my $sequence = "";
	

my %limits;
my $nb = 1;
my $nb_pos = 0;
my $count_gene = 0;
my $seqname;
open(REF,$reference_file);
while(<REF>)
{
	my $line = $_;
	chomp($line);
	
	if (/>([^\s]+)\s*/)
	{
		my $val = $1;
		#$val =~s/_/-/g;
		if ($seq_name && $sequence)
		{
			$sequences{$seq_name} = $sequence;
			
			$nb_pos += length($sequence);
			if ($nb_pos > $max_nb_pos)
			{
				$nb++;
				$nb_pos = 0;
				$count_gene = 0;
			}
			
			
		}
		$seqname = $val;
		$seq_name = $val;
		$seq_name =~s/\///g;
		$seq_name =~s/\\//g;
		$sequence = "";
		
		$count_gene++;
		$limits{$nb}{$count_gene} = $seqname;
	}
	else
	{
		$sequence .= $line;
	}
}
close(REF);
if ($seq_name && $sequence)
{
	$sequences{$seq_name} = $sequence;
}


##############################################################
# split input files
##############################################################
my $header_vcf = "";
open(V1,$vcf_file);
while(<V1>)
{
	$header_vcf .= $_;
	if (/#CHROM/)
	{
		last;
	}
}
close(V1);
chomp($header_vcf);

open(V,$vcf_file);
open(D,$depth_file);
my $last_line_vcf = "";
my $last_line_depth = "";
foreach my $nb(sort {$a <=> $b} keys(%limits))
{
	my $r = $limits{$nb};
	my %hash = %$r;

	###############################################
	# VCF
	###############################################
	
	
	open(V2,">$vcf_file.$nb");
	print V2 $header_vcf;
	print V2 "\n";
	print V2 $last_line_vcf;
	my %go;
	my $previous = 0;

	while(<V>)
	{
		my @infos = split(/\t/,$_);
		my $gene = $infos[0];
		#$gene =~s/_/-/g;
		
		if (!$go{$gene})
		{
			for (my $i=$previous;$i<=scalar keys(%hash);$i++)
			{
				my $current_gene = $hash{$i};

				if ($current_gene eq $gene)
				{
					$go{$gene} = 1;
					$previous = $i;
				}
			}
		}

		
		if (!/^#/ && $gene)
		{

			if ($go{$gene})
			{
				print V2 $_;
			}
			elsif ($previous)
			{
				$last_line_vcf = $_;
				last;
			}
			else
			{
				next;
			}
		}
	}
	
	close(V2);
	
	###############################################
	# depth
	###############################################
	my $header_depth = `head -1 $depth_file`;
	chomp($header_depth);
	
	open(D2,">$depth_file.$nb");
	print D2 $header_depth;
	print D2 "\n";
	print V2 $last_line_depth;
	my %go2;
	my $previous = 0;
	
	while(<D>)
	{
		my @infos = split(/\t/,$_);
		my @infos2 = split(/:/,$infos[0]);
		my $gene = $infos2[0];
		#$gene =~s/_/-/g;
		if (!$go2{$gene})
		{
			for (my $i=$previous;$i<=scalar keys(%hash);$i++)
			{
				my $current_gene = $hash{$i};

				if ($current_gene eq $gene)
				{
					$go2{$gene} = 1;
					$previous = $i;
				}
			}
		}

			
		if (!/^Locus/ && $gene)
		{

			if ($go2{$gene})
			{
				print D2 $_;
			}
			elsif ($previous)
			{
				$last_line_depth = $_;
				last;
			}
			else
			{
				next;
			}
		}
	}
	
	close(D2);
}
close(V);
close(D);




my %flanking;
open(VCF_LIGHT,">$tmpdir.heterozygote_positions.xls");

########################################################
# analysis for each splitted files
########################################################

open(LS,"ls $vcf_file.* |");
while(<LS>)
{
	my $file = $_;
	chomp($file);

	if ($file =~/(.*)\.(\d+)/)
	{
		my $particule = $1;
		my $num = $2;
		
		
		$vcf_file = $file;
		my $depth_file_tmp = "$depth_file.$num";
		
		
		my $gene;
		my %genes;
		my %genotyping;
		my %snp_stats;
		my %num_genes;
		my %columns;
		
		my $seqname;
		
		
		
		my $nb_snps = 0;
		my $num_gene = 0;
		my %indiv_order;
		
		open(VCF,"$vcf_file");
		while(<VCF>)
		{
			my $line = $_;
			chomp($line);
			if (/#CHROM/)
			{
				my @infos = split(/\t/,$line);
				print VCF_LIGHT "CHROM	POS	";
				for (my $j=9;$j<=$#infos;$j++)
				{
					my $individu = $infos[$j];
					print VCF_LIGHT "$individu	";
					$indiv_order{$j} = $individu;
					$individus{$individu}=1;
				}
				print VCF_LIGHT "\n";
			}
			elsif (!/^#/)
			{
				my @infos = split(/\t/,$line);
			
				my $ref_allele = $infos[3];
				my $alternate_allele = $infos[4];
				
				my $pos = $infos[1];
				$gene = $infos[0];
				my $filter = $infos[6];
				if ($filter ne 'PASS')
				{
					next;
				}
				
				$gene =~s/\///g;
				$gene =~s/\\//g;
		
				###################################
				# initialization of a new region
				###################################
				if ($gene && !$genes{$gene})
		        {
		        	$num_gene++;
			        $genes{$gene} = $num_gene;
			        $num_genes{$num_gene} = $gene;
					$nb_snps = 0;
		        }
		
		        my $snp_type = "[$ref_allele/$alternate_allele]";
				
				switch($snp_type)
				{
						case '[A/G]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "R";}
						case '[G/A]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "R";}
						case '[C/T]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "Y";}
						case '[T/C]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "Y";}
						case '[T/G]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "K";}
						case '[G/T]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "K";}
						case '[C/G]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "S";}
						case '[G/C]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "S";}
						case '[A/T]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "W";}
						case '[T/A]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "W";}
						case '[A/C]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "M";}
						case '[C/A]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "M";}
						case '[C/A/T]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "H";}
						case '[A/T/C]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "H";}
						case '[A/C/T]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "H";}
						case '[C/T/A]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "H";}
						case '[T/C/A]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "H";}
						case '[T/A/C]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "H";}
						
						case '[C/A/G]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "V";}
						case '[A/G/C]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "V";}
						case '[A/C/G]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "V";}
						case '[C/G/A]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "V";}
						case '[G/C/A]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "V";}
						case '[G/A/C]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "V";}
						
						case '[C/T/G]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "B";}
						case '[T/G/C]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "B";}
						case '[T/C/G]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "B";}
						case '[C/G/T]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "B";}
						case '[G/C/T]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "B";}
						case '[G/T/C]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "B";}
						
						case '[T/A/G]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "D";}
						case '[A/G/T]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "D";}
						case '[A/T/G]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "D";}
						case '[T/G/A]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "D";}
						case '[G/T/A]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "D";}
						case '[G/A/T]' {$snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"} = "D";}
				}
					
				$nb_snps++;
				
					
				$genotyping{$num_gene}{$pos}{$name_ref}{"genotype"} = "$ref_allele:$ref_allele";
				$genotyping{$num_gene}{$pos}{$name_ref}{"code_genotype"} = $ref_allele;
					
				my $heterozygote_content = 0;
					
				###################################
				# get genotypes for each sample
				###################################
				for (my $j=9;$j<=$#infos;$j++)
				{
					if (!$indiv_order{$j})
		            {
		            	next;
		            }
					my @fields = split(/:/,$infos[$j]);
					my $genotype = $fields[0];
					my $global_depth = $fields[2];
					my $heterozygous = 0;
					if ($genotype eq "1/0" or $genotype eq "0/1")
                                        {
						$heterozygous = 1;
					}
					$genotype =~s/0/$ref_allele/g;
					

					
					# if low coverage, do not consider the position (consider as missing data)
					if ($global_depth <= $min_coverage)
					{
						$genotype = "./.";
					}
					if ($heterozygous)
					{ 
						my ($allele1_depth,$allele2_depth) = split(",",$fields[1]);
						my $heterozygosity_frequency = 100;
                                        	if ($allele1_depth <= $allele2_depth && $global_depth > 0)
                                        	{
                                                	$heterozygosity_frequency = ($allele1_depth / $global_depth) * 100;
                                        	}
	                                        elsif ($allele1_depth > $allele2_depth && $global_depth > 0)
        	                                {
                	                                $heterozygosity_frequency = ($allele2_depth / $global_depth) * 100;
                        	                }

						# if frequency of alleles (heterozygote) lower than the threshold, consider as missing data 
						# (TODO: consider as homozygote (with majority allele))
                	                        if ($heterozygosity_frequency <= $min_frequency_for_heterozygote)
                        	                {
                                	                $genotype = "./.";
                                        	}
					}
						
					# several alternate alleles
					if ($alternate_allele =~/,/)
					{
						my @alternate_alleles = split(",",$alternate_allele);
						$genotype =~s/1/$alternate_alleles[0]/g;
						$genotype =~s/2/$alternate_alleles[1]/g;
						if ($alternate_alleles[2])
						{
							$genotype =~s/3/$alternate_alleles[2]/g;
						}
					}
					else
					{
						$genotype =~s/1/$alternate_allele/g;
					}
						
					$genotype =~s/\./0/g;
					my @alleles;
					if ($genotype =~/\//)
					{
						@alleles = split(/\//,$genotype);
					}
					elsif ($genotype =~/\|/)
		            {
		            	@alleles = split(/\|/,$genotype);
		            }
					my $individu = $indiv_order{$j};
						
					
					$genotyping{$num_gene}{$pos}{$individu}{"genotype"} = join(":",@alleles);
					$genotype = $genotyping{$num_gene}{$pos}{$individu}{"genotype"};
					my $code_genotype = "";
					
					
					switch ($genotype) 
					{
						case 'A:A' {$code_genotype = "A";}
						case 'T:T' {$code_genotype = "T";}
						case 'C:C' {$code_genotype = "C";}
						case 'G:G' {$code_genotype = "G";}
						case 'A:G' {$code_genotype = "R";}
						case 'G:A' {$code_genotype = "R";}
						case 'C:T' {$code_genotype = "Y";}
						case 'T:C' {$code_genotype = "Y";}
						case 'G:T' {$code_genotype = "K";}
						case 'T:G' {$code_genotype = "K";}
						case 'C:G' {$code_genotype = "S";}
						case 'G:C' {$code_genotype = "S";}
						case 'A:T' {$code_genotype = "W";}
						case 'T:A' {$code_genotype = "W";}
						case 'A:C' {$code_genotype = "M";}
						case 'C:A' {$code_genotype = "M";}
						case '0:0' {$code_genotype = "0";}
					}
					
					if ($code_genotype eq 'R' or $code_genotype eq 'Y' or $code_genotype eq 'K' or $code_genotype eq 'S' or $code_genotype eq 'W' or $code_genotype eq 'M')
					{
						$heterozygote_content = 1;
					}
		
					$genotyping{$num_gene}{$pos}{$individu}{"code_genotype"} = $code_genotype;
					
						
					foreach my $allele(@alleles)
					{
						if ($allele && $allele ne '0')
						{
							if ($columns{$j})
							{
								if ($columns{$j} !~ $allele)
								{
									$columns{$j} .= "/" . $allele;
								}
							}
							else
							{
								$columns{$j} = $allele;
							}
						}
					}
				}
				
				
				############################################################
				# get depth in case of heterozygosity
				############################################################
				if ($heterozygote_content)
				{
					my $new_gene = $gene;
					$new_gene =~s/_/-/g;
					print VCF_LIGHT "$new_gene	$pos	";
					for (my $j=9;$j<=$#infos;$j++)
					{
						my @fields = split(/:/,$infos[$j]);
						my $genotype = $fields[0];
						$genotype =~s/0/$ref_allele/g;
						my $depth = $fields[1];
						$depth=~s/\s//g;
							
						# several alternate alleles
						if ($alternate_allele =~/,/)
						{
							my @alternate_alleles = split(",",$alternate_allele);
							$genotype =~s/1/$alternate_alleles[0]/g;
							$genotype =~s/2/$alternate_alleles[1]/g;
							if ($alternate_alleles[2])
							{
								$genotype =~s/3/$alternate_alleles[2]/g;
							}
						}
						else
						{
							$genotype =~s/1/$alternate_allele/g;
						}
							
						$genotype =~s/\./0/g;
						print VCF_LIGHT "$genotype:$depth	";
					}
					print VCF_LIGHT "\n";
				}
			}
		}
		
		
		
		my %depth;
		my %depth_sequences;
		my %intervals;
		#####################################################
		# get depth of coverage information
		#####################################################
		if ($depth_file_tmp && -e "$depth_file_tmp")
		{
			my $line = 0;
			my %order_indiv;
			open(DEPTH,"$depth_file_tmp");
			while(<DEPTH>)
			{
				$line++;
				my @infos = split("\t",$_);
	
				if ($line == 1)
				{
					for (my $i = 3; $i <= $#infos; $i++)
					{
						if ($infos[$i] =~/Depth_for_(.*)/)
						{
							$order_indiv{$i} = $1;
						}
					}
				}
				else
				{
					my ($gene,$position) = split(":",$infos[0]);
					for (my $i = 3; $i <= $#infos; $i++)
					{
						my $individu = $order_indiv{$i};
						$gene =~s/\///g;
	                    $gene =~s/\\//g;
	                    my $current_depth = $infos[$i];
	                    $current_depth =~s/\s//g;
	#                    if ($depth_sequences{$gene}{$individu})
	#                    {
	#                    	$depth_sequences{$gene}{$individu} .= " " . $current_depth;
	#                    }
	#                    else
	#                    {
	#                    	$depth_sequences{$gene}{$individu} = $current_depth;
	#                    }
						$depth{$gene}{$individu}{$position} = $current_depth;
					}
				}
			}
			close(DEPTH);
		}
		
		foreach my $gene(keys(%genes))
		{
			if ($sequences{$gene})
			{
				my $seq = $sequences{$gene};
				my $num_gene = $genes{$gene};
				my $seq_length = length($seq);
				
				my $consensus_seq = $sequences{$gene};
					
				my $nb_polymorphisms = $snp_stats{$num_gene}{"positions"};
				my %positions = %$nb_polymorphisms;
				my @sorted_positions = sort {$a <=> $b } keys(%positions);
					
	
				############################################################
				# get flanking sequences
				############################################################
				my $previous_pos = 0;
				my $previous_position = 0;
				my $consensus_seq2 = "";
					
				foreach my $pos(@sorted_positions)
				{
					my $position;
					if ($pos =~/(\d+)\.*\d*/)
					{
						$position = $1;
					}
					my $size = $position - $previous_position - 1;
					my $subseq = substr($consensus_seq,$previous_position,$size);
						
	
					$previous_position = $position;
					$previous_pos = $pos;
			
					$flanking{$num_gene}{"positions"}{$pos}{"seq_before"} = $subseq;
					$consensus_seq2 .= $flanking{$num_gene}{"positions"}{$pos}{"seq_before"};
					$consensus_seq2 .= $snp_stats{$num_gene}{"positions"}{$pos}{"code_snp"};
				}
				if ($previous_pos)
				{
					my $size = length($consensus_seq) - $previous_position;
					my $subseq = substr($consensus_seq,$previous_position,$size);
					$flanking{$num_gene}{"final_flank"} = $subseq;
						
					$consensus_seq2 .= $flanking{$num_gene}{"final_flank"};
				}
				
					
					
				#####################################################
				# create FASTA alignment
				#####################################################


				if (!-d "$tmpdir/alignment$num")
				{
					mkdir "$tmpdir/alignment$num";
				}
				if (!-d "$tmpdir/depth$num")
				{
					mkdir "$tmpdir/depth$num";
				}
				open(ALIGNMENT,">$tmpdir/alignment$num/$gene.align.fas");
				open(DEPTH_FOR_SEQUENCE,">$tmpdir/depth$num/$gene.depth");
				foreach my $ind(keys(%individus))
				{
					my $new_gene = $gene;
					$new_gene =~s/_/-/g;
					print ALIGNMENT ">$ind" . "_" . $new_gene . "\n";
					print DEPTH_FOR_SEQUENCE ">$ind" . "_" . $new_gene . "\n";
					#######################################################################
					# if depth of coverage, reconstitution of sequence for each individual
					#######################################################################
					if (-e "$depth_file" && keys(%depth) > 0)
					{
						my @nts = split("",$consensus_seq);
						my $position_in_sequence = 0;
						foreach my $nt(@nts)
						{
							$position_in_sequence++;
							my $is_snp_position = 0;
							foreach my $pos(@sorted_positions)
							{
								if ($position_in_sequence == $pos)
								{
									$is_snp_position = 1;
								}
							}
							my $current_depth = $depth{$gene}{$ind}{$position_in_sequence};
							
							if ($ind eq $name_ref)
							{
								print DEPTH_FOR_SEQUENCE "1 ";
							}
							else
							{
								print DEPTH_FOR_SEQUENCE "$current_depth ";
							}
							
							if ($is_snp_position)
							{
								if ($genotyping{$num_gene}{$position_in_sequence}{$ind}{"code_genotype"} eq '0')
								{
									print ALIGNMENT "?";
								}
								else
								{
									if ($ind eq $name_ref)
									{
										print ALIGNMENT $nt;
									}
									else
									{
										print ALIGNMENT $genotyping{$num_gene}{$position_in_sequence}{$ind}{"code_genotype"};
									}
									
									
								}
							}
							else
							{
								if ($current_depth >= $min_coverage or $ind eq $name_ref)
								{
									print ALIGNMENT $nt;
								}
								else
								{
									print ALIGNMENT "?";
								}
							}
						}
						print ALIGNMENT "\n";
						print DEPTH_FOR_SEQUENCE "\n";
					}
					#######################################################################
					# without depth of coverage, use of flanking regions and insertion between each SNP
					#######################################################################
					else
					{
						foreach my $pos(@sorted_positions)
						{
							print ALIGNMENT $flanking{$num_gene}{"positions"}{$pos}{"seq_before"};
								
							if ($genotyping{$num_gene}{$pos}{$ind}{"code_genotype"} eq '0')
							{
								print ALIGNMENT "?";
							}
							else
							{
								print ALIGNMENT $genotyping{$num_gene}{$pos}{$ind}{"code_genotype"};
							}	
						}
						print ALIGNMENT $flanking{$num_gene}{"final_flank"} . "\n";
					}	
				}
				close(ALIGNMENT);
				close(DEPTH_FOR_SEQUENCE);
			}
		}
		system("rm -rf $vcf_file");
		system("rm -rf $depth_file_tmp");
		print "$depth_file_tmp\n";
	}
}
close(LS);

close(VCF_LIGHT);



my $nb_genes = 0;
my $num = 0;



#########################################################
# if reference file has been provided, manage sequences
#########################################################


my $j = 1;
while(-d "$tmpdir/alignment$j")
{
	chdir("$tmpdir/alignment$j");
	system("zip alignments$j.zip ./*.fas >>zip.log 2>&1");
	system("mv alignments$j.zip ../$tmpdir.alignments$j.zip");
	chdir("../..");
	$j++
}

$j = 1;
while(-d "$tmpdir/depth$j")
{
	chdir("$tmpdir/depth$j");
	system("zip depths$j.zip ./*.depth >>zip.log 2>&1");
	system("mv depths$j.zip ../$tmpdir.depths$j.zip");
	chdir("../..");
	$j++;
}
chdir("$tmpdir");
system("zip $tmpdir.results.zip *.zip >>zip.log 2>&1");
system("mv  $tmpdir.results.zip ../");
chdir("..");
if (-d "$tmpdir" && $tmpdir=~/\w+/)
{
	system("rm -rf ./$tmpdir");
}
