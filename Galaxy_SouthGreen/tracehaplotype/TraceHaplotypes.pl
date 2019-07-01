#!/usr/bin/perl

use strict;
use Getopt::Long;

my $usage = qq~Usage:$0 <args> [<opts>]

where <args> are:

    -v, --vcf          <VCF input>
    -h, --hybrid_form  <Comma-separated form: hybrid,parent1,parent2>
    -o, --out          <output basename> 	
	
      <opts> are:

	-w, --window       <window size. Default:500000 bp)
    -s, --summarize    <Number of windows to be summarized. Default: 4)
    -f, --filter       <Use only variants located in genes (no/yes)(if any annotation). Default: no
~;
$usage .= "\n";

my ($vcf,$hybrid_formula,$out,$window,$summarize,$filter);
GetOptions(
"vcf=s"        => \$vcf,
"out=s"        => \$out,
"window=s"     => \$window,
"summarize=s"  => \$summarize,
"filter=s"     => \$filter,
"hybrid_form=s"=> \$hybrid_formula
);


die $usage
  if ( !$vcf || !$hybrid_formula || !$out);

my $window_size = 500000;

my $summarized = 4;
if ($window){
	$window_size = $window;
}
if ($summarize){
	$summarized = $summarize;
}
my $filtering = "no";
if ($filter eq "yes"){
	$filtering = "yes";
}

my ($hybride,$parent1,$parent2) = split(",",$hybrid_formula);




my %colors = (
	1 => "#013CAB",
	2 => "#DACA08",
);

my %max;
my %positions;



my $n1 = 0;
my $n2 = 0;
my %order_genes;
my %hash;
my $nb = 0;
open(F,$vcf);
while(<F>){
	my $line = $_;
	$line =~s/\n//g;$line =~s/\r//g;
	my @i = split(/\t/,$line);
	if (/#CHROM/){
		for (my $k = 9; $k <= $#i; $k++){
			my $val = $i[$k];
			if ($val eq "$hybride"){
				$hybride = $k;
			}
			if ($val eq "$parent1"){
				$parent1 = $k;
			}
			if ($val eq "$parent2"){
                                $parent2 = $k;
                        }
		}
	}
	elsif (/^#/){

	}
	else{
		my $chr = $i[0];	
		my $pos = $i[1];
		my $ref = $i[3];
		my $alt = $i[4];
		
		my $genename = $i[7];
		if ($genename =~/(Ciclev10\d+m.g)/){
			$genename = $1;
		}
		else{
			if ($filtering eq "yes"){
				next;
			}
		}
		
		
		$chr =~s/scaffold_//g;
		
		my $gene = int($pos / $window_size);
		

		my $hybride_geno = $i[$hybride];
		my $parent1_geno = $i[$parent1];
		my $parent2_geno = $i[$parent2];
		
		#if ($hybride_geno =~/1/ or $parent1_geno =~/1/  or $parent2_geno =~/1/  ){
		#	$positions_found{$pos}=1;
		#}
		if (!$hash{$chr}{$gene}){$order_genes{$chr}.=",$gene";}
		#if ($gene !~/Ciclev10011363m.g/){next;}
	
		# incompatibility for parent1
		if ($parent1_geno !~/[01]/){
		}
		# incompatibility for parent2
		elsif ($parent2_geno !~/[01]/){
		
		}
		if ($hybride_geno eq "1\|1"){
			my $nb_compat1 = 0;
			my $nb_compat2 = 0;
			if ($parent1_geno =~/1\/1/ && $parent2_geno =~/1\/1/){
				$nb_compat2++;
				$nb_compat1++;
			}
			if ($parent1_geno =~/\.\/\./ && $parent2_geno =~/1\/1/){
                                $nb_compat1++;
                        }
			if ($parent2_geno =~/\.\/\./ && $parent1_geno =~/1\/1/){
                                $nb_compat2++;
                        }
			if ($parent1_geno =~/0\/1/ && $parent2_geno =~/1\/1/){
                                $nb_compat1++;
                        }
			if ($parent1_geno =~/1\/1/ && $parent2_geno =~/0\/1/){
                                $nb_compat2++;
                        }

		#	if ($nb_compat2 or $nb_compat1){print "$pos  $nb_compat1 $nb_compat2\n";}
			$hash{$chr}{$gene}{"nb_incompatiblity_modele1"}+=$nb_compat1;
			$hash{$chr}{$gene}{"nb_incompatiblity_modele2"}+=$nb_compat2;
		}
		elsif ($hybride_geno eq ".\|."){
			my $nb_compat1 = 0;
			my $nb_compat2 = 0;
			#if ($parent1_geno =~/\.\/\./ && $parent2_geno =~/1\/1/){
			#	$nb_compat1++;
			#}
			#if ($parent1_geno =~/0\/0/ && $parent2_geno =~/1\/1/){
                        #        $nb_compat1++;
                        #}
			if ($parent2_geno =~/1\/1/){
				$nb_compat1++;
			}
			if ($parent1_geno =~/1\/1/){
				$nb_compat2++;
			}
		#	if ($nb_compat2 or $nb_compat1){print "$pos  $nb_compat1 $nb_compat2\n";}
			$hash{$chr}{$gene}{"nb_incompatiblity_modele1"}+=$nb_compat1;
			$hash{$chr}{$gene}{"nb_incompatiblity_modele2"}+=$nb_compat2;
		}
		elsif ($hybride_geno eq "0\|0"){
			my $nb_compat1 = 0;
			my $nb_compat2 = 0;
			if ($parent1_geno =~/1\/1/ && $parent2_geno =~/0\/0/){
				$nb_compat2++;
				$nb_compat1++;
			}
			if ($parent1_geno =~/0\/0/ && $parent2_geno =~/1\/1/){
				$nb_compat2++;
				$nb_compat1++;
			}
			if ($parent1_geno =~/1\/1/ && $parent2_geno =~/0\/1/){
                                $nb_compat2++;
                                $nb_compat1++;
                        }
			if ($parent2_geno =~/1\/1/ && $parent1_geno =~/0\/1/){
                                $nb_compat2++;
                                $nb_compat1++;
                        }
			if ($parent1_geno =~/1\/1/ && $parent2_geno =~/\.\/\./){
                                $nb_compat2++;
                                $nb_compat1++;
                        }
			if ($parent2_geno =~/1\/1/ && $parent1_geno =~/\.\/\./){
                                $nb_compat2++;
                                $nb_compat1++;
                        }
		#	if ($nb_compat2 or $nb_compat1){print "$pos  $nb_compat1 $nb_compat2\n";}
			$hash{$chr}{$gene}{"nb_incompatiblity_modele1"}+=$nb_compat1;
			$hash{$chr}{$gene}{"nb_incompatiblity_modele2"}+=$nb_compat2;
		}
		elsif (($hybride_geno eq "0\|1" or $hybride_geno eq "1\|0") ){
			my $nb_compat1 = 0;
			my $nb_compat2 = 0;
			if ($parent1_geno =~/1\/1/ && $parent2_geno =~/0\/0/){
				$nb_compat2+=2;
			}
			if ($parent1_geno =~/1\/1/ && $parent2_geno =~/0\/1/){
                                $nb_compat2++;
                        }
			if ($parent1_geno =~/0\/1/ && $parent2_geno =~/0\/0/){
				$nb_compat2++;
			}
			
			# if ($parent2_geno eq "0/0" && $parent1_geno =~/1/){
				# print "$hybride_geno $parent1_geno $parent2_geno\n";
			# }
			
		
			

			if ($parent2_geno =~/1\/1/ && $parent1_geno =~/0\/0/){
                                $nb_compat1+=2;
                        }
                        if ($parent2_geno =~/1\/1/ && $parent1_geno =~/0\/1/){
                                $nb_compat1++;
                        }
			if ($parent2_geno =~/0\/1/ && $parent1_geno =~/0\/0/){
				$nb_compat1++;
			}

			if ($parent2_geno =~/0\/0/ && $parent1_geno =~/0\/0/){
				$nb_compat1++;
				$nb_compat2++;
			}
			if ($parent2_geno =~/1\/1/ && $parent1_geno =~/1\/1/){
                                $nb_compat1++;
                                $nb_compat2++;
                        }
		
			if ($parent1_geno =~/0\/1/ && $parent2_geno =~/0\/1/){
				# do nothing
			}
		#	if ($nb_compat2 or $nb_compat1){print "$pos $nb_compat1 $nb_compat2\n";}
			$hash{$chr}{$gene}{"nb_incompatiblity_modele1"}+=$nb_compat1;
			$hash{$chr}{$gene}{"nb_incompatiblity_modele2"}+=$nb_compat2;
		}		
		
	}
}
close(F);



#############################################################################
# incompatibility window by window
#############################################################################
my %genes_OK;
my %suites;
open(G,">$out.haplotypes.genes.txt");
foreach my $chr(sort {$a<=>$b}keys(%order_genes)){
	my @genelist = split(",",$order_genes{$chr});
	my $modele_region = 0;
	my %hashmodele;
	my $previous_modele = 0;
	my $num = 0;
	foreach my $gene(@genelist){
		if ($gene !~/\w+/){next;}

		my $modele = 0;
		if ($hash{$chr}{$gene}{"nb_incompatiblity_modele1"} > $hash{$chr}{$gene}{"nb_incompatiblity_modele2"}){
			$modele = 2;
		}
		elsif ($hash{$chr}{$gene}{"nb_incompatiblity_modele2"} > $hash{$chr}{$gene}{"nb_incompatiblity_modele1"}){
			$modele = 1;
		}
		else{
			$modele = 0;
		}

		
		if ($modele){
		
			$hashmodele{$modele}++;
			if ($hashmodele{$modele} >= 2){
				$modele_region = $modele;
			}
			
			if ($modele ne "0" && $modele != $previous_modele){
				$hashmodele{$modele}=0;
			}
		}

		my $start = $gene*$window_size;
		my $end = $start + $window_size;
		
		print G "$chr $gene $start $end $modele\n";
		if ($modele > 0){
			$suites{$chr}.= $modele;
			$genes_OK{$chr}{$num}="$start-$end";
			$num++;
		}
		$previous_modele = $modele;
	}
}
close(G);

open(H,">$out.haplotypes.genes.corrected.txt");
foreach my $chr(sort {$a<=>$b}keys(%order_genes)){
	my $suite = $suites{$chr};
	my @values = split("",$suite);
	my $previous_modele = 0;
	my $current;

	for (my $i = 0; $i <= $#values; $i++){
		my $value = $values[$i];
		my $nextvalue = $values[$i+1];
		my $modele;
		
		if ($value != $current){
			if ($nextvalue eq $value){
				$modele = $value;
			}
			else{
				$modele = $current;
			}
		}
		else{
			$modele = $value;
		}
		$current = $modele;
		my ($start,$end) = split("-",$genes_OK{$chr}{$i});
		print H "$chr $start $end $current\n";
	}
}
close(H);


#############################################################################
# incompatibility by sliding window (sum for N genes)
#############################################################################
open(G2,">$out.haplotypes.window.txt");
foreach my $chr(sort {$a<=>$b}keys(%order_genes)){
	my @genelist = split(",",$order_genes{$chr});
	my $modele_region = 0;
	my %hashmodele;
	my $previous_modele = 0;
	my $num = 0;
	for (my $i = 0; $i <= $#genelist; $i++){
		my $gene1 = $genelist[$i];
		if ($gene1 !~/\w+/){next;}
		my $sum1 = 0;
		my $sum2 = 0;
		for (my $k=0; $k < $summarized; $k++){
			my $gene = $genelist[$i+$k];
			$sum1+=$hash{$chr}{$gene}{"nb_incompatiblity_modele1"};
			$sum2+=$hash{$chr}{$gene}{"nb_incompatiblity_modele2"};
		}
		my $modele = 0;
		if ($sum1 > $sum2){
			$modele = 2;
		}
		elsif ($sum1 < $sum2){
			$modele = 1;
		}
		else{
			$modele = 0;
		}
		my $start = $gene1*$window_size;
		my $end = $start + $window_size;
		if ($modele > 0){
			$suites{$chr}.= $modele;
			$genes_OK{$chr}{$num}="$start-$end";
			$num++;
		}
		print G2 "$chr $start $end $modele\n";
	}
}
close(G2);


open(K,">$out.haplotypes.blocks.txt");
my $previous_modele = 0;
my $previous_chr;
my $start_block = 1;
my $end_block;
my $chrom;
my %blocks;
my @infolines;
open(G,"$out.haplotypes.genes.corrected.txt");
#open(G,"$out.haplotypes.window.txt");
while(<G>){
	my $line = $_; $line =~s/\n//g;$line =~s/\r//g;
	my ($chr,$start,$end,$modele) = split(" ",$line);
	
	if ($previous_chr && $chr ne $previous_chr && $end_block > 0){
	
		if ($end_block > $max{$chrom}){
			$max{$chrom} = $end_block;
		}
		
		print K "$previous_chr 0 $start_block $end_block $colors{$previous_modele}\n";
		$blocks{$previous_chr}{"$start_block-$end_block"} = $previous_modele;
		push(@infolines,"$previous_chr 0 $start_block $end_block $colors{$previous_modele}");
		$start_block = 1;
	}
	
	$chrom = $chr;
	
	if ($modele ne $previous_modele && $previous_modele != 0 && $start > 0){
			my $end_block = $start;
			
			if ($end_block > $max{$chrom}){
				$max{$chrom} = $end_block;
			}
		
			print K "$chr 0 $start_block $end_block $colors{$previous_modele}\n";
			$blocks{$chr}{"$start_block-$end_block"} = $previous_modele;
			push(@infolines,"$chr 0 $start_block $end_block $colors{$previous_modele}");
			$start_block = $end_block;
	}
	$end_block = $end;
	$previous_modele = $modele;
	$previous_chr = $chr;
}
if ($end_block > $max{$chrom}){
	$max{$chrom} = $end_block;
}
$blocks{$chrom}{"$start_block-$end_block"} = $previous_modele;
print K "$chrom 0 $start_block $end_block $colors{$previous_modele}\n";
push(@infolines,"$chrom 0 $start_block $end_block $colors{$previous_modele}");

foreach my $l(@infolines){
	#if ($l =~/^(\w+) (\d) (\d+) (\d+) (.*)\s$/){
	if ($l =~/^(\w+) (\d) (\d+) (\d+) #(\w+)/){

		my $chr = $1;
		my $haplo = $2;
		my $start = $3;
		my $end = $4;
		my $modele = "#".$5;
		
		if ($modele eq $colors{"1"}){print K "$chr 1 $start $end ".$colors{"2"}."\n";}
		if ($modele eq $colors{"2"}){print K "$chr 1 $start $end ".$colors{"1"}."\n";}
	}
}
close(G);
close(K);

open(C,">$out.chrom_length.txt");
foreach my $chr(sort {$a<=>$b} keys(%order_genes)){
	print C "$chr $max{$chr} AA\n";
}
close(C);

#############################################################################
# reintroduce haplotype information in the VCF
#############################################################################
open(FO,">$vcf.rephased.vcf");
open(F,$vcf);
while(<F>){
	my $line = $_;
	$line =~s/\n//g;$line =~s/\r//g;
	my @i = split(/\t/,$line);
	if (/#CHROM/){
		for (my $k = 9; $k <= $#i; $k++){
			my $val = $i[$k];
			if ($val eq "$hybride"){
				$hybride = $k;
			}
			if ($val eq "$parent1"){
				$parent1 = $k;
			}
			if ($val eq "$parent2"){
				$parent2 = $k;
			}
		}
		print FO "$line\n";
	}
	elsif (/^#/){
		print FO "$line\n";
	}
	else{
		my $chr = $i[0];	
		my $pos = $i[1];
		my $ref = $i[3];
		my $alt = $i[4];
		my $hybride_geno = $i[$hybride];
		
		if ($hybride_geno =~/0\|1/){
			print FO $chr;
			$chr =~s/scaffold_//g;
			my $blocks_of_chrom = $blocks{$chr};
			my %hashblocks = %$blocks_of_chrom;
			my $model;
			my $inter;
			foreach my $interval(keys(%hashblocks)){
			
				my ($startb,$endb) = split("-",$interval);
				if ($pos <= $endb && $pos >= $startb){
					$model = $blocks{$chr}{"$interval"};
					$inter = $interval;
				}
				
			}
			for (my $k = 1; $k <= $#i; $k++){
				my $val = $i[$k];
				if ($k == $hybride){
					if ($model == 1){
						$val = "0|1";
					}
					elsif ($model == 2){
						$val = "1|0";
					}
					$val =~s/\//\|/g;	
				}
				print FO "\t$val";
			}
			print FO "\n";
		}
		else{
			print FO "$line\n";
		}
	}
}
close(F);
close(FO);



