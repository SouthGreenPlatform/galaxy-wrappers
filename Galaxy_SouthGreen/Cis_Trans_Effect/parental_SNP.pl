#!/usr/bin/perl

###################
#
#
#
# Intellectual property belongs to IRD and SouthGreen developpement plateform
#
# Written by Simon Foulquier 
#
#################### 

use strict;
use warnings;
 

use feature qw{ say };
use Getopt::Long; 

# Declaration variables
my ($vcf_in,$list,$nom_parent1,$nb_parent1,$nom_parent2,$nb_parent2);

if (scalar @ARGV < 6) 
{
	print "\n\nUsage : parental_SNP.pl -i XXX.vcf -o YYY -s1 CC -s2 FG -n1 4 -n2 5\n\n
	************************** commande obligatoire !!!! ****************************\n\n
	-i vcf file : chemin d acces contenant le fichier VCF obtenue a l'aide des BAM des parents (exemple: /home/foulquie/palmier/BC.VCF) \n\n
	-o output_snp_list : le nom que l'on veut attribuer a la liste de SNP parentaux \n\n
	-s1 le nom du 1er parent\n\n
	-s2 le nom du 2eme parent\n\n
	-n1 le nombre d\'individus du parent 1\n\n
	-n2 le nombre d\'individus du parent 2\n\n
	**************** Option suplementaire (par defaut si absente) !!!! ************** \n\n
	-m1 le nombre maximal de valeurs manquantes pour le parent 1(par default 1)\n\n
	-m2 le nombre maximal de valeurs manquantes pour le parent 2(par default 1)\n\n
	-g le nombre minimal d'occurrences pour qu'un genotype soit pris en consideration (par parent) (par default 1)\n\n
	-p la profondeur de couverture minimale (par individu) soutenant un SNP (par default 10)\n\n
	-q la qualite min du Phred-scale quality score(40 correspond a une probabilite de 1/100000 de se tromper de genotype a un locus donne (40 par default)\n\n
	**************************************************************************************** \n\n";
	exit;
}
# valeur par default
my $max_missing_value1= 1;
my $max_missing_value2= 1;
my $min_same_genotype_value = 1;
my $min_depth_ind_value = 10;
my $qual_min = 40;

# Gestion des options placées lors lancement du script
GetOptions("i|input=s"=>\$vcf_in,
			"o|output=s"=>\$list,
			"s1|nameP1=s"=>\$nom_parent1,
			"s2|nameP2=s"=>\$nom_parent2,
			"n1|nbP1=s"=>\$nb_parent1,
			"n2|nbP2=s"=>\$nb_parent2,
			"m1|max_miss_val_P1=s"=>\$max_missing_value1,#default value : 1
			"m2|max_miss_val_P2=s"=>\$max_missing_value2,#default value : 1
			"g|occurence_genotype=s"=>\$min_same_genotype_value,#default value : 1
			"p|profondeur_min=s"=>\$min_depth_ind_value,#default value : 10
			"q|phredlike_quality_score_min=s"=>\$qual_min, #default value : 40
);

#mise en place des valeurs par default si les option ne sont pas tapé par l'utilisateur
if ($max_missing_value1) {
    say $max_missing_value1;
}
if ($max_missing_value2) {
    say $max_missing_value2;
}
if ($min_same_genotype_value) {
    say $min_same_genotype_value;
}
if ($min_depth_ind_value) {
    say $min_depth_ind_value;
}
if ($qual_min) {
    say $qual_min;
}

print "le nombre maximal de valeurs manquantes pour le parent ".$nom_parent1." : ".$max_missing_value1."\n";
print "le nombre maximal de valeurs manquantes pour le parent ".$nom_parent2." : ".$max_missing_value2."\n";
print "le nombre minimal d'occurrences pour qu'un genotype soit pris en consideration   : ".$min_same_genotype_value."\n";
print "la profondeur de couverture minimale (par individu) soutenant un SNP   : ".$min_depth_ind_value."\n";
print "la qualite min du Phred-scale quality score   : ".$qual_min."\n";
open(HOMEOLIST,">$list") or die("Fail to open $list : $!");

open(TRACE,">trace") or die("Fail to open trace");
open(VCFIN,"<$vcf_in") or die("Fail to open $vcf_in : $!");

my%traceancestor1;
my%traceancestor2;
my $indivu_total = $nb_parent1 + $nb_parent2;

print HOMEOLIST "CDS\tPOS\t$nom_parent1 allele1\t$nom_parent1 allele2\t$nom_parent1 allele3\t$nom_parent2 allele1\t$nom_parent2 allele2\t$nom_parent2 allele3\n";
print TRACE "ancestor\tchromosome\tposition\tallele\n";

my $cmpt_parent1 = 8 + $nb_parent1;  ##### compteur pour les colonnes correspondant au parent 1 ###### 
my $cmpt_parent2 = 1 + $cmpt_parent1;  ##### compteur pour les colonnes correspondant au parent 2 ###### 


while(my $vcf_line = <VCFIN>)
{
	next if ($vcf_line =~ m/^#/);
	chomp($vcf_line);
	my @infos_line = split(m/\t/,$vcf_line);
	my $col_fich = scalar (@infos_line);
	if (scalar(@infos_line) != (9 + $indivu_total))
	{
		print scalar(@infos_line). " " . (9 + $indivu_total);
		die "\n \n \n La somme des individus ne correspond pas au nombre individus total dans le vcf !!!!!! \n \n \n";
	}
	my @parent1_genotypes = @infos_line[9..$cmpt_parent1];
	my @parent2_genotypes = @infos_line[$cmpt_parent2..$#infos_line];
	my @parent_genotypes = (\@parent1_genotypes, \@parent2_genotypes);
		
	my $individus_pris_en_cmpt = 0;
	my $cds = $infos_line[0];
        my $pos = $infos_line[1];
        my $qual = $infos_line[5];
	
	next if ($qual < $qual_min);
	
	#print "ICI1 : $vcf_line\n$cds $pos et qualite : $qual\n";
	
	my $outer=0;

	my %genotypes_parent1_array;
	my %genotypes_parent2_array;
	my %genotypes_parent1_less_depth_array;
	my %genotypes_parent2_less_depth_array;
	my @parent1_attributed_nucleotide;
	my @parent2_attributed_nucleotide;
	my $indiv_cmpt = 0;
	my $nb_indiv1 = 0 ;
	my $nb_indiv2 = 0 ;
	
	OUTER:
	foreach my $parent_genotypes_ref (@parent_genotypes)
	{
		my @parent_genotypes = @$parent_genotypes_ref;
		my $missing_value_count_p1=0;
		my $missing_value_count_p2=0;
		my $verif_cpt = 0;
		my $genotype_homogene;
		my %table_en_cours;
		my $nb =0;
		
		foreach my $parent_genotype (@parent_genotypes)
		{	
			my @parent_genotype_split = split(m/:/,$parent_genotype);
			my $real_genotype = $parent_genotype_split[0];
			my $individu_depth = $parent_genotype_split[2];
			
			#print "ICI2 : genotype -> $real_genotype\n";
			#print "ICI2 : profondeur_par_individu -> $individu_depth\n";
				
			########### si la valeure est absente ##############
			if ($real_genotype =~ m/\.\/\./)
			{
				if ($indiv_cmpt < $nb_parent1) 
				{
					$missing_value_count_p1++;
					$indiv_cmpt++;
					#print "ICI4 : missing value\n";
					if ($missing_value_count_p1 > $max_missing_value1)  ##### par default $max_missing_value1 = 1;#####
					{
						$outer = 1;
						last OUTER;
					}
				}
				else
				{
					$missing_value_count_p2++;
					#print "ICI4 : missing value\n";
					if ($missing_value_count_p2 > $max_missing_value2)  ##### par default $max_missing_value2 = 1;#####
					{
						$outer = 1;
						last OUTER;
					}
				}
			}
			elsif ($indiv_cmpt < $nb_parent1) ###### si 1er parent ######
			{
				$genotypes_parent1_array{$real_genotype}{$individu_depth}++;
				$table_en_cours{$real_genotype}++;
				$indiv_cmpt++;
				$nb_indiv1 ++; ##### nb individu parent 1 retenu ######
			}
			else ###### sinon c'est le 2eme parent  ########
			{	
				foreach my $k (keys(%genotypes_parent1_array)) 
				{
					if (verif_allele_spe($k,$real_genotype) == 0 )
					{
						$outer = 1;
						last OUTER;					
					}
				}
				$genotypes_parent2_array{$real_genotype}{$individu_depth}++; ###### si pas sortie de ligne -----> 2eme parent ######
				$table_en_cours{$real_genotype}++;
				$nb_indiv2++;
			}
		}
		foreach my $val (values(%table_en_cours))  ####### verifie que l'on possede au le min d'un meme genotype ########
		{
			if ($val >= $min_same_genotype_value)
			{
				$verif_cpt = 1;
				last;
			}
		}
		if ($verif_cpt == 0)
		{
			$outer = 1;
			last OUTER;
		}
		foreach my $real_geno(keys %genotypes_parent1_array)
		{
			foreach my $depth (keys %{$genotypes_parent1_array{$real_geno}})
			{
				if ($depth < $min_depth_ind_value)######### si la profondeur n'est pas suffisante ###########
				{
					my $count = scalar ($genotypes_parent1_array{$real_geno}{$depth});
					$missing_value_count_p1 += $count;
					delete( $genotypes_parent1_array{$real_geno}{$depth} );
					$indiv_cmpt++;
					#print "ICI4 : missing value cause depth in parent 1 \n";
					if ($missing_value_count_p1 > $max_missing_value1)  ##### par default $max_missing_value1 = 1;#####
					{
						$outer = 1;
						last OUTER;
					}
				}
			}
		}
		foreach my $real_geno2(keys %genotypes_parent2_array)
		{
			foreach my $depth2 (keys %{$genotypes_parent2_array{$real_geno2}})
			{
				if ($depth2 < $min_depth_ind_value)######### si la profondeur n'est pas suffisante ###########
				{
					my $count2 = scalar ($genotypes_parent2_array{$real_geno2}{$depth2});
					$missing_value_count_p2 += $count2;
					delete( $genotypes_parent2_array{$real_geno2}{$depth2} );
					if ($missing_value_count_p2 > $max_missing_value2)  ##### par default $max_missing_value2 = 1;#####
					{
						$outer = 1;
						last OUTER;
					}
				}
			}
		}	
	}
	next if ($outer == 1);
	my $allele_ref = $infos_line[3];
	my $allele_alt = $infos_line[4];
	
				##### nouvelle hash sans la profondeur #######
	%genotypes_parent1_less_depth_array = delete_depth_hash(%genotypes_parent1_array);
	%genotypes_parent2_less_depth_array = delete_depth_hash(%genotypes_parent2_array);	
	
	my @array_allele_freq_p1 = hash_allele_freq_parent(%genotypes_parent1_less_depth_array);
	my @array_allele_freq_indivi_p1 = hash_allele_freq_parent_v2(%genotypes_parent1_less_depth_array);
	my $allele_vide = "0";
	for (my $r = 1; $r <= 3; $r++) ###### boucle permettant de mettre 0 en frequence au allele "X" ######
	{
		if (scalar(@array_allele_freq_indivi_p1) == $r)
		{
			for ( my $d = $r; $d <= 2; $d++) 
			{
				push(@array_allele_freq_indivi_p1, $allele_vide); 
			}
		}			
    }

	my @array_allele_freq_p2 = hash_allele_freq_parent(%genotypes_parent2_less_depth_array);
	my @array_allele_freq_indivi_p2 = hash_allele_freq_parent_v2(%genotypes_parent2_less_depth_array);
	my $allele_vide2 = "0";
	for (my $y = 1; $y <= 3; $y++) 
	{
		if (scalar(@array_allele_freq_indivi_p2) == $y)
		{
			for ( my $d = $y; $d <= 2; $d++) 
			{
				push(@array_allele_freq_indivi_p2, $allele_vide2);
			}
		}			
    }
	
	foreach my $ale (@array_allele_freq_p1)
	{
		my $ale_correspond;
		$ale_correspond = corresponding_allele($ale,$allele_ref,$allele_alt);
		push(@parent1_attributed_nucleotide, $ale_correspond);
	}
	
	foreach my $ale2 (@array_allele_freq_p2)
	{
		my $ale_correspond;
		$ale_correspond = corresponding_allele($ale2,$allele_ref,$allele_alt);
		push(@parent2_attributed_nucleotide, $ale_correspond);
	}
	@parent1_attributed_nucleotide = transform_list_allele(@parent1_attributed_nucleotide);
	@parent2_attributed_nucleotide = transform_list_allele(@parent2_attributed_nucleotide);

	print HOMEOLIST $cds . "\t" . $pos . "\t"  
	.$parent1_attributed_nucleotide[0] . "\t" . $parent1_attributed_nucleotide[1] . "\t" . $parent1_attributed_nucleotide[2] . "\t" 
	.$parent2_attributed_nucleotide[0] . "\t" . $parent2_attributed_nucleotide[1] . "\t" . $parent2_attributed_nucleotide[2] . "\n";

	$traceancestor1{$nom_parent1}{$cds}{$pos} = $parent1_attributed_nucleotide[0];
	$traceancestor2{$nom_parent2}{$cds}{$pos} = $parent2_attributed_nucleotide[0];
}

foreach my$parent(sort keys %traceancestor1){
	foreach my$cds(sort keys %{$traceancestor1{$parent}}){
		foreach my$pos(sort keys %{$traceancestor1{$parent}{$cds}}){
			print TRACE "$parent\t$cds\t$pos\t$traceancestor1{$parent}{$cds}{$pos}\n";
		}
	}
}
foreach my$parent(sort keys %traceancestor2){
	foreach my$cds(sort keys %{$traceancestor2{$parent}}){
		foreach my$pos(sort keys %{$traceancestor2{$parent}{$cds}}){
			print TRACE "$parent\t$cds\t$pos\t$traceancestor2{$parent}{$cds}{$pos}\n";
		}
	}
}

close(VCFIN);
close(HOMEOLIST);

##################################################  FUNCTIONS  ############################################################################
sub verif_allele_spe  #### permet de verifier entre 2 genotypes si il y au moins une allele commune 
{
	my ($geno_espe1,$geno_espe2) = @_;
	my $verif;
	my @allele_esp1 = split('/',$geno_espe1);
	my @allele_esp2 = split('/',$geno_espe2);
	for my $val (@allele_esp2)
	{
		if( ($allele_esp1[0] == $val) or  ($allele_esp1[1] == $val)) 
		{
			$verif = 0;
			last;
		} 
		else 
		{
			$verif = 1;
		}
	}	
	return $verif;
}
sub delete_depth_hash
{
	my (%hash_genotype_depth) = @_;
	my %new_hash;
	
	foreach my $geno(keys %hash_genotype_depth)
	{
		foreach my $dp(keys %{$hash_genotype_depth{$geno}})
		{
			my $freq = scalar($hash_genotype_depth{$geno}{$dp});
			$new_hash{$geno} +=$freq;
		}
	}
	return %new_hash;
}
sub hash_allele_freq_parent ##### sort les allele par ordre decroissant de frequence #######
{
	my (%hash_genotype) = @_;
	my $corresp_allele;
	my @tri_allele;
	my %hash_allele_freq;

	foreach my $k ( sort ({ $hash_genotype{$a} <=> $hash_genotype{$b} } keys %hash_genotype)) #### trie  les genotypes par frequence ######
	{
		my @single_allele = split('/',$k);
		if ($single_allele[0] == $single_allele[1])
		{
			my $freq = $hash_genotype{$k} * 2;
			$hash_allele_freq{$single_allele[0]} += $freq;
		}
		else
		{
			my $freq = $hash_genotype{$k};
			$hash_allele_freq{$single_allele[0]} += $freq;
			$hash_allele_freq{$single_allele[1]} += $freq;
		}
	}
	foreach my $ke ( reverse sort ({ $hash_allele_freq{$a} <=> $hash_allele_freq{$b} } keys %hash_allele_freq)) #### trie  les alleles par frequence ######
	{
		push(@tri_allele, $ke);
	}
	return @tri_allele;
}
sub hash_allele_freq_parent_v2 ##### sort la frequence allele par ordre decroissant #######
{
	my (%hash_genotypev2) = @_;
	my $corresp_allelev2;
	my @freq_allelev2;
	my %hash_allele_freqv2;

	foreach my $kv2 ( sort ({ $hash_genotypev2{$a} <=> $hash_genotypev2{$b} } keys %hash_genotypev2)) #### trie  les genotypes par frequence ######
	{
		my @single_allelev2 = split('/',$kv2);
		if ($single_allelev2[0] == $single_allelev2[1])
		{
			my $freqv2 = $hash_genotypev2{$kv2} * 2;
			$hash_allele_freqv2{$single_allelev2[0]} += $freqv2;
		}
		else
		{
			my $freqv2 = $hash_genotypev2{$kv2};
			$hash_allele_freqv2{$single_allelev2[0]} += $freqv2;
			$hash_allele_freqv2{$single_allelev2[1]} += $freqv2;
		}
	}
	foreach my $kev2 ( reverse sort ({ $hash_allele_freqv2{$a} <=> $hash_allele_freqv2{$b} } keys %hash_allele_freqv2)) #### trie  les alleles par frequence ######
	{
		push(@freq_allelev2, $hash_allele_freqv2{$kev2});
	}
	return @freq_allelev2;
}

sub corresponding_allele ##### retrouve l'allele corespondant au genotipe GATK #####
{
	my ($signle_al_num,$ref_allele,$alt_allele) = @_;
	my $corresp_allele;
	if ($signle_al_num =~ m/0/)
	{
		$corresp_allele = $ref_allele;
	}
	elsif ($signle_al_num =~ m/1/)
	{
		$corresp_allele = substr($alt_allele,0,1);
	}
	elsif ($signle_al_num =~ m/2/)
	{
		$corresp_allele = substr($alt_allele,2,1);
	}
	elsif ($signle_al_num =~ m/3/)
	{
		$corresp_allele = substr($alt_allele,4,1);
	}
	return $corresp_allele;
}

sub transform_list_allele  ##### met des X si pas de 2eme ou 3eme allele ####
{
	my (@liste_al) = @_;
	my $allele_vide = "X";
	for (my $i = 1; $i <= 3; $i++) 
	{
		if (scalar(@liste_al) == $i)
		{
			for ( my $d = $i; $d <= 2; $d++) 
			{
				push(@liste_al, $allele_vide);
			}
		}			
    }
	return @liste_al;
}
exit;
