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
 
use Data::Dumper;
#use Text::Table;
use feature qw{ say };
use Getopt::Long; 

# Declaration variables
my ($bam_file,$homeoSNP_file,$output);
if (scalar @ARGV < 3)

{
	print "\n\nUsage : perl final_geno_cds_nb_reads.pl -i -s -o\n
	************************** commande obligatoire !!!! ****************************\n\n
	-i chemin d acces contenant le fichier BAM (exemple: -i /home/foulquie/palmier/BC.bam) \n\n
	-s chemin d acces a la liste de SNP parentaux (exemple: -s /home/foulquie/palmier/liste/list_SNP) \n\n
	-o nom donné a la matrice de sortie (exemple: -o meso411254) \n\n
	**************** Option suplementaire (par default si absente) !!!! ************* \n\n
	-ns le nombre de SNP minimum par reads necessaire pour etre valide (par default 1 SNP)\n\n
	-mp le score de mapping minimum par read (par default 1 ) \n\n
		***************************************************************************************************\n\n";
	exit;
}

# valeur par default
my $seuil_mapq= 1;
my $min_snp_read = 1;

# Gestion des options placées lors lancement du script
GetOptions("i|input=s"=>\$bam_file,
			"s|liste_SNP=s"=>\$homeoSNP_file,
			"o|output=s"=>\$output,
			"ns|seuil minimum de snp par read =s"=>\$min_snp_read, #default value : 1
			"mp|le score de mapping minimum par read =s"=>\$seuil_mapq);#default value : 1

			
#mise en place des valeurs par default si les option ne sont pas tapé par l'utilisateur

if ($seuil_mapq) {
    say $seuil_mapq;
}
if ($min_snp_read) {
    say $min_snp_read;
}

print "le nombre minimun de SNP par reads : ".$min_snp_read."\n";
print "le score de mapping minimum par read  : ".$seuil_mapq."\n";


open(SAM," samtools view $bam_file |") or die ("Cannot open $bam_file \n");
open(HOMEOSNP,"<$homeoSNP_file") or die ("Cannot open $homeoSNP_file\n");
open(OUTPUT,">$output") or die ("Cannot open $output\n");


#Fichier homeoSNP chargé en mémoire
my $parent1;
my $parent2;
my %homeo_SNPs;
my @attributs;
my %conflict_reads;
while (defined(my $homeolist_line = <HOMEOSNP>))
{
	if ($homeolist_line =~ m/^CDS/)
	{
		@attributs= split(m/\t/,$homeolist_line);
		my $x = scalar (@attributs);
		$parent1 = $attributs[2];
		my $lgp1 = length $parent1;
		$parent1 = substr($parent1,0,($lgp1-8));
		$parent2 = $attributs[5];
		my $lgp2 = length $parent2;
		$parent2 = substr($parent2,0,($lgp2-8));
	}
	chomp($homeolist_line);
	my @homeoSNPs_infos = split(m/\t/,$homeolist_line);
	my $cds = $homeoSNPs_infos[0];
	my $position = $homeoSNPs_infos[1];
	my $parent1_attributed_nucleotide1 = $homeoSNPs_infos[2];
	my $parent1_attributed_nucleotide2 = $homeoSNPs_infos[3];
	my $parent1_attributed_nucleotide3 = $homeoSNPs_infos[4];
	my $parent2_attributed_nucleotide1 = $homeoSNPs_infos[5];
	my $parent2_attributed_nucleotide2 = $homeoSNPs_infos[6];
	my $parent2_attributed_nucleotide3 = $homeoSNPs_infos[7];
	
	if ($parent1_attributed_nucleotide3 ne "X")
	{
		$homeo_SNPs{$cds}{$position}{"parent1_attribution"} = [$parent1_attributed_nucleotide1 , $parent1_attributed_nucleotide2, $parent1_attributed_nucleotide3];
	}
	if ($parent1_attributed_nucleotide2 ne "X")
	{
		$homeo_SNPs{$cds}{$position}{"parent1_attribution"} = [$parent1_attributed_nucleotide1 , $parent1_attributed_nucleotide2];
	}
	else
	{
		$homeo_SNPs{$cds}{$position}{"parent1_attribution"} = [$parent1_attributed_nucleotide1];
	}


	if ($parent2_attributed_nucleotide3 ne "X")
	{
		$homeo_SNPs{$cds}{$position}{"parent2_attribution"} = [$parent2_attributed_nucleotide1 , $parent2_attributed_nucleotide2 , $parent2_attributed_nucleotide3];
	}
	if ($parent2_attributed_nucleotide2 ne "X")
	{
		$homeo_SNPs{$cds}{$position}{"parent2_attribution"} = [$parent2_attributed_nucleotide1 , $parent2_attributed_nucleotide2];
	}
	else
	{
		$homeo_SNPs{$cds}{$position}{"parent2_attribution"} = [$parent2_attributed_nucleotide1];
	}
}

close(HOMEOSNP);



#Parcours du fichier bam et classement des reads

my %count_reads;
my %deja_vu;
my $nombre_de_read_categoriser_total=0;
my $nombre_de_read_avec_conflit_total=0;
my %stats;
my $nb_SNP;


while(defined ( my $line = <SAM>))
{
	next if ($line =~ m/^@/); #Header
	chomp($line);
	my @sam_infos = split(/\t/,$line);
	my $mapq = $sam_infos[4];
	next if ($mapq < $seuil_mapq); #discard read that is mapped to several spots
	my $cds = $sam_infos[2];


	if (! $deja_vu{$cds}) #Initialiser le count read pour le cds si ce n'est pas deja fait
	{
		$count_reads{$cds}{"$parent1"} = 0;
		$count_reads{$cds}{"$parent2"} = 0;
		$count_reads{$cds}{"ambigu"} = 0;
		$count_reads{$cds}{"non_classe"} = 0;	

		$stats{$cds}{"nb_read_classes"}=0;

		#print $cds . "\n";
		$deja_vu{$cds}++;
	}
	my $read_name = $sam_infos[0];

	#print "\nOn traite le read : " . $read_name . " placé sur le cds : " . $cds . "\n";
	#
	#Verifier si le cds sur lequel est mappé le read est bien dans le hash hoemo_SNPs
	if (! exists($homeo_SNPs{$cds}))
	{
		#$count_reads{$cds}{"non_classe"}++;
		next;
	}
	
	#Calcul de l'intervalle de position que couvre le read grâce au cigar string et à la position de début de match du read
	#plus cigar string avec insertions d'étoiles s'il il a des deletion (D) ou des bases "sautées" (N)
	my $positon_debut_read = $sam_infos[3];
	my $position_fin_read = $positon_debut_read-1;
	my $cigar_string = $sam_infos[5];

	#print "Voila la cigar string d'origine : " . $cigar_string . "\n";

	my @cigar_list;
	my $read_string = $sam_infos[9];
	my $read_string_new="";
	my $compteur=0;
	while ($cigar_string =~ m/(\d+)(\D{1})/) #Un nombre \d+ et une non numérique \D (lettre)
	{
		#Append to the cigar list 
		my $nombre = $1;
		my $lettre = $2;
		my $cigar_part = $nombre . ":" . $lettre;
		push(@cigar_list,$cigar_part);
		my $number = length($1)+length($2);
		$cigar_string = substr($cigar_string,$number);

		#Increment position fin read
		if ($lettre eq "M" or $lettre eq "D" or $lettre eq "N")
		{
			$position_fin_read+= $nombre;
		}

		#String read new
		if ($lettre eq "D" or $lettre eq "N")
		{
			$read_string_new = $read_string_new . "*"x$nombre;
		}
		elsif ($lettre ne "H" and $lettre ne "P")
		{
			$read_string_new = $read_string_new . substr($read_string,$compteur,$nombre);
			$compteur+=$nombre;
		}	
	}
	
	
	#print "Voila le cigar string new : " . $read_string_new . "\n";

	#print "Position debut read = " . $positon_debut_read . " ; Position fin read = " . $position_fin_read . "\n";

	
	#Récupérer les positions attribués du CDS en question dans le hash homeo_SNPs{$cds} qui sont sur le read
	my @positions_attribuables;
	my $ref_pos = $homeo_SNPs{$cds};
	my %positions = %$ref_pos;
	foreach my $pos(sort ( {$a <=> $b } keys(%positions)))
	{
		if ($pos >= $positon_debut_read and $pos <= $position_fin_read)
		{
			push(@positions_attribuables,$pos);
		}
	}
	if (scalar(@positions_attribuables) eq 0)
	{
		$count_reads{$cds}{"non_classe"}++;
		next;
	}
	
	#Pour chaque position pouvant être attribuée, trouver sa position au sein du read pour ensuite trouver le nucleotide voulu
	my $parent1_attribution_count=0;
	my $parent2_attribution_count=0;
	foreach my $position(@positions_attribuables)
	{
		#print "Position attribuable : " . $position . "\n";

		my $nb_match_theorique = $position - $positon_debut_read + 1;
		my $penality_number=0; #For insertion (I) or soft clipping (S) when at the beginning
		my $increment_number=0; #When the ratio increment number / nb match theorique >= 1, we can stop look over the read
		 
		foreach my $cigar_section(@cigar_list)
		{
			my @cigar_section_infos = split(/:/,$cigar_section);
			my $cigar_section_number = $cigar_section_infos[0];
			my $cigar_section_letter = $cigar_section_infos[1];
			
			if ($cigar_section_letter eq "M" or $cigar_section_letter eq "N" or $cigar_section_letter eq "D")
			{
				$increment_number+=$cigar_section_number;
			}
			elsif ($cigar_section_letter eq "I" or $cigar_section_letter eq "S")
			{
				$penality_number+=$cigar_section_number;
			}

			last if ( ($increment_number/$nb_match_theorique) >= 1);
		}
		
		my $place_on_read = $nb_match_theorique + $penality_number - 1;

		#print "Position au sein du read : (en base Perl) : " . $place_on_read . "\n";


		my $nucleotide = substr($read_string_new,$place_on_read,1);
		

		#print "Nucleotide trouve : " . $nucleotide . "\n";

		foreach my $val(values @{$homeo_SNPs{$cds}{$position}{"parent1_attribution"}})
		{
			if ($val eq $nucleotide ) 
			{
				$parent1_attribution_count+=1;
			}
		}
	
		foreach my $val(values @{$homeo_SNPs{$cds}{$position}{"parent2_attribution"}})
		{
			if ($val eq $nucleotide )
			{
				$parent2_attribution_count+=1;
			}
		}
		
	}
	
	my $difference = $parent1_attribution_count - $parent2_attribution_count;
	
	my $nb_snp_found_read = $parent1_attribution_count + $parent2_attribution_count; #nombre de snp trouvé pour ce read
	next if ($nb_snp_found_read < $min_snp_read ); #si le nombre de snp sur le read n'est pas au minimun celui attendu on pass au reads suivant
	
	
	#print "\n===> Pour ce read : " . $read_name . "parent1 total count attribution = " . $parent1_attribution_count . " eug total count attribution = " . $parent2_attribution_count . "\n";
	if ($parent1_attribution_count != 0 and $parent2_attribution_count != 0)######------modifier------########
	{
		$count_reads{$cds}{"ambigu"}++;
		$conflict_reads{$cds}{$read_name}++;
		$nombre_de_read_avec_conflit_total++;
	}
	elsif ($difference > 0)
	{
		$count_reads{$cds}{"$parent1"}++;
		$nombre_de_read_categoriser_total++;
		$stats{$cds}{"nb_read_classes"}++;
		#print "count_read cds parent1 +1\n";
	}
	elsif ($difference < 0)
	{
		$count_reads{$cds}{"$parent2"}++;
		$nombre_de_read_categoriser_total++;
		$stats{$cds}{"nb_read_classes"}++;
		#print "count_reads cds eug +1\n";
	}
	else
	{
		$count_reads{$cds}{"non_classe"}++;
	}
}



my $string_fin = "NB_READS_";
my $string_fin1 = $string_fin.$parent1;
my $string_fin2 = $string_fin.$parent2;

print OUTPUT "CDS\t%R1*RT\t%R2*RT\t$string_fin1\t$string_fin2\tNB_READ_AMBIGU\tNB_READ_CLASSE\tNB_READ_NON_CLASSE\tNB_READ_TOTAL\n";
foreach my $cds(sort( {$a cmp $b} keys(%count_reads)))
{
	my $RT = $count_reads{$cds}{"non_classe"}+ $count_reads{$cds}{"ambigu"}+$count_reads{$cds}{"$parent2"} + $count_reads{$cds}{"$parent1"}; 
	my $RC = $count_reads{$cds}{"$parent2"} + $count_reads{$cds}{"$parent1"};
	my $R1 =0;
	my $R2 =0;
	if ($RC != 0)
	{
		$R1 = ($count_reads{$cds}{"$parent1"}/$RC )* $RT; #reads P1 : R1 = %reads P1 * reads totaux
		$R2 = ($count_reads{$cds}{"$parent2"}/$RC )* $RT; #reads P2 : R2 = %reads P2 * reads totaux
	}
	#tableau de resultats pour l'echantillon
	print OUTPUT $cds . "\t" . $R1 . "\t" . $R2 . "\t" . $count_reads{$cds}{"$parent1"}. "\t" . $count_reads{$cds}{"$parent2"}. "\t" . $count_reads{$cds}{"ambigu"}."\t" . $RC."\t" .$count_reads{$cds}{"non_classe"}. "\t" . $RT . "\t" . "\n";

}

close(SAM);
close(OUTPUT);


