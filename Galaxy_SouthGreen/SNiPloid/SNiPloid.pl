#!/usr/bin/perl

use Getopt::Long;
use Switch;
use Tie::File;

	   #####################################################
	   #                                                   #
	 #	    @@@@  @   @   @   @@@@   @         @      @      #
	#	   @	  @@  @       @   @  @                @       #
	# 	    @@@   @ @ @   @   @@@@   @   @@@   @   @@@@       #
	#          @  @  @@   @   @      @  @   @  @  @   @       #
	 #	   @@@@   @   @   @   @      @   @@@   @   @@@       #
	   #                                                   #
	   #####################################################

###############################################################################################################
#
#		SNiPloid
#		Author : Marine PERALTA
#
###############################################################################################################
#
#		Galaxy Version  
#
###############################################################################################################

#___________________________________
# Samples names
#-----------------------------------
$polyploidName = "" ;
$polyploid2Name = "" ; #
$genome1Name = "" ;
$genome2Name = "" ;
#___________________________________
# VCF files
#-----------------------------------
$VCFpolyploid = "" ;
$VCFpolyploid2 = "" ; #
$VCFgenome1 = "" ;
$VCFgenome2 = "" ;
$merged_VCF = "" ; # Polyploid + Genome1 + Genome 2
#___________________________________
# Depth of Coverage File
#-----------------------------------
$DOCpolyploid = "" ;
$DOCpolyploid2 = "" ; #
$DOCgenome1 = "" ;
$DOCgenome2 = "" ;
$merged_DOC = "" ; # Polyploid + Genome1 + Genome 2
#___________________________________
# Depth for each sample
#-----------------------------------
$depthPolyploid = 0 ;
$depthPolyploid2 = 0 ; #
$depthGenome1 = 0 ;
$depthGenome2 = 0 ;
#___________________________________
# Output Files
#-----------------------------------
$SNP_csv 		= "SNP_tab.txt";
$SNP_html 		= "SNP_view.html";
$SNP_count 		= "SNP_synthesis_tab.html";
$SNP_count_csv 	= "SNP_synthesis_tab.txt";
#___________________________________
# Other parameters
#-----------------------------------
$enableLowQuality = 0 ;  #default value for enable quality SNP = only PASS SNP are considered
$ref = 0 ; # default parameter = extern

$filtre_ouPas = 0 ;
$value_filter_p1 = 0 ;
$value_filter_p2 = 0 ;

$REPimages = "img_sniploid/";

$poly_poly_analysis = 0 ;


my $usage = qq~
Basic usage

For comparison between a polyploid and its parental diploid genomes:

    $0 --vp <VCF_polyploid> --vg1 <VCF_diploid> --cpp <depth_polyploid> --cg1 <depth_diploid> --dp <min_depth_polyploid> --dg1 <min_depth_diploid> --ref 1
    
For comparison between 2 polyploids:

    $0 --vp <VCF_polyploid1> --vp2 <VCF_polyploid2> --cpp <depth_polyploid1> --cpp2 <depth_polyploid2> --dp <min_depth_polyploid1> --dp2 <min_depth_polyploid2>
  
Usage:$0 <args>
where <args> are:

    --vp           <VCF file for polyploid>
    --vp2          <VCF file for polyploid 2>
    --vg1          <VCF file for diploid genome 1>
    --vg2          <VCF file for diploid genome 2>
    
    --cpp          <Depth file for polyploid>
    --cpp2         <Depth file for polyploid 2>
    --cg1          <Depth file for diploid genome 1>
    --cg2          <Depth file for diploid genome 2>
    
    --dp           <Minimum read depth at a position to make a call for polyploid>
    --dp2          <Minimum read depth at a position to make a call for polyploid 2>
    --dg1          <Minimum read depth at a position to make a call for diploid genome 1>
    --dg2          <Minimum read depth at a position to make a call for diploid genome 2>
    
    --oc           <Output file name for SNP list in csv>
    --oh           <Output file name for SNP list in HTML>
    --ocs          <Output file name for SNP count per gene in csv>
    --ohs          <Output file name for SNP count per gene in HTML>
    
    --vfp1         <Minimul allele frequency to consider as variant for polyploid 1 (in %). Default: 0>
    --vfp2         <Minimul allele frequency to consider as variant for polyploid 2 (in %). Default: 0>
    
    --elq          <Enable low quality SNP tag. Default: 0>
    --gn2          <Specify a name for diploid genome 2>
    --ref          <The reference must be included in the analysis as diploid genome. Default: 0>
~;
$usage .= "\n";


=pod
Add option for "Heterozygosity"
 Enable "heterozygosity" for genome 1 (reference intern) - not necessary...
 Enable "heterozygosity" for genome 1 and genome 2 (reference extern)
=cut

GetOptions (
	# "pn=s"	=>	\$polyploidName,
	# "pn2=s" =>	\$polyploid2Name, #
	# "gn1=s"	=>	\$genome1Name,  
	"gn2=s"		=>	\$genome2Name,
	"vp=s"		=>	\$VCFpolyploid,
	"vp2=s"		=>	\$VCFpolyploid2, #
	"vg1=s"		=>	\$VCFgenome1,
	"vg2=s"		=>	\$VCFgenome2,
	"vm=s"		=>	\$merged_VCF,
	"cpp=s"		=>	\$DOCpolyploid,
	"cpp2=s"	=>	\$DOCpolyploid2, #
	"cg1=s"		=>	\$DOCgenome1,
	"cg2=s"		=>	\$DOCgenome2,
	"cm=s"		=>	\$merged_DOC,
	"dp=i"		=>	\$depthPolyploid,
	"dp2=i"		=>	\$depthPolyploid2, #
	"dg1=i"		=>	\$depthGenome1,
	"dg2=i"		=>	\$depthGenome2,
	"oc=s"		=>	\$SNP_csv,
	"oh=s"		=>	\$SNP_html,
	"ohs=s"		=>	\$SNP_count,
	"ocs=s"		=>	\$SNP_count_csv,
	"elq=i"		=>	\$enableLowQuality,
	"ref=i"		=>	\$ref,
	#"fop=i" 	=>	\$filtre_ouPas,
	"vfp1=i" 	=>	\$value_filter_p1,
	"vfp2=i" 	=>	\$value_filter_p2,
	"img=s" 	=>	\$REPimages
# h = i	= >	\ $heterozygosity ,
);


# Validation - Samples names


die $usage
  if ( (!$VCFgenome1 || !$DOCgenome1 )  && (!$VCFpolyploid  || !$DOCpolyploid)   || (!$VCFpolyploid2 || !$DOCpolyploid2 )  && (!$VCFpolyploid  || !$DOCpolyploid));



%intervalle1 ;
%intervalle2 ;
%snp = () ;
my %snp_final ;
my %five = () ;
my %phased_regions = () ;

$nbTotGenes = 0 ;
$nbTotGenesVal = 0 ;
$nbTotGenesAna = 0 ;


if ($VCFpolyploid2 ne "") {
	$poly_poly_analysis = 1 ;
}

# if ($polyploidName eq "") {
	# print STDOUT "*** /!\\ ERROR: Missing name for polyploid - You have to specify a name for the polyploid species [--pn \"polyploid_name\"] $!" ;
	# die ("*** /!\\ ERROR: Missing name for polyploid - You have to specify a name for the polyploid species [--pn \"polyploid_name\"] $!") ;
# }

if ($poly_poly_analysis == 1) {
	print STDOUT "\nAnalysis Type: Polyploid vs Polyploid\n---------------------------------------";
	# print STDOUT "\nPolyploid 1: ".$polyploidName ;
	# print STDOUT "\nPolyploid 2:".$polyploid2Name ;
	# if ($polyploid2Name eq "") {
		# print STDOUT "*** /!\\ ERROR: Missing name for polyploid 2 - You have to specify a name for the polyploid species 2 [--pn2 \"polyploid_2_name\"] $!" ;
		# die ("*** /!\\ ERROR: Missing name for polyploid - You have to specify a name for the polyploid species 2 [--pn \"polyploid_2_name\"] $!") ;
	# }
}
else {
	print STDOUT "\nAnalysis Type: Polyploid vs Parental Genomes\n---------------------------------------";
	# print STDOUT "\nPolyploid: ".$polyploidName ;
	# print STDOUT "\nGenome 1: ".$genome1Name ;
	# print STDOUT "\nGenome 2: ".$genome2Name ;
	# if ($genome1Name eq "") {
		# die ("*** /!\\ ERROR: Missing name for genome 1 - You have to specify a name for the genome 1 species") ;
	# }
	# if ($genome2Name eq "") {
		# die ("*** /!\\ ERROR: Missing name for genome 2 - You have to specify a name for the genome 2 species") ;
	# }
	# Validation - depth
	if ($depthPolyploid == 0) {	
		die ("*** /!\\ ERROR: Missing depth information for polyploid");
	}
	if ($depthGenome1 == 0) {
		die ("*** /!\\ ERROR: Missing depth information for genome 1");
	}
	if ($ref == 0 && $depthGenome2 == 0) {
		die ("*** /!\\ ERROR: Missing depth information for genome 2");
	}
}




$time = time ;


################################################################
# 1) Polyploid vs Polyploid analysis
################################################################
if ($poly_poly_analysis == 1) {
	#print STDOUT "\n PASS";
	&Intervall_part1($DOCpolyploid) ;
	&Intervall_part2($DOCpolyploid2,$depthPolyploid2) ;
	
	&VCF_Analysis($VCFpolyploid);
	&VCF_Analysis($VCFpolyploid2);
	# CSS, titles, img, etc.
	&intro_output ;
	&poly_poly_output ;
}
################################################################
# 2) Polyploid vs Parental Diploid Genomes Analysis
################################################################
else {

# PART 1 : CREATING COMMON INTERVALS 
	
	&Intervall_part1($DOCpolyploid) ;
	&Intervall_part2($DOCgenome1,$depthGenome1) ;
	if ($ref == 0) { # genome2 => no parental genome as reference
		&Intervall_part2($DOCgenome2,$depthGenome2) ;
	}

# PART 2 and 3 : CREATING SNP TAB AND OUTPUTS

# VCF_Analysis : Create SNP hash and phasing
	
	&VCF_Analysis($VCFpolyploid);
	if ($ref == 1) { # Reference = one of two parental genomes
		&VCF_Analysis($VCFgenome1);
		# CSS, titles, img, etc.
		&intro_output ;
		# SNP Comparison and display
		&int_output ;
	}	
	else { # Extern Reference
		&VCF_Analysis($VCFgenome1);
		&VCF_Analysis($VCFgenome2);
		# CSS, titles, img, etc.
		&intro_output ;
		# SNP Comparison and display
		&ext_output ;
	}		
}




sub Intervall_part1 {
	my(@args) = @_;
	#print STDOUT "\nTEST ::: ".$args[0] ;
	open (TABSNP, $args[0]) or die ("Pbm a l'ouverture du fichier : $args[0]");
		@DOC = <TABSNP> ;
	close TABSNP ;

	$rec = 0 ;
	$position_pre ;
	$val_deb = "";
	$val_fin = "";
	$name_pre = "";
	
		
	foreach $line(@DOC) {
		if ($line ne $DOC[0]) {
			@ligne = split(/\t/ , $line);
			@position = split(/:/ , $ligne[0]);
			$name_gene = $position[0] ;

			if ($merged == 0) { # 1st File	 - Polyploid
				$depthcov = $ligne[1] ;
				if ($name_gene){
					if ($name_gene ne $gene_pre)
					{
						if ($rec == 1) {
							$position_fin = $position_pre ;
							$val_fin = $val_deb.$position_fin ; # Intervalle end position
							$intervalle1{$gene_pre}{$val_fin} = "ok" ;
						}
						$rec = 0;
					}	
					if ($depthcov >= $depthPolyploid){
						if ($rec == 0) {
							$position_deb = $position[1] ;
							$val_deb = $position_deb."-"; # Intervalle start position
						}
						$rec = 1 ;
					}
					if ($depthcov < $depthPolyploid){
						if ($rec == 1) {
							$position_fin = $position_pre ;
							$val_fin = $val_deb.$position_fin ; # Intervalle end position
							$intervalle1{$gene_pre}{$val_fin} = "ok" ;
						}
						$rec = 0;
					}
					
				}
			}
			else { # Merged files (2 or 3 species)
				if ($ref == 0) { # 3 species
					$depthcov1 = $ligne[$indiceGenome2] ;
					$depthcov2 = $ligne[$indicePolyploid1] ;
					$depthcov3 = $ligne[$indiceGenome1] ;				
					if ($name_gene){
						if (($depthcov1 >= $depthGenome2) && ($depthcov2 >= $depthPolyploid)&& ($depthcov3 >= $depthGenome1)){
							if ($rec == 0) {
								$position_deb = $position[1] ;
								$val_deb = $position_deb."-";
							}
							$rec = 1 ;
						}
						if (($depthcov1 < $depthGenome2) || ($depthcov2 < $depthPolyploid) || ($depthcov3 < $depthGenome1)){
							if ($rec == 1) {
								$position_fin = $position_pre ;
								$val_fin = $val_deb.$position_fin ;
								$intervalle1{$gene_pre}{$val_fin} = "ok" ;
							}
							$rec = 0 ;
						}
					}
				}
				else { # 2 species
					$depthcov1 = $ligne[$indicePolyploid1] ;
					$depthcov2 = $ligne[$indiceGenome1] ;
					if ($name_gene){
						if (($depthcov1 >= $depthPolyploid) && ($depthcov2 >= $depthGenome1)){
							if ($rec == 0) {
								$position_deb = $position[1] ;
								$val_deb = $position_deb."-";
							}
							$rec = 1 ;
						}
						if (($depthcov1 < $depthPolyploid) || ($depthcov2 < $depthGenome1)){
							if ($rec == 1) {
								$position_fin = $position_pre ;
								$val_fin = $val_deb.$position_fin ;
								$intervalle1{$gene_pre}{$val_fin} = "ok" ;
							}
							$rec = 0 ;
						}
					}
				}
			}	
			$position_pre = $position[1] ;
			$gene_pre = $name_gene ;
		}		
	}
	return (%intervalle1) ;
	
}
sub Intervall_part2 {
	
	my(@args) = @_;
	#print "\nintervall part 2 : $args[1]";
	
	open (TABSNP, $args[0]) or die ("Pbm a l'ouverture du fichier : $args[0]");
	#print STDOUT "\n$args[0]";
	@DOC = <TABSNP> ;
	my %tab ;
	foreach $li(@DOC) {
		if ($li =~ /^(.+):(.+)\t(.+)\t.+\t.+$/) {
   			$tab{$1}{$2} = $3;
		}
	}
	close TABSNP ;

	
	
	$rec = 0 ;
	$position_pre ;
	$val_deb = "";
	$val_fin = "";

	foreach my $interval(sort (keys(%intervalle1))){ 
		
		my $ref = $intervalle1{$interval};
		my %intervalls = %$ref;
		$name_gene = $interval ;
	
		foreach my $intervall(sort (keys(%intervalls))){ 
			$final = 2 ;
			$rec = 0 ;
			($debut,$fin) = split(/-/,$intervall);
			for ($i=$debut; $i <=$fin; $i++) {
				if ($tab{$interval}{$i} >= $args[1]){
					if ($rec == 0) {
						$position_deb = $i ;
						$val_deb = $position_deb."-";
					}
					$rec = 1 ;
					$final = 0 ;
				}
				if ($tab{$interval}{$i} < $args[1]){
					$final = 1 ;
					if ($rec == 1) {
						$position_fin = $i-1 ;
						$val_fin = $val_deb.$position_fin ;
						$intervalle2{$name_gene}{$val_fin} = "ok" ;
					}
					$rec = 0 ;
				}
			}
			if ($final == 0) {
				$val_fin = $val_deb.$fin ;
				$intervalle2{$name_gene}{$val_fin} = "ok" ;
			}
		}
	}
	if ($VCFgenome2 ne ""){
		%intervalle1 = %intervalle2 ;
	}
	
	foreach my $interval(sort (keys(%intervalle2))){
		my $ref = $intervalle2{$interval};
		my %intervalls = %$ref;
		$name_gene = $interval ;
	}	
	return (%intervalle2) ;
}

sub VCF_Analysis {
	%snp_final = () ;
	$compt_phasing = 0 ;
	$compt_five = 0 ;
	my(@args) = @_;
	
	open (TABSNP, "$args[0]") or die ("ERROR : file $args[0] don't exists");
		@VCF = <TABSNP> ;
	close TABSNP ;
	
	###########################################
	# test if VCF was filtered
	###########################################
	my $vcf_file = $args[0];
	my $grep_pass = `grep -c 'PASS' $vcf_file`;
	chomp($grep_pass);
	my $pass = "PASS";
	if (defined $grep_pass && $grep_pass == 0)
	{
		$pass = ".";
	}
	
	#print "$pass $grep_pass $vcf_file\n";
	
	foreach $line(@VCF){
		if ($line =~ /^#CHROM.+FORMAT\t(.+)$/) {
			$name_record = $1 ;
		}
		if ($line !~ /^#/){
			@infos_line = split(/\t/,$line) ;
			$gene = $infos_line[0];
			$position = $infos_line[1];
			$ref_allele = $infos_line[3];
			$alt_allele = $infos_line[4];
			
			if ($ref_allele =~/\w\w/ or $alt_allele =~/\w\w/)
			{
				next;
			}
			
			$snp_code = "[$ref_allele/$alt_allele]";
			$quality_of_snp = $infos_line[6];
			$depth_recuperation = $infos_line[7];
			$alleles = $infos_line[9];
		
			($GT,$AD,$FDP,$GQ,$PL) = split(":",$alleles);
			
			
			# PHASING
			
			if (($GT =~ /\|/) && ($previous_GT =~ /\//)) { # initialisation région
				$compt_phasing ++ ;
				$phased_regions{$gene}{$compt_phasing}{$previous_position} = $previous_GT ;
				$phased_regions{$gene}{$compt_phasing}{$position} = $GT ;
			}
			if (($GT =~ /\|/) && ($previous_GT =~ /\|/)) { # extension région
				$phased_regions{$gene}{$compt_phasing}{$position} = $GT ;
			}
			
			
			# $FDP = Filtered Depth
			# $DP = Total Depth
			
			my $DP;
			my @tags = split(";",$depth_recuperation);
			foreach my $tag(@tags)
			{
				if ($tag =~/DP=/)
				{
					$DP = $tag;
				}
			}

			#($sub1,$sub2) = split(",",$AD);
			#$somme = $sub1 + $sub2 ;
			
			$somme = 0;
			my @depth_of_alleles = split(",",$AD);
			my $sub1 = $depth_of_alleles[0];
			foreach my $depth_of_allele(@depth_of_alleles)
			{
				$somme += $depth_of_allele;
			}
			
			

			if ($somme == 0 ) {
				print STDOUT "ERROR : Cannot calculate ratio for ".$gene." [pos:".$position."]\n\"".$line."\"";
				die ("ERROR : Cannot calculate ratio for ".$gene." [pos:".$position."]\n\"".$line."\"");
			}
			else {
				$ratio = ($sub1/$somme)*100;
				$ratio = sprintf("%.0f", $ratio);
			}
			
			@DP = split ("=",$DP) ;

			$test_inside_interval = 0 ;
		
			my $ref = $intervalle2{$gene};
			my %hash = %$ref;
			
			foreach my $interval(keys(%hash)){
				my @pos = split(/-/,$interval) ;
				if ($position >= $pos[0] && $position <= $pos[1]) {
					$test_inside_interval = 1 ;
					last ;
				}
			}
			# ENABLE LOW_QUALITY SNP
			if ($enableLowQuality == 1) {
				if ($test_inside_interval == 1 ){ #
					if ($args[0] eq $VCFpolyploid) { # Polyploid
						$polyploidName = $name_record ;
						$snp{$gene}{$position} = $snp_code."\t".$AD."\t".$GT."\t".$DP[1]."-".$FDP ;
					}
					else { 
						if ($args[0] eq $VCFgenome1)  { # genome1
							$genome1Name = $name_record ;
							if (exists $snp{$gene}{$position}) { # if polyploid SNP
								$snp{$gene}{$position} = $snp{$gene}{$position}."\t".$snp_code."\t".$GT ;
								
								($code_snp,$ratio,$GT_poly,$DP_P,$code_G1,$GT_G1) = split(/\t/,$snp{$gene}{$position});
								@recupAlleles = split(/\[/,$code_snp);
								@recupAlleles = split(/\]/,$recupAlleles[1]);
								($alRef,$alAltP) = split(/\//,$recupAlleles[0]);
								@recupAlleles = split(/\[/,$code_G1);
								@recupAlleles = split(/\]/,$recupAlleles[1]);
								($alRef,$code_G1) = split(/\//,$recupAlleles[0]);
								
								#print "\nINFOS\n".$GT_poly."\t";
								#print $GT_G1."\t";
								#print $code_G1."\t";
								#print $alAltP."\n";
								if ((($GT_poly =~ /^0.1$/)||($GT_poly =~ /^1.0$/)) && (($GT_G1 =~ /^1.1$/)) && ($code_G1 eq  $alAltP)) {
									$five{$gene}{$position} = $GT_poly ; 
								}	
							}
							else { # if no polyploid SNP, key is empty
								$snp{$gene}{$position} = $ref_allele."\t\t\t\t".$snp_code."\t".$GT ;
							}
						}
						else { # genome2
							if ($args[0] eq $VCFgenome2) {
								$genome2Name = $name_record ;
								if (exists $snp{$gene}{$position}) { # if polyploid SNP
									$snp{$gene}{$position} = $snp{$gene}{$position}."\t".$snp_code."\t".$GT;
								}
								else { # if no polyploid SNP and no genome1, key is empty
									$snp{$gene}{$position} = $ref_allele."\t\t\t\t".$ref_allele."\t\t".$snp_code."\t".$GT ;
								}
							}
						}
						if ($args[0] eq $VCFpolyploid2) { # polyploid2
							$polyploid2Name = $name_record ;
							if (exists $snp{$gene}{$position}) { # if polyploid SNP
								$snp{$gene}{$position} = $snp{$gene}{$position}."\t".$snp_code."\t".$AD."\t".$GT."\t".$DP[1]."-".$FDP ;
							}
							else { # if no polyploid SNP, key is empty
								$snp{$gene}{$position} = $ref_allele."\t\t\t\t".$snp_code."\t".$AD."\t".$GT."\t".$DP[1]."-".$FDP ;
							}
						}
					}
				}
			}
			# ONLY PASS SNP CONSIDERED
			else {
				if (($test_inside_interval == 1 ) && ($quality_of_snp eq $pass) && ($snp{$gene}{$position} ne "LQ")){ #
					if ($args[0] eq $VCFpolyploid) { # Polyploid
						$polyploidName = $name_record ;
						$snp{$gene}{$position} = $snp_code."\t".$AD."\t".$GT."\t".$DP[1]."-".$FDP ;
					}
					else { 
						if ($args[0] eq $VCFgenome1)  { # genome1
							$genome1Name = $name_record ;
							if (exists $snp{$gene}{$position}) { # if polyploid SNP
								$snp{$gene}{$position} = $snp{$gene}{$position}."\t".$snp_code."\t".$GT ;
								
								($code_snp,$ratio,$GT_poly,$DP_P,$code_G1,$GT_G1) = split(/\t/,$snp{$gene}{$position});
								@recupAlleles = split(/\[/,$code_snp);
								@recupAlleles = split(/\]/,$recupAlleles[1]);
								($alRef,$alAltP) = split(/\//,$recupAlleles[0]);
								@recupAlleles = split(/\[/,$code_G1);
								@recupAlleles = split(/\]/,$recupAlleles[1]);
								($alRef,$code_G1) = split(/\//,$recupAlleles[0]);
								
								#print "\nINFOS\n".$GT_poly."\t";
								#print $GT_G1."\t";
								#print $code_G1."\t";
								#print $alAltP."\n";
								if ((($GT_poly =~ /^0.1$/)||($GT_poly =~ /^1.0$/)) && (($GT_G1 =~ /^1.1$/)) && ($code_G1 eq  $alAltP)) {
									$five{$gene}{$position} = $GT_poly ; 
								}	
							}
							else { # if no polyploid SNP, key is empty
								$snp{$gene}{$position} = $ref_allele."\t\t\t\t".$snp_code."\t".$GT ;
							}
						}
						else { # genome2
							if ($args[0] eq $VCFgenome2) {
								$genome2Name = $name_record ;
								if (exists $snp{$gene}{$position}) { # if polyploid SNP
									$snp{$gene}{$position} = $snp{$gene}{$position}."\t".$snp_code."\t".$GT;
								}
								else { # if no polyploid SNP and no genome1, key is empty
									$snp{$gene}{$position} = $ref_allele."\t\t\t\t".$ref_allele."\t\t".$snp_code."\t".$GT ;
								}
							}
						}
						if ($args[0] eq $VCFpolyploid2) { # polyploid2
							$polyploid2Name = $name_record ;
							if (exists $snp{$gene}{$position}) { # if polyploid SNP
								$snp{$gene}{$position} = $snp{$gene}{$position}."\t".$snp_code."\t".$AD."\t".$GT."\t".$DP[1]."-".$FDP ;
							}
							else { # if no polyploid SNP, key is empty
								$snp{$gene}{$position} = $ref_allele."\t\t\t\t".$snp_code."\t".$AD."\t".$GT."\t".$DP[1]."-".$FDP ;
							}
						}
					}
				}
				else {
					if ($quality_of_snp ne $pass) {  
						$snp{$gene}{$position} = "LQ";
					}
				}
			}
			################################################################################################################################	
		}
		$previous_GT = $GT ;
		$previous_position = $position ;
	}
	foreach my $s(sort(keys(%snp))){
		my $ref = $snp{$s};
		my %hash = %$ref;
		foreach my $snip(keys(%hash)){
			if ($snp{$s}{$snip} ne "LQ"){
				$snp_final{$s}{$snip} = $snp{$s}{$snip} ;
			}
		}
	}		
	return (%snp_final) ;
}

sub intro_output {

###########################################################
#			ANALYSE - CREATION FICHIERS DE SORTIE		  #
###########################################################

# Ouverture des fichiers
open (HTMLSNP, ">$SNP_html");
open (TABSNP, ">$SNP_csv");
open (HTMLCOUNT, ">$SNP_count");
open (TABCOUNT, ">$SNP_count_csv");

print HMTL "<html>\n";
print HTMLCOUNT "<html>\n";

print HTMLSNP "<head>\n";
print HTMLCOUNT "<head>\n";

#####################################################
#                       CSS                         #
#####################################################
print HTMLSNP "<style type=\"text/css\">\n";

print HTMLSNP "th {\n";
print HTMLSNP " border-color:black;\n";
print HTMLSNP " border-style:solid; \n";
print HTMLSNP " border-width:3px;\n";
print HTMLSNP " font-family: calibri;\n";
print HTMLSNP " }\n";

print HTMLSNP "body {text-align:center;}\n";

print HTMLSNP "table {\n";
print HTMLSNP " border-color:black;\n";
print HTMLSNP " margin:auto;\n";
print HTMLSNP " border-collapse: collapse;\n";
print HTMLSNP " border-width:3px; \n";
print HTMLSNP " border-style:solid; \n";
print HTMLSNP " }\n";

print HTMLSNP ".bord1 { \n";

print HTMLSNP " font-size: 11pt;\n";
print HTMLSNP " font-family: calibri;\n";
print HTMLSNP " border-width:1px;\n";
print HTMLSNP " border-top:3px;\n";
print HTMLSNP " border-left:3px;\n";
print HTMLSNP " border-right:3px;\n";
print HTMLSNP " border-style:solid; \n";
print HTMLSNP " border-color:black;\n";
print HTMLSNP " background-color : #c6c3bd; \n";
print HTMLSNP " }\n";

print HTMLSNP ".bord2 { \n";

print HTMLSNP " font-size: 11pt;\n";
print HTMLSNP " font-family: calibri;\n";
print HTMLSNP " border-width:1px;\n";
print HTMLSNP " border-top:3px;\n";
print HTMLSNP " border-left:3px;\n";
print HTMLSNP " border-right:3px;\n";
print HTMLSNP " border-style:solid; \n";
print HTMLSNP " border-color:black;\n";
print HTMLSNP " background-color : #c6c3ee; \n";
print HTMLSNP " }\n";


print HTMLSNP "td { \n";
print HTMLSNP " border-color:black;\n";
print HTMLSNP " }\n";

print HTMLSNP ".tdm { \n";
print HTMLSNP " border-color:black;\n";
print HTMLSNP " border-left:3px;\n";
print HTMLSNP " }\n";


print HTMLSNP ".td1 { \n";
print HTMLSNP " border-color:black;\n";
print HTMLSNP " font-size: 11pt;\n";
print HTMLSNP " font-family: calibri;\n";
print HTMLSNP " border-width:1px;\n";
print HTMLSNP " border-left:3px;\n";
print HTMLSNP " border-right:3px;\n";
print HTMLSNP " border-style:solid; \n";
print HTMLSNP " background-color : #c6c3bd; \n";
print HTMLSNP " }\n";

print HTMLSNP ".td2 { \n";
print HTMLSNP " border-color:black;\n";
print HTMLSNP " font-size: 11pt;\n";
print HTMLSNP " font-family: calibri;\n";
print HTMLSNP " border-width:1px;\n";
print HTMLSNP " border-left:3px;\n";
print HTMLSNP " border-right:3px;\n";
print HTMLSNP " border-style:solid; \n";
print HTMLSNP " background-color : #c6c3ee; \n";
print HTMLSNP " }\n";

print HTMLSNP ".ted { \n";
print HTMLSNP " border-color:black;\n";
print HTMLSNP " font-weight : bold;\n";
print HTMLSNP " background-color : #A19EED; \n";
print HTMLSNP " }\n";

print HTMLSNP ".ted2 { \n";
print HTMLSNP " border-color:black;\n";
print HTMLSNP " font-weight : bold;\n";
print HTMLSNP " background-color : #9A9D7C; \n";
print HTMLSNP " }\n";

print HTMLSNP ".tedG { \n";
print HTMLSNP " border-left:3px;\n";
print HTMLSNP " border-style:solid; \n";
print HTMLSNP " border-color:black;\n";
print HTMLSNP " font-weight : bold;\n";
print HTMLSNP " background-color : #A19EED; \n";
print HTMLSNP " }\n";

print HTMLSNP ".tedG2 { \n";
print HTMLSNP " border-left:3px;\n";
print HTMLSNP " border-style:solid; \n";
print HTMLSNP " border-color:black;\n";
print HTMLSNP " font-weight : bold;\n";
print HTMLSNP " background-color : #9A9D7C; \n";
print HTMLSNP " }\n";

print HTMLSNP ".final { \n";
print HTMLSNP " border-left:3px;\n";
print HTMLSNP " border-right:0px;\n";
print HTMLSNP " border-top:0px;\n";
print HTMLSNP " border-bottom:0px;\n";
print HTMLSNP " border-style:solid; \n";
print HTMLSNP " border-color:black;\n";
print HTMLSNP " background-color : white; \n";
print HTMLSNP " }\n";

print HTMLSNP ".auto-style1 {";
print HTMLSNP "	font-weight: normal;";
print HTMLSNP "	font-size: x-small;";
print HTMLSNP "}";


print HTMLSNP "</style>\n";

print HTMLCOUNT "<style type=\"text/css\">\n";

print HTMLCOUNT "th {\n";
print HTMLCOUNT " border-style:solid; \n";
print HTMLCOUNT " border-color:black;\n";
print HTMLCOUNT " border-width:3px;\n";
print HTMLCOUNT " font-family:calibri;\n";
print HTMLCOUNT " }\n";

print HTMLCOUNT "table {\n";
print HTMLCOUNT " margin:auto;\n";
print HTMLCOUNT " border-collapse: collapse;\n";
print HTMLCOUNT " border-width:3px; \n";
print HTMLCOUNT " border-style:solid; \n";
print HTMLCOUNT " border-color:black;\n";
print HTMLCOUNT " }\n";

print HTMLCOUNT ".th {\n";
print HTMLCOUNT " font-weight : normal;\n";
print HTMLCOUNT " border-style:solid; \n";
print HTMLCOUNT " border-color:white;\n";
print HTMLCOUNT " border-width:0px;\n";
print HTMLCOUNT " font-family:consolas;\n";
print HTMLCOUNT " }\n";

print HTMLCOUNT ".tab2 {\n";
print HTMLCOUNT " margin:auto;\n";
print HTMLCOUNT " border-collapse: collapse;\n";
print HTMLCOUNT " border-style:solid; \n";
print HTMLCOUNT " border-width:3px; \n";
print HTMLCOUNT " border-color:white;\n";
print HTMLCOUNT " }\n";

print HTMLCOUNT ".tab {\n";
print HTMLCOUNT " margin:auto;\n";
print HTMLCOUNT " border-collapse: collapse;\n";
print HTMLCOUNT " border-width:3px;\n ";
print HTMLCOUNT " border-style:solid;\n ";
print HTMLCOUNT " border-color:black;\n";
print HTMLCOUNT " }\n";

print HTMLCOUNT ".td1 { \n";
print HTMLCOUNT " border-color:black;\n";
print HTMLCOUNT " font-size: 11pt;\n";
print HTMLCOUNT " font-family: calibri;\n";
print HTMLCOUNT " border-width:1px;\n";
print HTMLCOUNT " border-left:3px;\n";
print HTMLCOUNT " border-right:3px;\n";
print HTMLCOUNT " border-style:solid; \n";
print HTMLCOUNT " background-color : #c6c3bd; \n";
print HTMLCOUNT " }\n";

print HTMLCOUNT ".td2 { \n";
print HTMLCOUNT " border-color:black;\n";
print HTMLCOUNT " font-size: 11pt;\n";
print HTMLCOUNT " font-family: calibri;\n";
print HTMLCOUNT " border-width:1px;\n";
print HTMLCOUNT " border-left:3px;\n";
print HTMLCOUNT " border-right:3px;\n";
print HTMLCOUNT " border-style:solid; \n";
print HTMLCOUNT " background-color : #c6c3ee; \n";
print HTMLCOUNT " }\n";

print HTMLCOUNT ".td3 { \n";
print HTMLCOUNT " border-color:black;\n";
print HTMLCOUNT " font-size: 11pt;\n";
print HTMLCOUNT " font-weight: bold;\n";
print HTMLCOUNT " font-family: calibri;\n";
print HTMLCOUNT " border-width:3px;\n";
print HTMLCOUNT " border-left:3px;\n";
print HTMLCOUNT " border-right:3px;\n";
print HTMLCOUNT " border-style:solid; \n";
print HTMLCOUNT " background-color : white; \n";
print HTMLCOUNT " }\n";


print HTMLCOUNT ".ted { \n";
print HTMLCOUNT " border-color:black;\n";
print HTMLCOUNT " font-weight : bold;\n";
print HTMLCOUNT " background-color : #A19EED; \n";
print HTMLCOUNT " }\n";

print HTMLCOUNT ".ted2 { \n";
print HTMLCOUNT " border-color:black;\n";
print HTMLCOUNT " font-weight : bold;\n";
print HTMLCOUNT " background-color : #9A9D7C; \n";
print HTMLCOUNT " }\n";

print HTMLCOUNT ".ted3 { \n";
print HTMLCOUNT " border-color:black;\n";
print HTMLCOUNT " font-family: calibri;\n";
print HTMLCOUNT " color: white;\n";
print HTMLCOUNT " background-color : #333333; \n";
print HTMLCOUNT " }\n";

print HTMLCOUNT ".auto-style1 {";
print HTMLCOUNT "	font-weight: normal;";
print HTMLCOUNT "	font-size: x-small;";
print HTMLCOUNT "}";

print HTMLCOUNT "</style>\n";

###################################################################################################################################################################################

print HTMLSNP "</head>\n";

print HTMLSNP "<center><img src=\"".$REPimages."SNiPloid7.png\" WIDTH=250></center>";
if ($poly_poly_analysis == 0) {
	print HTMLSNP "<center><img src=\"".$REPimages."arbre.png\" WIDTH=400></center>";
}
print HTMLSNP "<p>\n";
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print HTMLCOUNT "</head>\n";
print HTMLCOUNT "<table><tr><td class=\"tab2\"><img src=\"".$REPimages."SNiPloid7.png\" WIDTH=250>";

print HTMLCOUNT "<h3><font face=\"calibri\">Synthesis of the analysis</font></h3></td>";
if ($poly_poly_analysis == 0) {
	
	print HTMLCOUNT "<td class=\"tab2\" width=50></td>";
	print HTMLCOUNT "<td class=\"tab2\"><table border=\"1\" border cellpadding=\"5\" style=\"text-align:center;\">";
	print HTMLCOUNT "<tr><th>Diploids</th><th>Polyploid</th><th>Identity</th><th>Interpretation</th></tr>";
	print HTMLCOUNT "<tr><td>[1/2]</td><td>[1]</td><td>!=</td><td><img src=\"".$REPimages."1.png\" height=30></td></tr>";
	print HTMLCOUNT "<tr><td>[1/2]</td><td>[2]</td><td>!=</td><td><img src=\"".$REPimages."2.png\" height=30></td></tr>";
	print HTMLCOUNT "<tr><td>[1]</td><td>[1/2]</td><td>!=</td><td><img src=\"".$REPimages."3.png\" height=30></td></tr>";
	print HTMLCOUNT "<tr><td>[2]</td><td>[1/2]</td><td>!=</td><td><img src=\"".$REPimages."4.png\" height=30></td></tr>";
	print HTMLCOUNT "<tr><td>[1/2]</td><td>[1/2]</td><td>=</td><td><img src=\"".$REPimages."5v.png\" height=30></td></tr>";
	print HTMLCOUNT "<tr><td>[1]</td><td>[2]</td><td>!=</td><td><img src=\"".$REPimages."other.png\" height=30></td></tr>";
	print HTMLCOUNT "<tr><td>[1]</td><td>[2/3]</td><td>!=</td><td><img src=\"".$REPimages."other.png\" height=30></td></tr>";
	print HTMLCOUNT "<tr><td>[1/2]</td><td>[2/3]</td><td>!=</td><td><img src=\"".$REPimages."other.png\" height=30></td></tr>";
	print HTMLCOUNT "<tr><td>[1/2]</td><td>[1/2]</td><td>!=</td><td><img src=\"".$REPimages."other.png\" height=30><img src=\"".$REPimages."HG1.png\" height=30></td></tr>";
	print HTMLCOUNT "<tr><td>[1]</td><td>[1/2]</td><td>!=</td><td><img src=\"".$REPimages."other.png\" height=30><img src=\"".$REPimages."HG1.png\" height=30></td></tr>";
	print HTMLCOUNT "</table></td>";
	print HTMLCOUNT "<td class=\"tab2\" width=50></td>";
	print HTMLCOUNT "<td class=\"tab2\"><center><img src=\"".$REPimages."arbre.png\" WIDTH=400></center></td>";
	#print HTMLCOUNT "<td><table border=\"1\" border cellpadding=\"5\" style=\"text-align:center;\"><tr><th>Diploids</th><th>Polyploid</th><th>Identity</th><th>Interpretation</th></tr></table></td>";
	print HTMLCOUNT "</tr></table>";
}
print HTMLCOUNT "<p>\n";
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

print HTMLSNP "<body>\n";

if ($poly_poly_analysis == 1) {
	print HTMLSNP "<center><h3><font face=\"calibri\">Result of SNP comparison of two Polyploids</font></h3></center>";
}
else {
	print HTMLSNP "<center><h3><font face=\"calibri\">Result of SNP comparison of a Polyploid and its Parental Genomes (Genome 1 and Genome 2 as reference)</font></h3></center>";
}

print HTMLSNP "<p>\n";

# COLUMNS - HTMLSNP SNP VIEW
print HTMLSNP "<table border=\"1\" border cellpadding=\"5\" style=\"text-align:center;\"> \n";
print HTMLSNP "<tr>\n";
print HTMLSNP "<th>Gene</th>"; 																						# (1) Gene
print HTMLSNP "<th>Position</th>"; 																					# (2) Position

if ($poly_poly_analysis == 1) {
	print HTMLSNP "<th>REF<br></th>";
	print HTMLSNP "<th>Polyploid 1<br><span class=\"auto-style1\">".$polyploidName."</span></th>"; 											# (3) Polyploid
	print HTMLSNP "<th>Polyploid 2<br><span class=\"auto-style1\">".$polyploid2Name."</span></th>"; 										# (4) Polyploid 2
	
	print HTMLSNP "<th>[Filtered/Total] Depth<br>Polyploid 1<br><span class=\"auto-style1\">".$polyploidName."</span></th>"; 			# (8) Filtered Depth
	print HTMLSNP "<th>[Filtered/Total] Depth<br>Polyploid 2<br><span class=\"auto-style1\">".$polyploid2Name."</span></th>"; 			# (9) Total Depth
	
	# Entête fichier SNP VIEW TAB
	print TABSNP "Gene\t";																							# (1) Gene
	print TABSNP "Position\t";																						# (2) Position
	print TABSNP "REF\t";																						# (2) Position
	print TABSNP "Polyploid 1: ".$polyploidName."\t";																	# (3) Polyploid
	print TABSNP "Polyploid 2: ".$genome1Name."\t";																	# (4) Genome 1
	print TABSNP "P1 Filtered Depth\t";																				# (9) Filtered Depth
	print TABSNP "P1 Total Depth\t";																					# (10) Total Depth
	print TABSNP "P2 Filtered Depth\t";																				# (9) Filtered Depth
	print TABSNP "P2 Total Depth\n";																					# (10) Total Depth
}
else {
	 																		# (3) Reference
	print HTMLSNP "<th>Polyploid<br><span class=\"auto-style1\">".$polyploidName."</span></th>"; 					# (3) Polyploid
	print HTMLSNP "<th>Genome 1<br><span class=\"auto-style1\">".$genome1Name."</span></th>"; 						# (4) Genome 1
	print HTMLSNP "<th>Genome 2<br><span class=\"auto-style1\">".$genome2Name."</span></th>"; 						# (5) Genome 2
	print HTMLSNP "<th>Validation</th>";
	print HTMLSNP "<th>Ratio (%)<br><span class=\"auto-style1\">".$genome2Name." : ".$genome1Name."</span></th>";	# (7) Ratio
	print HTMLSNP "<th>Filtered<br>Depth</th>"; 																	# (8) Filtered Depth
	print HTMLSNP "<th>Total<br>Depth</th>"; 																		# (9) Total Depth
	print HTMLSNP "<th>SNP Class</th>"; 																		# (9) Total Depth
	
	# Entête fichier SNP VIEW TAB
	print TABSNP "Gene\t";																							# (1) Gene
	print TABSNP "Position\t";																						# (2) Position
	print TABSNP "Polyploid: ".$polyploidName."\t";																	# (3) Polyploid
	print TABSNP "Genome 1: ".$genome1Name."\t";																	# (4) Genome 1
	print TABSNP "Genome 2: ".$genome2Name."\t";																	# (5) Genome 2
	print TABSNP "Validation\t";																					# (6) Validation
	print TABSNP "Ratio (%) ".$genome2Name."\t";																	# (7) Ratio Genome 1
	print TABSNP "Ratio (%) ".$genome1Name."\t";																	# (8) Ratio Genome 2
	print TABSNP "Filtered Depth\t";																				# (9) Filtered Depth
	print TABSNP "Total Depth\t";																					# (10) Total Depth
	print TABSNP "SNP Class\n";																					# (10) Total Depth
}
print HTMLSNP "</tr>\n";



#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print HTMLCOUNT "<body>\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";
print HTMLCOUNT "<p>\n";
# COLUMNS - HTML SYNTHESIS
print HTMLCOUNT "<table border=\"1\" border cellpadding=\"5\" style=\"text-align:center;\"> \n";
print HTMLCOUNT "<tr>\n";
print HTMLCOUNT "<th>Gene</th>";																								# (1) Gene
print HTMLCOUNT "<th>Interval Size<br>Analysed (pb)</th>";																		# (2) Interval Size Analyzed
print HTMLCOUNT "<th>nb Positions<br>with SNP</th>";	
if ($poly_poly_analysis == 0) {																									# (3) nB Positions With SNP
	print HTMLCOUNT "<th><img src=\"".$REPimages."1.png\" height=30></th>";														# (4) [1]
	print HTMLCOUNT "<th><img src=\"".$REPimages."2.png\" height=30></th>";														# (5) [2]
	print HTMLCOUNT "<th><img src=\"".$REPimages."3ou4.png\" height=30></th>";													# (6) [3 or 4]
	print HTMLCOUNT "<th><img src=\"".$REPimages."3.png\" height=30></th>";													# (6) [3 or 4]
	print HTMLCOUNT "<th><img src=\"".$REPimages."4.png\" height=30></th>";													# (6) [3 or 4]
	print HTMLCOUNT "<th><img src=\"".$REPimages."5v.png\" height=30></th>";														# (11) [5]
	print HTMLCOUNT "<th><img src=\"".$REPimages."other.png\" height=30></th>";																							# (7) [other]
	print HTMLCOUNT "<th><img src=\"".$REPimages."HG1.png\" height=30><br><span class=\"auto-style1\">".$genome1Name."</span></th>";	# (8) Heterozygosity For Genome 1	
	print HTMLCOUNT "<th>SNP Intra-Diploids <br><img src=\"".$REPimages."5v.png\" height=30> + <img src=\"".$REPimages."1.png\" height=30> + <img src=\"".$REPimages."2.png\" height=30> + <img src=\"".$REPimages."other.png\" height=30></th>";																					# (9) SNP Diploids
	print HTMLCOUNT "<th>SNP Intra-Polyploid <br> <img src=\"".$REPimages."5v.png\" height=30> + <img src=\"".$REPimages."3ou4.png\" height=30> + <img src=\"".$REPimages."other.png\" height=30></th>";																					# (10) SNP Polyploid
	print HTMLCOUNT "<th>Ratio (%)<br><span class=\"auto-style1\">".$genome2Name." : ".$genome1Name."</span></th>\n\n\n\n\n\n\n\n";				# (12) Ratio %
	
	print HTMLCOUNT "</tr>\n";
# Entête HTMLSNP Synthesis
	print TABCOUNT "Gene\t";																									# (1) Gene
	print TABCOUNT "Interval Size Analysed (pb)\t";																				# (2) Interval Size Analysed (pb)
	print TABCOUNT "nb SNP positions\t";																						# (3) nb SNP positions
	print TABCOUNT "1\t";																										# (4) [1]
	print TABCOUNT "2\t";																										# (5) [2]
	print TABCOUNT "3 or 4\t";																									# (6) [3 or 4]
	print TABCOUNT "3\t";																										# (6) [3 or 4]
	print TABCOUNT "4\t";																										# (6) [3 or 4]
	print TABCOUNT "5\t";																										# (11) [5]
	print TABCOUNT "other\t";																									# (7) [other]
	print TABCOUNT "SNP Heterozygosity Genome 1\t";																				# (8) SNP Heterozygosity Genome 1
	print TABCOUNT "SNP Diploids\t";																							# (9) SNP Diploids
	print TABCOUNT "SNP Polyploid\t";																							# (10) SNP Polyploid
	print TABCOUNT "Ratio (%) ".$genome2Name."\t";																				# (12) Ratio (%) Genome 2
	print TABCOUNT "Ratio (%) ".$genome1Name."\n";	 																			# (13) Ratio (%) Genome 1
}
else {
	print HTMLCOUNT "<th>P1 = P2<br><span style=\"font-weight: normal\"><span style=\"background:#DE8A8A\">[1/2]</span> vs <span style=\"background:#DE8A8A\">[1/2]</span></span></th>";					# (4) [1]
	print HTMLCOUNT "<th>P1 = P2<br><span style=\"font-weight: normal\"><span style=\"background:#5CAAD2\">[1]</span> vs <span style=\"background:#5CAAD2\">[1]</span></span></th>";							# (4) [1]
	print HTMLCOUNT "<th>SNP<br>interpolyploids<br>P1 &ne; P2<br><span style=\"auto-style1\"></span></th>";			# (4) [1] #DE8A8A
	print HTMLCOUNT "<th>P1 &ne; P2<br>2 Alleles<br><span style=\"font-weight: normal\"><span style=\"background:#5CAAD2\">[1]</span> vs <span style=\"background:#5CAAD2\">[2]</span></span></th>";		# (4) [1]
	print HTMLCOUNT "<th>P1 &ne; P2<br>2 Alleles<br><span style=\"font-weight: normal\"><span style=\"background:#DE8A8A\">[1/2]</span> vs <span style=\"background:#5CAAD2\">[1]</span> or <span style=\"background:#5CAAD2\">[2]</span></span></th>";		# (4) [1]
	print HTMLCOUNT "<th>P1 &ne; P2<br>3 Alleles<br><span style=\"font-weight: normal\"><span style=\"background:#DE8A8A\">[1/2]</span> vs <span style=\"background:#5CAAD2\">[3]</span></span></th>";						# (4) [1]
	print HTMLCOUNT "<th>P1 &ne; P2<br>3 Alleles<br><span style=\"font-weight: normal\"><span style=\"background:#DE8A8A\">[1/2]</span> vs <span style=\"background:#DE8A8A\">[1/3]</span></span></th>";				# (4) [1]
	print HTMLCOUNT "<th>SNP<br>intra P1<br><span class=\"auto-style1\">".$polyploidName."</span></th>";						# (4) [1]
	print HTMLCOUNT "<th>SNP<br>intra P2<br><span class=\"auto-style1\">".$polyploid2Name."</span></th>";						# (4) [1] 
	
	print HTMLCOUNT "</tr>\n";
		# Entête HTMLSNP Synthesis
	print TABCOUNT "Gene\t";									# (1) Gene
	print TABCOUNT "Interval Size Analysed (pb)\t";				# (2) Interval Size Analysed (pb)
	print TABCOUNT "nb positions with SNP\t";					# (3) nb SNP positions
	print TABCOUNT "P1 = P2 [1/2] vs [1/2]\t";					# (4) [1]
	print TABCOUNT "P1 = P2 [1] vs [1]\t";						# (4) [1]
	print TABCOUNT "SNP interpolyploids P1 diff P2\t";			# (4) [1] #DE8A8A
	print TABCOUNT "P1 diff P2 2 Alleles [1] vs [2]\t";			# (4) [1]
	print TABCOUNT "P1 diff P2 2 Alleles [1/2] vs [1] or [2]\t";	# (4) [1]
	print TABCOUNT "P1 diff P2 3 Alleles [1/2] vs [3]\t";			# (4) [1]
	print TABCOUNT "P1 diff P2 3 Alleles [1/2] vs [1/3]\t";		# (4) [1]
	print TABCOUNT "SNP intra P1 ".$polyploidName."\t" ;				# (4) [1]
	print TABCOUNT "SNP intra P2 ".$polyploid2Name."\n" ;			# (4) [1]

}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
$ligneInter = 0 ;
$totalSNP = 0 ;
$totalSize = 0 ;

$total11 = 0 ;
$total22 = 0 ;
$total3ou4 = 0 ;
$total5 = 0 ;
$total512 = 0 ;
$total534 = 0 ;
$totalOther = 0 ;
$totalGenome2 = 0 ;

$total3 = 0 ;
$total4 = 0 ;

$totalNbPolyploid1 = 0 ;				# SNP heterozygosity for P1
$totalNbPolyploid2 = 0 ;				# SNP heterozygosity for P2
$totalNbCommuns =  0 ;					# SNP heterozygosity [P1] = [P2]
$totalNbCommunsHomo = 0 ;			# SNP homozygosity [P1] = [P2]
$totalNbDifferent =  0 ;				# [P1] ne [P2]
$totalNbHomoDiff =  0 ;		# Example : P1 = [A] ; P2 = [G]
$totalNbAlleleCommun =  0 ;		# Example : P1 = [A/G] ; P2 = [A]
$totalAlleleDifferent =  0 ;	# Example : P1 = [A/G] ; P2 = [C] or [T]
$totalAlleleCommunH =  0 ;		# Example : P1 = [A/G] ; P2 = [A/C]
		
}

sub int_output { # intern reference
	
	foreach my $s(sort(keys(%snp_final))){
		#######################################
		if ($ligneInter == 0) {
			print HTMLSNP "<tr class=\"bord1\">\n";
			print HTMLCOUNT "<tr class=\"td1\ border-width =\"3px\">\n";
		}
		else {
			print HTMLSNP "<tr class=\"bord2\">\n";
			print HTMLCOUNT "<tr class=\"td2\">\n";
		}
		#######################################
	
		my $ref = $snp_final{$s};
		my %hash = %$ref;
		$taille = keys(%hash);
		
		# taille de la ligne gene
		print HTMLSNP "<td  rowspan=\"".$taille."\"><b>".$s."</b></td>";
		
		$ligneOK = 0 ;
		
		$nbPolyploid = 0;
		$nbGenomes = 0 ;
		$nbGenome2 = 0 ;  		# SNP chez les diploides dans un cas "Other" avec genome 1 heterozygote
		$nbCommuns = 0;
		$nbPoly_only = 0 ;
		$nbSub_only = 0 ;
		
		$case5 = 0 ;
		$case1 = 0 ;
		$case2 = 0 ;
		$case3ou4 = 0;
		$case3 = 0 ;
		$case4 = 0 ;
		$caseOther = 0; 
		$casePolyplother = 0 ;	# SNP chez le polyploide dans un cas "Other"
		$caseDiplother = 0 ; 	# SNP chez les diploides dans un cas "Other"
	
	
	
		# Moyenne Ponderee
		$moyenneSNPindep1 = 0 ;
		$moyenneSNPindep2 = 0 ;
		
		@tabTrie = sort ({ $a <=> $b }keys %hash);

		#foreach my $c(sort ({$hash{$a} <=> $hash{$b}} keys %hash)) {
		foreach my $c(@tabTrie) {
			if ($ligneOK == 1) {
				if ($ligneInter == 0) {
					print HTMLSNP "<tr class=\"td1\">\n";
				}
				else {
					print HTMLSNP "<tr class=\"td2\">\n";
				}
			}
			################################################################################
			### Recuperation des informations ###
			
			($code_snp,$ratio,$GT_poly,$DP_P,$code_G1,$GT_G1) = split(/\t/,$snp_final{$s}{$c});
			($DP_P, $FDP) = split(/-/, $DP_P);
			#print STDOUT "\n($code_snp:$GT_poly) - ($code_G1:$GT_G1)" ;
			if ($GT_poly ne "") { # Polyploide = [0.0] ou [0.1] ou [1.1]
				@recupAlleles = split(/\[/,$code_snp);
				@recupAlleles = split(/\]/,$recupAlleles[1]);
				($alRef,$alAltP) = split(/\//,$recupAlleles[0]);
				# Attribution des alleles au polyploide si pas de SNP
				if ($GT_poly =~ /^0.0$/) { $code_snp = $alRef ; }
				if ($GT_poly =~ /^1.1$/) { $code_snp = $alAltP ; }
				# Attribution des alleles au genome 1 si pas de SNP
				if (($GT_G1 eq "") || ($GT_G1 =~ /^0.0$/)) { $code_G1 = $alRef ; }
				if ($GT_G1 =~ /^1.1$/) {
					@recupAlleles = split(/\[/,$code_G1);
					@recupAlleles = split(/\]/,$recupAlleles[1]);
					($alRef,$alAlt) = split(/\//,$recupAlleles[0]);
					$code_G1 = $alAlt;
				}
			}
			elsif ($GT_G1 ne "") { # pas de SNP polyploide dans le fichier 1 (fichiers non mergés) -> equivalent de [0.0]
				@recupAlleles = split(/\[/,$code_G1);
				@recupAlleles = split(/\]/,$recupAlleles[1]);
				($alRef,$alAlt) = split(/\//,$recupAlleles[0]);
				# Attribution des Alleles au genome 1
				if ($GT_G1 =~ /^1.1$/) { $code_G1 = $alAlt ; }
				if ($GT_G1 =~ /^0.0$/) { $code_G1 = $alRef ; }
			}

			################################################################################
			$noSNPpoly = "ok" ;
			
			
			
			
			if ((($GT_poly =~ /^0.1$/)||($GT_poly =~ /^1.0$/)) && (($GT_G1 =~ /^1.1$/)) && ($code_G1 eq  $alAltP)) {
				if ($ligneInter == 0) {
					print HTMLSNP "<td class=\"tedG2\">".$c."</td>";
					print HTMLSNP "<td class=\"ted2\">".$code_snp."</td>";
					print HTMLSNP "<td class=\"ted2\">".$code_G1."</td>";
					print HTMLSNP "<td class=\"ted2\">".$alRef."</td>"; #REF
					print HTMLSNP "<td class=\"ted2\">OK</td>";
					
				}
				else {
					print HTMLSNP "<td class=\"tedG\">".$c."</td>";
					print HTMLSNP "<td class=\"ted\">".$code_snp."</td>";
					print HTMLSNP "<td class=\"ted\">".$code_G1."</td>";
					print HTMLSNP "<td class=\"ted\">".$alRef."</td>"; #REF
					print HTMLSNP "<td class=\"ted\">OK</td>";
				}
			}
			else {
				print HTMLSNP "<td style=\"border-left:3px solid black\">".$c."</td>";
				print HTMLSNP "<td>".$code_snp."</td>";
				print HTMLSNP "<td>".$code_G1."</td>";
				print HTMLSNP "<td>".$alRef."</td>"; #REF
				print HTMLSNP "<td>not OK</td>";
			}	
			print TABSNP $s."\t".$c."\t".$code_snp."\t".$code_G1."\t".$alRef."\t";		
	
	
	
			$tailleImg = 35 ;
			
			
			if (($GT_poly =~ /^0.1$/)||($GT_poly =~ /^1.0$/)) { #	SNP POLYPLOID - [0/1] [0|1] [1|0]
				# Moyenne du Ratio -----------------------------------
				($sub1,$sub2) = split(",",$ratio);
				$somme = $sub1 + $sub2 ;

				if ($somme == 0 ) {
					print STDOUT "ERROR : Cannot calculate ratio for ".$gene." [pos:".$position."]\n\"".$line."\"";
					die ("ERROR : Cannot calculate ratio for ".$gene." [pos:".$position."]\n\"".$line."\"");
				}
				else {
					$ratio = ($sub1/($sub1+$sub2))*100;
					$ratio = sprintf("%.0f", $ratio);
				}
				#$ratio = sprintf("%.0f", $ratio);
				$ratio2 = 100-$ratio ;

				#-----------------------------------------------------
				if (($GT_G1 =~ /^1.1$/)){ # Pas de SNP Genome1  [1/1] [1|1] SNP entre genome1 et genome2
					if ($code_G1 eq $alAltP) {							# 5
						$moyenneSNPindep1 = $moyenneSNPindep1 + $ratio ;
						$moyenneSNPindep2 = $moyenneSNPindep2 + $ratio2 ;
						if ($ligneInter == 0) {
							print HTMLSNP "<td class=\"ted2\">".$ratio.":".$ratio2."<br>";							# RATIO %
							print HTMLSNP "<img src=\"".$REPimages."r1.png\" height=10 width=".$ratio.">";			# IMG RATIO 1
							print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=10 width=".$ratio2."></td>";	# IMG RATIO 2
							print HTMLSNP "<td class=\"ted2\">".$FDP."</td>";
							print HTMLSNP "<td class=\"ted2\">".$DP_P."</td>";
							print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5v.png\" height=".$tailleImg."></td>";
							
						}
						else {
							print HTMLSNP "<td class=\"ted\">".$ratio.":".$ratio2."<br>";							# RATIO %
							print HTMLSNP "<img src=\"".$REPimages."r1.png\" height=10 width=".$ratio.">";			# IMG RATIO 1
							print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=10 width=".$ratio2.">";			# IMG RATIO 2
							print HTMLSNP "<td class=\"ted\">".$FDP."</td>";
							print HTMLSNP "<td class=\"ted\">".$DP_P."</td>";
							print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5v.png\" height=".$tailleImg."></td>";

						}
						print TABSNP "OK\t".$ratio."\t".$ratio2."\t".$FDP."\t".$DP_P."\t5";
						$case5 ++ ;
					}
					else { # Other 0.1 - 1.1 (O GA A)					# Other [SNP DIPLO + SNP POLY]
						print HTMLSNP "<td>#</td><td>#</td><td>#</td>";
						print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."other.png\" height=".$tailleImg."></td>";
						print TABSNP "not OK\t#\t#\t#\t#\tother";
						$caseOther ++ ;
						$casePolyplother ++ ;
						$caseDiplother ++ ;
					}
				}
				else {
					if (($GT_G1 =~ /^0.0$/)||($GT_G1 =~ /^$/)){			# 3 ou 4
	
						###############################################################################################
						# PHASING
						###############################################################################################
						$phasedornot = 0 ;
						my $is_3 = 0;
						my $is_4 = 0;
						if ($GT_poly =~/\|/){
							
							#print STDOUT $s."\n";
							my $ref = $phased_regions{$s};
							my %hash = %$ref ;
							foreach my $num_reg(sort(keys(%hash))){
								if  (exists $phased_regions{$s}{$num_reg}{$c}) {
									$genotype = $phased_regions{$s}{$num_reg}{$c} ;
									my $ref2 = $phased_regions{$s}{$num_reg};
									my %hash2 = %$ref2 ;
									foreach my $pos(sort(keys(%hash2))){
										if (exists $five{$s}{$pos}){
											if (($five{$s}{$pos} =~ /0.1/ && $GT_poly =~ /0.1/) or ($five{$s}{$pos} =~ /1.0/ && $GT_poly =~ /1.0/)){
												print HTMLSNP "<td class=\"ted2\">".$ratio.":".$ratio2."<br>";                                                  # RATIO %
					                                                        print HTMLSNP "<img src=\"".$REPimages."r1.png\" height=10 width=".$ratio.">";                  # IMG RATIO 1
                                        					                print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=10 width=".$ratio2."></td>";    # IMG RATIO 2
					                                                        print HTMLSNP "<td class=\"ted2\">".$FDP."</td>";
                                        					                print HTMLSNP "<td class=\"ted2\">".$DP_P."</td>";
												print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."3.png\" height=".$tailleImg."></td>";
												$case3 ++ ;
												$is_3 = 1;
											}
											else {
												print HTMLSNP "<td class=\"ted2\">".$ratio.":".$ratio2."<br>";                                                  # RATIO %
                                                                                                print HTMLSNP "<img src=\"".$REPimages."r1.png\" height=10 width=".$ratio.">";                  # IMG RATIO 1
                                                                                                print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=10 width=".$ratio2."></td>";    # IMG RATIO 2
                                                                                                print HTMLSNP "<td class=\"ted2\">".$FDP."</td>";
                                                                                                print HTMLSNP "<td class=\"ted2\">".$DP_P."</td>";

												print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."4.png\" height=".$tailleImg."></td>";
												$case4 ++ ;
												$is_4 = 1;
											}
											$phasedornot = 1 ;
											last ;
										}
									}
								}
							}
							#print STDOUT "\n";		
						}
						if ($phasedornot == 0) {
							print HTMLSNP "<td class=\"ted2\">".$ratio.":".$ratio2."<br>";                                                  # RATIO %
                                                        print HTMLSNP "<img src=\"".$REPimages."r1.png\" height=10 width=".$ratio.">";                  # IMG RATIO 1
                                                        print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=10 width=".$ratio2."></td>";    # IMG RATIO 2
                                                        print HTMLSNP "<td class=\"ted2\">".$FDP."</td>";
                                                        print HTMLSNP "<td class=\"ted2\">".$DP_P."</td>";
							print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."3ou4.png\" height=".$tailleImg."></td>";
							$case3ou4 ++ ;
						}
						
						###############################################################################################
						
	
	
						print TABSNP "not OK\t".$ratio."\t".$ratio2."\t".$FDP."\t".$DP_P."\t";	
						if ($is_3){print TABSNP "3";}
						elsif ($is_4){print TABSNP "4";}
						else{print TABSNP "3or4";}
						
					}
					else { #0/1
										# heterozygosity G1
						print HTMLSNP "<td>#</td><td>#</td><td>#</td>";
						print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."other.png\" height=".$tailleImg."><img src=\"".$REPimages."HG1.png\" height=".$tailleImg."></td>";
						print TABSNP "not OK\t#\t#\t#\t#\tother,heterozygosity for genome 1";
						$nbGenome2 ++ ;
						$casePolyplother ++ ;
						$caseOther ++ ;
					}
				}
			}
				
			if (($GT_poly =~ /^1.1$/)) { #	POLYPLOID NE REFERENCE - [1/1]
				if ($GT_G1 && ($GT_G1 !~ /^0.0$/) && ($GT_G1 !~ /^1.1$/)){ # SNP Genome1 intra [0/1] [0|1] [1|0]
					print HTMLSNP "<td>#</td><td>#</td><td>#</td>";
					print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."other.png\" height=".$tailleImg."><img src=\"".$REPimages."HG1.png\" height=".$tailleImg."></td>";	
					print TABSNP "not OK\t#\t#\t#\t#\tother,heterozygosity for genome 1";
					$nbGenome2 ++ ;
					$caseOther ++ ;
				}
				elsif (!$GT_G1){ # POLYPLOID A/A DIPLOIDS T/T
					print HTMLSNP "<td>#</td><td>#</td><td>#</td>";
					print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."other.png\" height=".$tailleImg."></td>";	
					print TABSNP "not OK\t#\t#\t#\t#\tother";
					$caseOther ++ ;
				}
				
				else { # Pas de SNP Genome1 [0/0] [0|0] [1/1] [1|1] SNP entre genome1 et genome2
					if ($GT_G1 =~ /^0.0$/){									# Other [NOTHING]
						print HTMLSNP "<td>#</td><td>#</td><td>#</td>";
						print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."other.png\" height=".$tailleImg."></td>";
						print TABSNP "not OK\t#\t#\t#\t#\tother";
						$caseOther ++ ;
					}
					else {
						if ($code_G1 eq $alAltP) {							# 2
							print HTMLSNP "<td>#</td><td>#</td><td>#</td>";
							print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."2.png\" height=".$tailleImg."></td>";
							print TABSNP "not OK\t#\t#\t#\t#\t2";
							$case2 ++ ;
						}
						else {											# Other [SNP DIPLO]
							print HTMLSNP "<td>#</td><td>#</td><td>#</td>";
							print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."other.png\" height=".$tailleImg."></td>";
							print TABSNP "not OK\t#\t#\t#\t#\tother";
							$caseOther ++ ;
							$caseDiplother ++ ;
						}
					}
				}
			}
			
			if (($GT_poly =~ /^0.0$/)||($GT_poly =~ /^$/)) { #POLYPLOID == REFERENCE - [0|0]
				####################################
				if (($GT_G1 !~ /^0.0$/) && ($GT_G1 !~ /^1.1$/)){ # SNP Genome1 intra [0/1] [0|1] [1|0]
					print HTMLSNP "<td>#</td><td>#</td><td>#</td>";
					print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."other.png\" height=".$tailleImg."><img src=\"".$REPimages."HG1.png\" height=".$tailleImg."></td>";
					print TABSNP "not OK\t#\t#\t#\t#\tother,heterozygosity for genome 1";
					$nbGenome2 ++ ;
					$caseOther++;
				}
				else { # Pas de SNP Genome1 [0/0] [0|0] [1/1] [1|1] SNP entre genome1 et genome2
					if ($GT_G1 =~ /^1.1$/){
						print HTMLSNP "<td>#</td><td>#</td><td>#</td>";
						print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."1.png\" height=".$tailleImg."></td>";
						print TABSNP "not OK\t#\t#\t#\t#\t1";
						$case1 ++ ;
					}	
				}	
			}
		
			print TABSNP "\n";
			print HTMLSNP "</tr>\n";
			$ligneOK = 1 ;
		}
	
		if ($ligneInter == 0) {	
			print HTMLCOUNT "<td class=\"ted2\" style=\"border-right:3px solid black\">".$s."</td>";
		}
		else {
			print HTMLCOUNT "<td class=\"ted\" style=\"border-right:3px solid black\">".$s."</td>";
		}
		print TABCOUNT $s."\t";	
			
		#######################################
		if ($ligneInter == 0) { $ligneInter = 1 ; }
		else { $ligneInter = 0 ; }
		#######################################
		
		# Calcul des intervalles #
		##########################
		$taille_totale = 0 ;
		my $ref = $intervalle2{$s};
		my %hash = %$ref;
		
		foreach my $interval(keys(%hash)){
			my @pos = split(/-/,$interval);
			$taille_inter = $pos[1]-$pos[0]+1 ;
			$taille_totale = $taille_totale + $taille_inter;
		}
		$total1 = $case5 + $case1 + $case2 + $caseDiplother + $nbGenome2;
		$total2 = $case5 + $case3ou4 + $casePolyplother;
		

		
		# SYNTHESIS
		print HTMLCOUNT "<td>".$taille_totale."</td><td style=\"border-left:3px solid black\">".$taille."</td></td>";	
		print HTMLCOUNT "<td style=\"border-left:3px solid black\">";
		print HTMLCOUNT $case1."</td><td>".$case2."</td><td>".$case3ou4."</td><td>";
		print HTMLCOUNT $case3."</td><td>".$case4."</td><td>".$case5."</td><td>";
		print HTMLCOUNT $caseOther."</td><td style=\"border-left:3px solid black\">".$nbGenome2."</td><td style=\"border-left:3px solid black\">".$total1."</td>";
		print HTMLCOUNT "<td style=\"border-left:3px solid black\">".$total2."</td>";
		print TABCOUNT $taille_totale."\t".$taille."\t";
		print TABCOUNT $case1."\t".$case2."\t".$case3ou4."\t".$case3."\t".$case4."\t".$case5."\t".$caseOther."\t".$nbGenome2."\t".$total1."\t".$total2."\t";

		$nbTotGenesAna ++ ;
	
		if ($case5 != 0) {
			$nbTotGenesVal ++ ;
			if ($ligneInter == 1) {	
				print HTMLCOUNT "<td class=\"ted2\" style=\"border-left:3px solid black\">";
				print HTMLCOUNT sprintf("%.0f", $moyenneSNPindep1/$case5).":".sprintf("%.0f", $moyenneSNPindep2/$case5)."<br>";
				print HTMLCOUNT "<img src=\"".$REPimages."r1.png\" height=10 width=".sprintf("%.0f", $moyenneSNPindep1/$case5).">";
				print HTMLCOUNT "<img src=\"".$REPimages."r2.png\" height=10 width=".sprintf("%.0f", $moyenneSNPindep2/$case5)."></td>";
			}
			else {
				print HTMLCOUNT "<td class=\"ted\" style=\"border-left:3px solid black\">";
				print HTMLCOUNT sprintf("%.0f", $moyenneSNPindep1/$case5).":".sprintf("%.0f", $moyenneSNPindep2/$case5)."<br>";
				print HTMLCOUNT "<img src=\"".$REPimages."r1.png\" height=10 width=".sprintf("%.0f", $moyenneSNPindep1/$case5).">";
				print HTMLCOUNT "<img src=\"".$REPimages."r2.png\" height=10 width=".sprintf("%.0f", $moyenneSNPindep2/$case5)."></td>";
			}
			print TABCOUNT sprintf("%.0f", $moyenneSNPindep1/$case5)."\t".sprintf("%.0f", $moyenneSNPindep2/$case5)."\t";	
		}
		else {
			print HTMLCOUNT "<td style=\"border-left:3px solid black\"></td>";
			print TABCOUNT "\t";
		}
		print HTMLCOUNT "</tr>";
		print TABCOUNT "\n";
		
		$totalSize = $totalSize + $taille_totale ;
		$totalSNP = $totalSNP + $taille ;
		$total11 = $total11 + $case1 ;
		$total22 = $total22 + $case2 ;
		$total3ou4 = $total3ou4 + $case3ou4 ;
		$total3 = $total3 + $case3 ;
		$total4 = $total4 + $case4 ;
		$total5 = $total5 + $case5 ;
		$total512 = $total512 + $total1 ;
		$total534 = $total534 + $total2 ;
		$totalOther = $totalOther + $caseOther ;
		$totalGenome2 = $totalGenome2 + $nbGenome2 ;
	}
	########## MODIF DERNIERE MINUTE ################"
	print HTMLCOUNT "<tr class=\"td3\">\n<td>";

	print HTMLCOUNT $nbTotGenesAna."<td style=\"border-left:3px solid black\">";
	print HTMLCOUNT $totalSize."</td><td style=\"border-left:3px solid black\">";
	print HTMLCOUNT $totalSNP."</td><td style=\"border-left:3px solid black\">";
	print HTMLCOUNT $total11."</td><td>";
	print HTMLCOUNT $total22."</td><td>";
	print HTMLCOUNT $total3ou4."</td><td>";
	print HTMLCOUNT $total3."</td><td>";
	print HTMLCOUNT $total4."</td><td>";
	print HTMLCOUNT $total5."</td><td>";
	print HTMLCOUNT $totalOther."</td><td style=\"border-left:3px solid black\">";
	print HTMLCOUNT $totalGenome2."</td><td style=\"border-left:3px solid black\">";
	print HTMLCOUNT $total512."</td><td style=\"border:3px solid black\">";
	print HTMLCOUNT $total534."</td><td style=\"border:3px solid black\">";
	
	print HTMLCOUNT $nbTotGenesVal."</td>";

	print HTMLCOUNT "</tr>";


	print TABCOUNT "$nbTotGenesAna\t$totalSize\t$totalSNP\t$total11\t$total22\t$total3ou4\t$total3\t$total4\t$total5\t$totalOther\t$totalGenome2\t$total512\t$total534\t";
	print TABCOUNT "\n";

	####################################################				
	print HTMLSNP "</table>\n";
	print HTMLSNP "</html>\n";
	close HTMLSNP ;

	print HTMLCOUNT "</table>\n";
	print HTMLCOUNT "</html>\n";
	close HTMLCOUNT ;

	close TABSNP;
	close TABCOUNT ;

		tie @array, 'Tie::File', $SNP_count or die ;
	
		$array[113] = "<tr><td class=\"ted3\">".$nbTotGenesAna."</td><td class=\"ted3\" style=\"border-left:3px solid black\">".$totalSize."</td><td class=\"ted3\" style=\"border-left:3px solid black\">".$totalSNP."</td></td>";	
		$array[114] = "<td class=\"ted3\" style=\"border-left:3px solid black\">";
		$array[115] = $total11."</td><td class=\"ted3\">".$total22."</td><td class=\"ted3\">".$total3ou4."</td><td class=\"ted3\">";
		$array[116] = $total3."</td><td class=\"ted3\">".$total4."</td><td class=\"ted3\">".$total5."</td><td class=\"ted3\">";
		$array[117] = $totalOther."</td><td class=\"ted3\" style=\"border-left:3px solid black\">".$totalGenome2."</td><td class=\"ted3\" style=\"border-left:3px solid black\">".$total512."</td>";
		$array[118] = "<td class=\"ted3\" style=\"border-left:3px solid black\">".$total534."</td><td class=\"ted3\">".$nbTotGenesVal."</td></tr>";		
}




sub ext_output { # Extern reference

print TABCOUNT "Gene;Interval Size Analysed (pb);nb SNP;1;2;3 or 4;5;other;SNP Diploids;SNP Polyploid;Ratio (%) $genome2Name:$genome1Name\n";

foreach my $s(sort(keys(%snp))){
	
	#######################################
	if ($ligneInter == 0) {
		print HTMLSNP "<tr class=\"bord1\">\n";
		print HTMLCOUNT "<tr class=\"td1\ border-width =\"3px\">\n";
	}
	else {
		print HTMLSNP "<tr class=\"bord2\">\n";
		print HTMLCOUNT "<tr class=\"td2\">\n";
	}
	#######################################
	

	my $ref = $snp{$s};
	my %hash = %$ref;
	$taille = keys(%hash);
		# taille de la ligne gene
		print HTMLSNP "<td  rowspan=\"$taille\"><b>$s</b></td>";
		
	$ligneOK = 0 ;
	
	$nbPolyploid = 0;
	$nbGenomes = 0 ;
	$nbCommuns = 0;
	$nbPoly_only = 0 ;
	$nbSub_only = 0 ;
	
	$case5 = 0 ;
	$case1 = 0 ;
	$case2 = 0 ;
	$case3ou4 = 0;
	$caseOther = 0;
	$caseGenome2 = 0 ;
	
	#Moyenne Ponderee
	$moyenneSNPindep1 = 0 ;
	$moyenneSNPindep2 = 0 ;
	
	@tabTrie = sort ({ $a <=> $b }keys %hash);
	
	#foreach my $c(sort ({$hash{$a} <=> $hash{$b}} keys %hash)) {
	foreach my $c(@tabTrie) {
		$nb1 = 0 ;
		$nb2 = 0 ;
		$nb3 = 0 ;
		if ($ligneOK == 1) {
			if ($ligneInter == 0) {
				print HTMLSNP "<tr class=\"td1\">\n";
			}
			else {
				print HTMLSNP "<tr class=\"td2\">\n";
			}
		}
		### Recuperation des informations ###
		($code_snp,$ratio,$GT_poly,$GT_G1,$GT_G2,$DP_P) = split(/\t/,$snp{$s}{$c});
		@recupAlleles = split(/\[/,$code_snp);
		@recupAlleles = split(/\]/,$recupAlleles[1]);
		($alRef,$alAlt) = split(/\//,$recupAlleles[0]);
		$noSNPpoly = "ok" ;
		
		#print STDOUT "\n $c $GT_poly $GT_G1 $GT_G2";
		print HTMLSNP "<td>$c</td>";
		print TABSNP "$c;";
			#####################################################################
			#	SNP POLYPLOID - [0/1] [0|1] [1|0]
			#####################################################################
			if (($GT_poly =~ /^0.1$/) || ($GT_poly =~ /^1.0$/) ) {		# Polyploid [0.1] [1.0]
				$ratio = sprintf("%.0f", $ratio);
				$ratio2 = 100-$ratio ;
				$moyenneSNPindep1 = $moyenneSNPindep1 + $ratio ;
				$moyenneSNPindep2 = $moyenneSNPindep2 + $ratio2 ;
				print HTMLSNP "<td>$code_snp</td>";
				print TABSNP "$code_snp;";
				if (($GT_G1 =~ /^1.1$/)){ #  Genome1  [1/1] [1|1]
					print HTMLSNP "<td>$alAlt</td>";
					print TABSNP "$alAlt;";
					if($GT_G2 =~ /^1.1$/){ # Genome2 = Alt [1.1]
						print HTMLSNP "<td>$alAlt</td>";
						print HTMLSNP "<td></td>";
						print HTMLSNP "<td></td>";
						print HTMLSNP "<td></td>";
						print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."3ou4.png\" height=30></td>";
						print TABSNP "$alAlt;";
						print TABSNP ";";
						print TABSNP ";";
						print TABSNP ";";
						print TABSNP "3 or 4;";
					}
					else {
						if($GT_G2 =~ /^0.0$/){ # Genome2 = Ref [0.0]
							print HTMLSNP "<td>$alRef</td>";
							print HTMLSNP "<td>OK</td>";
							print HTMLSNP "<td>[$ratio/$ratio2]</td>";
							print HTMLSNP "<td></td>";
							print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5v.png\" height=30></td>";
							print TABSNP "$alRef";
							print TABSNP "OK";
							print TABSNP "[$ratio/$ratio2;]";
							print TABSNP "$DP_P;";
							print TABSNP "<td class=\"final\"><img src=\"".$REPimages."5v.png\" height=30></td>";
						}
						else { # Genome2 = SNP [0.1]
							print HTMLSNP "<td>$code_snp</td>";
							print HTMLSNP "<td>OK</td>";
							print HTMLSNP "<td>[$ratio/$ratio2]</td>";
							print HTMLSNP "<td></td>";
							print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5prime.png\" height=30></td>";		
							print TABSNP "$code_snp;";
							print TABSNP "OK;";
							print TABSNP "[$ratio/$ratio2];";
							print TABSNP "$DP_P;";
							print TABSNP "5';";		
						}
					}
				}
				else{
					if ($GT_G1 =~ /^0.0$/){  # # Genome1  [0/0] [0|0]
						print HTMLSNP "<td>$alRef</td>";
						print TABSNP "$alRef;";
						if($GT_G2 =~ /^1.1$/){ # Genome2 = Alt [1.1]
							print HTMLSNP "<td>$alAlt</td>";
							print HTMLSNP "<td>OK</td>";
							print HTMLSNP "<td>[$ratio/$ratio2]</td>";
							print HTMLSNP "<td></td>";
							print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5v.png\" height=30></td>";
							print TABSNP "$alAlt;";
							print TABSNP "OK;";
							print TABSNP "[$ratio/$ratio2];";
							print TABSNP ";";
							print TABSNP "5;";
						}
						else {
							if($GT_G2 =~ /^0.0$/){ # Genome2 = Ref [0.0]
								print HTMLSNP "<td>$alRef</td>";
								print HTMLSNP "<td></td>";
								print HTMLSNP "<td></td>";
								print HTMLSNP "<td></td>";
								print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."3ou4.png\" height=30></td>";
								print TABSNP "$alRef;";
								print TABSNP ";";
								print TABSNP ";";
								print TABSNP ";";
								print TABSNP "3 or 4;";
							}
							else { # Genome2 = SNP [0.1]
								print HTMLSNP "<td>$code_snp</td>";
								print HTMLSNP "<td>OK</td>";
								print HTMLSNP "<td>[$ratio/$ratio2]</td>";
								print HTMLSNP "<td></td>";
								print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5prime.png\" height=30></td>";
								print TABSNP "$code_snp;";
								print TABSNP "OK;";
								print TABSNP "[$ratio/$ratio2];";
								print TABSNP ";";
								print TABSNP "5';";		
								}
							}
						}		
						else { # SNP Genome1  [0/1] [0|1] [1|0]
							print HTMLSNP "<td>$code_snp</td>"; #REF
							print HTMLSNP "$code_snp;"; #REF
							if($GT_G2 =~ /^1.1$/){ # Genome2 = Alt [1.1]
								print HTMLSNP "<td>$alAlt</td>";
								print HTMLSNP "<td>OK</td>";
								print HTMLSNP "<td>[$ratio/$ratio2]</td>";
								print HTMLSNP "<td></td>";
								print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5prime.png\" height=30></td>";
								print TABSNP "$alAlt;";
								print TABSNP "OK;";
								print TABSNP "[$ratio/$ratio2];";
								print TABSNP ";";
								print TABSNP "5'";
							}
							else {
								if($GT_G2 =~ /^0.0$/){ # Genome2 = Ref [0.0]
									print HTMLSNP "<td>$alRef</td>";
									print HTMLSNP "<td>OK</td>";
									print HTMLSNP "<td>[$ratio/$ratio2]</td>";
									print HTMLSNP "<td></td>";
									print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5prime.png\" height=30></td>";
									print TABSNP "$alRef;";
									print TABSNP "OK;";
									print TABSNP "[$ratio/$ratio2];";
									print TABSNP ";";
									print TABSNP "5';";
								}
								else { # Genome2 = SNP [0.1]
									print HTMLSNP "<td>$code_snp</td>";
									print HTMLSNP "<td></td>";
									print HTMLSNP "<td></td>";
									print HTMLSNP "<td></td>";
									print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5second.png\" height=30></td>";	
									print TABSNP "$code_snp;";
									print TABSNP ";";
									print TABSNP ";";
									print TABSNP ";";
									print TABSNP "5''";		
								}
							}
						}
						print TABSNP "\n";
						$ligneOK = 1 ;
						print HTMLSNP "</tr>\n";
					}
			}
			#####################################################################
			#	POLYPLOID NE REFERENCE - [1/1]
			#####################################################################
			#print STDOUT "GTPOLY:$GT_poly@";
			if (($GT_poly =~ /^1.1$/) ) {	
				print HTMLSNP "<td>$alAlt</td>";
				print TABSNP "$alAlt;";
				####################################
				if (($GT_G1 !~ /^0.0$/) && ($GT_G1 !~ /^1.1$/)){ # Genome 1 [0/1] [0|1] [1|0]
					print HTMLSNP "<td>$code_snp</td>";
					print TABSNP "$code_snp";
					if ($GT_G2 =~ /^1.1$/) {
						print HTMLSNP "<td>$alAlt</td>";
						print HTMLSNP "<td></td>";
						print HTMLSNP "<td></td>";
						print HTMLSNP "<td></td>";
						print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5primeou1.png\" height=30></td>";
						print TABSNP "$alAlt;;;;5' or 1;";
					}
					else {
						if  ($GT_G2 =~ /^0.0$/) {
							print HTMLSNP "<td>$alRef</td>";
							print HTMLSNP "<td></td>";
							print HTMLSNP "<td></td>";
							print HTMLSNP "<td></td>";
							print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5primeou1.png\" height=30></td>";
							print TABSNP "$altRef;;;;5' or 1;";
						}
						else { # Genome2 = SNP [0.1]
							print HTMLSNP "<td>$code_snp</td>";
							print HTMLSNP "<td></td>";
							print HTMLSNP "<td></td>";
							print HTMLSNP "<td></td>";
							print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5second.png\" height=30></td>";		
							print TABSNP "$code_snp;;;;5'';";
						}
					}
				}
				else { 
					if ($GT_G1 =~ /^0.0$/){ #  Genome1 [0/0] [0|0]
					###################
						print HTMLSNP "<td>$alRef</td>";
						print TABSNP "$alRef";
						if ($GT_G2 =~ /^1.1$/) { # Genome2 = Alt [1.1]	
							print HTMLSNP "<td>$alAlt</td>";
							print HTMLSNP "<td></td>";
							print HTMLSNP "<td></td>";
							print HTMLSNP "<td></td>";
							print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."1.png\" height=30></td>";
							print TABSNP "$alAlt;;;;1;";
						}
						else {
							if  ($GT_G2 =~ /^0.0$/) {
								print HTMLSNP "<td>$alAlt</td>";
								print HTMLSNP "<td></td>"; 
								print HTMLSNP "<td></td>";
								print HTMLSNP "<td></td>";
								print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."3ou4.png\" height=30></td>";
								print TABSNP "$alAlt;;;;3 or 4;";
							}
							else { # Genome2 = SNP [0.1]
								print HTMLSNP "<td>$code_snp</td>";
								print HTMLSNP "<td></td>";
								print HTMLSNP "<td></td>";
								print HTMLSNP "<td></td>";
								print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5primeou2.png\" height=30></td>";
								print TABSNP "$code_snp;;;;5' or 2";		
							}	
						}
					}
					else { #  [1/1] [1|1] Genome 1
						print HTMLSNP "<td>$alAlt</td>";
						print TABSNP "$alAlt";
						if ($GT_G2 =~ /^1.1$/)  { # Genome2 = Alt [1.1]	
							print HTMLSNP "<td>$alAlt</td>";
							print HTMLSNP "<td>OK</td>";
							print HTMLSNP "<td>[$ratio/$ratio2]</td>";
							print HTMLSNP "<td>$DP_P</td>";
							print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5v.png\" height=30></td>";
							print TABSNP "$alAlt;OK;[$ratio/$ratio2];$DP_P;5;";
						}
						else {
							if  ($GT_G2 =~ /^0.0$/) {
								print HTMLSNP "<td>$alRef</td>";
								print HTMLSNP "<td></td>"; 
								print HTMLSNP "<td></td>";
								print HTMLSNP "<td></td>";
								print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."2.png\" height=30></td>";
								print TABSNP "$alRef;;;;2;";
							}
							else { # Genome2 = SNP [0.1]
								print HTMLSNP "<td>$code_snp</td>";
								print HTMLSNP "<td></td>";
								print HTMLSNP "<td></td>";
								print HTMLSNP "<td></td>";
								print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5primeou2.png\" height=30></td>";
								print TABSNP "$code_snp;;;;5' or 2" ;		
							}	
						}
					}
				}
				$nbPolyploid ++ ;
				print TABSNP "\n";
				$ligneOK = 1 ;
				print HTMLSNP "</tr>\n";
			
		}	
		#####################################################################
		#	POLYPLOID == REFERENCE - [0|0]
		#####################################################################
		if (($GT_poly =~ /^0.0$/) ) {	
			print HTMLSNP "<td>$alRef</td>";
			print TABSNP "$alRef;";
			if (($GT_G1 !~ /^0.1$/) && ($GT_G1 !~ /^1.0$/)){				 # Genome1 intra [0/1] [0|1] [1|0]
				print HTMLSNP "<td>$code_snp</td>";
				print TABSNP "$code_snp";
				if ($GT_G2 =~ /^1.1$/)  { 													# Genome2 = Alt [1.1]	
					print HTMLSNP "<td>$alAlt</td>";
					print HTMLSNP "<td></td>";
					print HTMLSNP "<td></td>";
					print HTMLSNP "<td></td>";
					print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5primeou1.png\" height=30></td>";
					print TABSNP "$alAlt;;;;5' or 1";
				}
				else {
					if  ($GT_G2 =~ /^0.0$/) { 													# Genome2 = Ref [0.0]
						print HTMLSNP "<td>$alRef</td>";
						print HTMLSNP "<td></td>";
						print HTMLSNP "<td></td>";
						print HTMLSNP "<td></td>";
						print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5primeou1.png\" height=30></td>";
						print TABSNP "$alRef;;;;5' or 1";
					}
					else { # Genome2 = SNP [0.1]
						if ( ($GT_G2 =~ /^0.1$/) || ($GT_G2 =~ /^1.0$/)){ 													# Genome2 = SNP [0.1]
							print HTMLSNP "<td>$code_snp</td>";
							print HTMLSNP "<td></td>";
							print HTMLSNP "<td></td>";
							print HTMLSNP "<td></td>";
							print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5second.png\" height=30></td>";	
							print TABSNP "$code_snp;;;;5''";
						}	
					}	
				}
			}
			
			else { 
				if ($GT_G1 =~ /^0.0$/){ # Pas de SNP Genome1 [0/0] [0|0]
					print HTMLSNP "<td>$alRef</td>";
					print TABSNP "$alRef";
					if ($GT_G2 =~ /^1.1$/) { # Genome2 = Alt [1.1]	
						print HTMLSNP "<td>$alAlt</td>";
						print HTMLSNP "<td></td>";
						print HTMLSNP "<td></td>";
						print HTMLSNP "<td></td>";
						print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."2.png\" height=30></td>";
						print TABSNP "$alAlt;;;;2;";
					}
					else {
						if ( ($GT_G2 =~ /^0.1$/) || ($GT_G2 =~ /^1.0$/)){
							print HTMLSNP "<td>$code_snp</td>";
							print HTMLSNP "<td>OK</td>";
							print HTMLSNP "<td>[$ratio/$ratio2]</td>";
							print HTMLSNP "<td>$DP_P</td>";
							print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5primeou2.png\" height=30></td>";
							print TABSNP "$code_snp;OK;[$ratio/$ratio2];$DP_P;5' or 2";
						}	
					}
				}
					
					else { #  [1/1] [1|1] SNP entre genome1 et genome2
						print HTMLSNP "<td>$alAlt</td>";
						print TABSNP "$alAlt";
						if ($GT_G2 =~ /^1.1$/)  { # Genome2 = Alt [1.1]	
							print HTMLSNP "<td>$alAlt</td>";
							print HTMLSNP "<td></td>";
							print HTMLSNP "<td></td>";
							print HTMLSNP "<td></td>";
							print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."other.png\" height=30></td>";
							print TABSNP "$alAlt;;;;other;";
						}
						else {
							if  ($GT_G2 =~ /^0.0$/) {
								print HTMLSNP "<td>$alRef</td>";
								print HTMLSNP "<td></td>";
								print HTMLSNP "<td></td>";
								print HTMLSNP "<td></td>";
								print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."1.png\" height=30></td>";
								print TABSNP "$alRef;;;;1;";
							}
							else { # Genome2 = SNP [0.1]
								print HTMLSNP "<td>$code_snp</td>";
								print HTMLSNP "<td></td>";
								print HTMLSNP "<td></td>";
								print HTMLSNP "<td></td>";
								print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5primeou2.png\" height=30></td>";
								print TABSNP "$code_snp;;;;5' or 2;";
							}	
						}
					}
						
				}		
			}
			$nbSub_only ++;
			print TABSNP "\n";
			$ligneOK = 1 ;
			print HTMLSNP "</tr>\n";
		}
	
	if ($ligneInter == 0) {
		print HTMLCOUNT "<td class=\"ted2\" style=\"border-right:3px solid black\">$s</td>";
	}
	else {
		print HTMLCOUNT "<td class=\"ted\" style=\"border-right:3px solid black\">$s</td>";
	}
	
	print TABCOUNT "Gene;Interval Size;Analysed (pb);nb SNP;1;2;3 or 4;5;5';5'';5' or 1;5'' or 2;other;SNP Diploids;SNP Polyploid;Ratio (%) $genome2Name:$genome1Name;";
	
	#######################################
	if ($ligneInter == 0) {
		$ligneInter = 1 ;
	}
	else {
		$ligneInter = 0 ;
	}
	#######################################
	
	# Calcul des intervalles
	$taille_totale = 0 ;
	my $ref = $intervalle2{$s};
	my %hash = %$ref;
		
		foreach my $interval(keys(%hash)){
			my @pos = split(/-/,$interval);
			$taille_inter = $pos[1]-$pos[0]-1 ;
			$taille_totale = $taille_totale + $taille_inter;
		}
	
	$total1 = $case5 + $case1 + $case2 ;
	$total2 = $case5 + $case3ou4 ;
	
	print HTMLCOUNT "<td>$taille_totale</td><td style=\"border-left:3px solid black\">$taille</td></td>";	
	print HTMLCOUNT "<td style=\"border-left:3px solid black\">$case1</td><td>$case2</td><td>$case3ou4</td><td>$case5</td><td>$caseGenome2</td><td>$caseOther</td><td style=\"border-left:3px solid black\">$total1</td><td style=\"border-left:3px solid black\">$total2</td>";
	if ($case5 != 0) {
		print HTMLCOUNT "<td style=\"border-left:3px solid black\">".sprintf("%.0f", $moyenneSNPindep1/$case5).":".sprintf("%.0f", $moyenneSNPindep2/$case5)."</td>";	
	}
	else {
		print HTMLCOUNT "<td style=\"border-left:3px solid black\"></td>";
	}
	print HTMLCOUNT "</tr>";
	
	$totalSize = $totalSize + $taille_totale ;
	$totalSNP = $totalSNP + $taille ;
	$total11 = $total11 + $case1 ;
	$total22 = $total22 + $case2 ;
	$total3ou4 = $total3ou4 + $case3ou4 ;
	$total5 = $total5 + $case5 ;
	$total512 = $total512 + $total1 ;
	$total534 = $total534 + $total2 ;
	$totalOther = $totalOther + $caseOther ;
	$totalGenome2 = $totalGenome2 + $caseGenome2 ;
}
########## MODIF DERNIERE MINUTE ################"
print HTMLCOUNT "<tr class=\"td3\">\n";
print HTMLCOUNT "<td></td><td style=\"border-left:3px solid black\">$totalSize</td><td style=\"border-left:3px solid black\">$totalSNP</td><td style=\"border-left:3px solid black\">$total11</td><td>$total22</td><td>$total3ou4</td><td>$total5</td><td>$totalGenome2</td><td>$totalOther</td><td style=\"border-left:3px solid black\">$total512</td><td style=\"border-left:3px solid black\">$total534</td><td style=\"border-left:3px solid black\"></td>";
print HTMLCOUNT "</tr>";
####################################################
print HTMLSNP "</table>\n";
print HTMLSNP "</html>\n";
print HTMLCOUNT "</table>\n";
print HTMLCOUNT "</html>\n";

close TABSNP;
close HTMLSNP ;
close HTMLCOUNT ;
	}

sub poly_poly_output {
	foreach my $s(sort(keys(%snp_final))){
		#######################################
		if ($ligneInter == 0) {
			print HTMLSNP "<tr class=\"bord1\">\n";
			print HTMLCOUNT "<tr class=\"td1\ border-width =\"3px\">\n";
		}
		else {
			print HTMLSNP "<tr class=\"bord2\">\n";
			print HTMLCOUNT "<tr class=\"td2\">\n";
		}
		#######################################
	
		my $ref = $snp_final{$s};
		my %hash = %$ref;
		$taille = keys(%hash);
		
		# taille de la ligne gene
		print HTMLSNP "<td  rowspan=\"".$taille."\"><b>".$s."</b></td>";
		
		$ligneOK = 0 ;
		
		$nbPolyploid1 = 0 ;		# SNP heterozygosity for P1
		$nbPolyploid2 = 0 ;		# SNP heterozygosity for P2
		$nbCommuns = 0 ; 		# SNP heterozygosity [P1] = [P2]
		$nbCommunHomo = 0 ;		# SNP homozygosity [P1] = [P2]
		$nbDifferent = 0 ;		# [P1] ne [P2]
		$alleleCommun = 0 ;		# Example : P1 = [A/G] ; P2 = [A]
		$alleleDifferent = 0 ;	# Example : P1 = [A/G] ; P2 = [C] or [T]
		$alleleCommunH = 0 ;	# Example : P1 = [A/G] ; P2 = [A/C]
		$nbHomoDiff = 0 ;
		
		@tabTrie = sort ({ $a <=> $b }keys %hash);
		
		$taille = 0;
		
		

		#foreach my $c(sort ({$hash{$a} <=> $hash{$b}} keys %hash)) {
		foreach my $c(@tabTrie) {
			
			if ($ligneOK == 1) {
				if ($ligneInter == 0) {
					print HTMLSNP "<tr class=\"td1\">\n";
				}
				else {
					print HTMLSNP "<tr class=\"td2\">\n";
				}
			}
		
			#print STDOUT "\n\n\n".$snp_final{$s}{$c} ;
			($code_snp,$AD,$GT_poly,$DP_P,$code_snp2,$AD_2,$GT_poly2,$DP_P2) = split(/\t/,$snp_final{$s}{$c});
			($DP_P, $FDP) = split(/-/, $DP_P);
			($DP_P2, $FDP2) = split(/-/, $DP_P2);
			
			#print STDOUT "\nALLELES :".$AD." - ".$AD_2;
			
			($sub1_1,$sub1_2) = split(",",$AD);
			($sub2_1,$sub2_2) = split(",",$AD_2);
			
			# $sub1_1 = sprintf("%.0f", $sub1_1);
			# $sub1_2 = sprintf("%.0f", $sub1_2);
			# $sub2_1 = sprintf("%.0f", $sub2_1);
			# $sub2_2 = sprintf("%.0f", $sub2_2);
			if ($DP_P > 0) {
				$SG1 = ($sub1_1/$DP_P*100) ;
				$SG2 = ($sub1_2/$DP_P*100) ;
			}
			if ($DP_P2 > 0) {
				$SG3 = ($sub2_1/$DP_P2*100) ;
				$SG4 = ($sub2_2/$DP_P2*100) ;
			}
			

			if ($GT_poly ne "") { # Polyploide 1 = [0.0] ou [0.1] ou [1.1]
				@recupAlleles = split(/\[/,$code_snp);
				@recupAlleles = split(/\]/,$recupAlleles[1]);
				($alRef,$alAltP) = split(/\//,$recupAlleles[0]);
				# Attribution des alleles au polyploide 1 si pas de SNP
				if ($GT_poly =~ /^0.0$/) { $code_snp = $alRef ; } 
				if ($GT_poly =~ /^1.1$/) { $code_snp = $alAltP ; }
				# Attribution des alleles au polyploide 2 si pas de SNP
				if (($GT_poly2 eq "") || ($GT_poly2 =~ /^0.0$/)) { $code_snp2 = $alRef ; }
				if ($GT_poly2 =~ /^1.1$/) {
					@recupAlleles = split(/\[/,$code_snp2);
					@recupAlleles = split(/\]/,$recupAlleles[1]);
					($alRef,$alAlt2) = split(/\//,$recupAlleles[0]);
					$code_snp2 = $alAlt2 ;
				}
			}
			elsif ($GT_poly2 ne "") { # pas de SNP polyploide 1 dans le fichier 1 (fichiers non mergés) -> equivalent de [0.0]
				@recupAlleles = split(/\[/,$code_snp2);
				@recupAlleles = split(/\]/,$recupAlleles[1]);
				($alRef,$alAlt2) = split(/\//,$recupAlleles[0]);
				# Attribution des Alleles au polyploide 2
				if ($GT_poly2 =~ /^1.1$/) { $code_snp2 = $alAlt2 ; }
				if ($GT_poly2 =~ /^0.0$/) { $code_snp2 = $alRef ; }
			}
			#print STDOUT "\n($code_snp:$GT_poly) - ($code_snp2:$GT_poly2)" ;
			$noSNPpoly = "ok" ;
			#____________________________________________________________________________________________________________________________________________			
			# [1] P1 = 0/1 ; P2 = 0/1		(2 alleles)
			#print STDOUT "\n".($code_snp2 eq $code_snp);
			if ((($GT_poly =~ /^0.1$/)||($GT_poly =~ /^1.0$/)) && (($GT_poly2 =~ /^0.1$/)||($GT_poly2 =~ /^1.0$/))) {
				
				if (($SG1 > $value_filter_p1) && ($SG2 > $value_filter_p1) && ($SG3 > $value_filter_p2) && ($SG4 > $value_filter_p2) && $code_snp2 eq $code_snp) {
					
					
					if ($ligneInter == 0) {
						
						print HTMLSNP "<td class=\"tedG2\">".$c."</td>";
						print HTMLSNP "<td class=\"tedG2\">".$alRef."</td>";
						print HTMLSNP "<td class=\"ted2\">".$code_snp."</td>";
						print HTMLSNP "<td class=\"ted2\">".$code_snp2."</td>";
						
						print TABSNP $s . "\t";
						print TABSNP $c . "\t";
						print TABSNP $alRef . "\t";
						print TABSNP $code_snp . "\t";
						print TABSNP $code_snp2 . "\t";
						print TABSNP $FDP."/".$DP_P;
						
						print HTMLSNP "<td class=\"ted2\">".$FDP."/".$DP_P;
						if (($DP_P) != 0) {
							print HTMLSNP "<br><img src=\"".$REPimages."r1.png\" height=5 width=".sprintf("%.0f", $FDP/$DP_P*100).">";
							print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=5 width=".sprintf("%.0f", 100-$FDP/$DP_P*100).">";
							print HTMLSNP "<br>".$sub1_1." - ".$sub1_2."</td>";
							print TABSNP "," . $sub1_1." - ".$sub1_2."\t";
						}
						else {
							print HTMLSNP "</td>";
							print TABSNP "\t";
						}
						
						print TABSNP $FDP2."/".$DP_P2;
						print HTMLSNP "<td class=\"ted2\">".$FDP2."/".$DP_P2;
						if (($DP_P2) != 0) {
							print HTMLSNP "<br><img src=\"".$REPimages."r1.png\" height=5 width=".sprintf("%.0f", $FDP2/$DP_P2*100).">";
							print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=5 width=".sprintf("%.0f", 100-$FDP2/$DP_P2*100).">";
							print HTMLSNP "<br>".$sub2_1." - ".$sub2_2."</td>";
							print TABSNP "," . $sub2_1." - ".$sub2_2."\t";
						}
						else {
							print HTMLSNP "</td>";
							print TABSNP "\t";
						}
						print TABSNP "\n";
						print HTMLSNP "</tr>\n";
					}
					else {
						#
							print HTMLSNP "<td class=\"tedG\">".$c."</td>";
							print HTMLSNP "<td class=\"tedG\">".$alRef."</td>";
							print HTMLSNP "<td class=\"ted\">".$code_snp."</td>";
							print HTMLSNP "<td class=\"ted\">".$code_snp2."</td>";
							print HTMLSNP "<td class=\"ted\">".$FDP." / ".$DP_P;
							print TABSNP $s . "\t";
							print TABSNP $c . "\t";
							print TABSNP $alRef . "\t";
							print TABSNP $code_snp . "\t";
							print TABSNP $code_snp2 . "\t";
							print TABSNP $FDP."/".$DP_P;
							if (($DP_P) != 0) {
								print HTMLSNP "<br><img src=\"".$REPimages."r1.png\" height=5 width=".sprintf("%.0f", $FDP/$DP_P*100).">";
								print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=5 width=".sprintf("%.0f", 100-$FDP/$DP_P*100).">";
								print HTMLSNP "<br>".$sub1_1." - ".$sub1_2."</td>";
								print TABSNP "," . $sub1_1." - ".$sub1_2."\t";
							}
							else {
								print HTMLSNP "</td>";
								print TABSNP "\t";
							}
							print TABSNP $FDP2."/".$DP_P2;
							print HTMLSNP "<td class=\"ted\">".$FDP2." / ".$DP_P2;
							if (($DP_P2) != 0) {
								print HTMLSNP "<br><img src=\"".$REPimages."r1.png\" height=5 width=".sprintf("%.0f", $FDP2/$DP_P2*100).">";
								print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=5 width=".sprintf("%.0f", 100-$FDP2/$DP_P2*100).">";
								print HTMLSNP "<br>".$sub2_1." - ".$sub2_2."</td>";
								print TABSNP "," . $sub2_1." - ".$sub2_2."\t";
							}
							else {
								print HTMLSNP "</td>";
								print TABSNP "\t";
							}
							print TABSNP "\n";
							print HTMLSNP "</tr>\n";
					}
					$nbPolyploid1 ++ ;		# SNP heterozygosity for P1
					$nbPolyploid2 ++ ;		# SNP heterozygosity for P2
					$nbCommuns ++ ; 		# SNP heterozygosity [P1] = [P2]
					$taille++;
					
				}
				else {
					if (($SG1 > $value_filter_p1) && ($SG2 > $value_filter_p1) && ($SG3 > $value_filter_p2) && ($SG4 > $value_filter_p2)) {
						if ($alAlt2 ne $alAlt) { # P1 [A/G] P2 [A/C]		(3 alleles)
							

					
							print HTMLSNP "<td style=\"border-left:3px solid black\">".$c."</td>";
							print HTMLSNP "<td style=\"border-left:3px solid black\">".$alRef."</td>";
							print HTMLSNP "<td>".$code_snp."</td>";
							print HTMLSNP "<td>".$code_snp2."</td>";
							print HTMLSNP "<td>".$FDP."/".$DP_P;
							print TABSNP $s . "\t";
							print TABSNP $c . "\t";
							print TABSNP $alRef . "\t";
							print TABSNP $code_snp . "\t";
							print TABSNP $code_snp2 . "\t";
							print TABSNP $FDP."/".$DP_P;
							if (($DP_P) != 0) {
								print HTMLSNP "<br><img src=\"".$REPimages."r1.png\" height=5 width=".sprintf("%.0f", $FDP/$DP_P*100).">";
								print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=5 width=".sprintf("%.0f", 100-$FDP/$DP_P*100).">";
								print HTMLSNP "<br>".$sub1_1." - ".$sub1_2."</td>";
								print TABSNP "," . $sub1_1." - ".$sub1_2."\t";
							}
							else {
								print HTMLSNP "</td>";
								print TABSNP "\t";
							}
							print TABSNP $FDP2."/".$DP_P2;
							print HTMLSNP "<td>".$FDP2."/".$DP_P2;
							if (($DP_P2) != 0) {
								print HTMLSNP "<br><img src=\"".$REPimages."r1.png\" height=5 width=".sprintf("%.0f", $FDP2/$DP_P2*100).">";
								print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=5 width=".sprintf("%.0f", 100-$FDP2/$DP_P2*100).">";
								print HTMLSNP "<br>".$sub2_1." - ".$sub2_2."</td>";
								print TABSNP "," . $sub2_1." - ".$sub2_2."\t";
							}
							else {
								print HTMLSNP "</td>";
								print TABSNP "\t";
							}
							print TABSNP "\n";
							print HTMLSNP "</tr>\n";
							$nbDifferent ++ ;		
							$alleleCommunH ++ ;		
							$taille++;
						}
					}
				}
			}
				
			else { # ALL
				# COMMON PART
				#print STDOUT "\nBOUM : ".$SG1." + ".$SG2." + ".$SG3." + ".$SG4 ;
				#if (($SG1> $value_filter_p1) && ($SG2> $value_filter_p1) && ($SG3> $value_filter_p2) && ($SG4> $value_filter_p2)) {
					
					# print HTMLSNP "<td></td><td></td>";
						
					# [5] P1 = 1/1 ; P2 = 1/1		(1 allele)		P1 [A] P2 [A]
					if (($GT_poly =~ /^1.1$/) && ($GT_poly2 =~ /^1.1$/)) {

						
						print HTMLSNP "<td style=\"border-left:3px solid black\">".$c."</td>";
						print HTMLSNP "<td class=\"border-left:3px solid black\">".$alRef."</td>";
						print HTMLSNP "<td>".$code_snp."</td>";
						print HTMLSNP "<td>".$code_snp2."</td>";
						print HTMLSNP "<td>".$FDP."/".$DP_P;
						print TABSNP $s . "\t";
						print TABSNP $c . "\t";
						print TABSNP $alRef . "\t";
						print TABSNP $code_snp . "\t";
						print TABSNP $code_snp2 . "\t";
						print TABSNP $FDP."/".$DP_P;
						if (($DP_P) != 0) {
							print HTMLSNP "<br><img src=\"".$REPimages."r1.png\" height=5 width=".sprintf("%.0f", $FDP/$DP_P*100).">";
							print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=5 width=".sprintf("%.0f", 100-$FDP/$DP_P*100).">";
							print HTMLSNP "<br>".$sub1_1." - ".$sub1_2."</td>";
							print TABSNP "," . $sub1_1." - ".$sub1_2."\t";
						}
						else {
							print HTMLSNP "</td>";
							print TABSNP "\t";
						}
						print TABSNP $FDP2."/".$DP_P2;
						print HTMLSNP "<td>".$FDP2."/".$DP_P2;
						if (($DP_P2) != 0) {
							print HTMLSNP "<br><img src=\"".$REPimages."r1.png\" height=5 width=".sprintf("%.0f", $FDP2/$DP_P2*100).">";
							print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=5 width=".sprintf("%.0f", 100-$FDP2/$DP_P2*100).">";
							print HTMLSNP "<br>".$sub2_1." - ".$sub2_2."</td>";
							print TABSNP "," . $sub2_1." - ".$sub2_2."\t";
						}
						else {
							print HTMLSNP "</td>";
							print TABSNP "\t";
						}
						print TABSNP "\n";
						print HTMLSNP "</tr>\n";
						#*********************************************************************************************************************
						if ($code_snp2 eq $code_snp) {
							$nbCommunHomo ++ ;
							$taille++;
						}
						else { # 						(2 alleles)		P1 [A] P2 [C]	
							$nbDifferent ++ ;
							$nbHomoDiff ++ ;
							$taille++;
						}
					}
					# [2] [4] P1 = 0/1 ; P2 = 1/1			 
					if (((($GT_poly =~ /^0.1$/) || ($GT_poly =~ /^0.1$/)) && ($GT_poly2 =~ /^1.1$/))) {
						if (($SG1> $value_filter_p1) && ($SG2> $value_filter_p1)) {
							
							
							
							print HTMLSNP "<td style=\"border-left:3px solid black\">".$c."</td>";
							print HTMLSNP "<td class=\"border-left:3px solid black\">".$alRef."</td>";
							print HTMLSNP "<td>".$code_snp."</td>";
							print HTMLSNP "<td>".$code_snp2."</td>";
							print HTMLSNP "<td>".$FDP."/".$DP_P;
							print TABSNP $s . "\t";
							print TABSNP $c . "\t";
							print TABSNP $alRef . "\t";
							print TABSNP $code_snp . "\t";
							print TABSNP $code_snp2 . "\t";
							print TABSNP $FDP."/".$DP_P;
							if (($DP_P) != 0) {
								print HTMLSNP "<br><img src=\"".$REPimages."r1.png\" height=5 width=".sprintf("%.0f", $FDP/$DP_P*100).">";
								print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=5 width=".sprintf("%.0f", 100-$FDP/$DP_P*100).">";
								print HTMLSNP "<br>".$sub1_1." - ".$sub1_2."</td>";
								print TABSNP "," . $sub1_1." - ".$sub1_2."\t";
							}
							else {
								print HTMLSNP "</td>";
								print TABSNP "\t";
							}
							print TABSNP $FDP2."/".$DP_P2;
							print HTMLSNP "<td>".$FDP2."/".$DP_P2;
							if (($DP_P2) != 0) {
								print HTMLSNP "<br><img src=\"".$REPimages."r1.png\" height=5 width=".sprintf("%.0f", $FDP2/$DP_P2*100).">";
								print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=5 width=".sprintf("%.0f", 100-$FDP2/$DP_P2*100).">";
								print HTMLSNP "<br>".$sub2_1." - ".$sub2_2."</td>";
								print TABSNP "," . $sub2_1." - ".$sub2_2."\t";
							}
							else {
								print HTMLSNP "</td>";
								print TABSNP "\t";
							}
							print TABSNP "\n";
							print HTMLSNP "</tr>\n";
							if ($alAlt2 ne $alAlt) { # 	(2 alleles)	P1 [A/G] P2 [G]
								$nbDifferent ++ ;
								$alleleCommun ++ ;
								$nbPolyploid1 ++ ;
								$taille++;
							}
							else { #					(3 alleles) P1 [A/G] P2 [C]
								$nbDifferent ++ ;
								$alleleDifferent ++ ;
								$nbPolyploid1 ++ ;
								$taille++;
							}
						}
					}
					if ((($GT_poly2 =~ /^0.1$/) || ($GT_poly2 =~ /^0.1$/)) && ($GT_poly =~ /^1.1$/)) {
						if (($SG3> $value_filter_p2) && ($SG4> $value_filter_p2)) {
							print HTMLSNP "<td style=\"border-left:3px solid black\">".$c."</td>";
							print HTMLSNP "<td class=\"border-left:3px solid black\">".$alRef."</td>";
							print HTMLSNP "<td>".$code_snp."</td>";
							print HTMLSNP "<td>".$code_snp2."</td>";
							print HTMLSNP "<td>".$FDP."/".$DP_P;
							print TABSNP $s . "\t";
							print TABSNP $c . "\t";
							print TABSNP $alRef . "\t";
							print TABSNP $code_snp . "\t";
							print TABSNP $code_snp2 . "\t";
							print TABSNP $FDP."/".$DP_P;
							if (($DP_P) != 0) {
								print HTMLSNP "<br><img src=\"".$REPimages."r1.png\" height=5 width=".sprintf("%.0f", $FDP/$DP_P*100).">";
								print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=5 width=".sprintf("%.0f", 100-$FDP/$DP_P*100).">";
								print HTMLSNP "<br>".$sub1_1." - ".$sub1_2."</td>";
								print TABSNP "," . $sub1_1." - ".$sub1_2."\t";
							}
							else {
								print HTMLSNP "</td>";
								print TABSNP "\t";
							}
							print TABSNP $FDP2."/".$DP_P2;
							print HTMLSNP "<td>".$FDP2."/".$DP_P2;
							if (($DP_P2) != 0) {
								print HTMLSNP "<br><img src=\"".$REPimages."r1.png\" height=5 width=".sprintf("%.0f", $FDP2/$DP_P2*100).">";
								print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=5 width=".sprintf("%.0f", 100-$FDP2/$DP_P2*100).">";
								print HTMLSNP "<br>".$sub2_1." - ".$sub2_2."</td>";
								print TABSNP "," . $sub2_1." - ".$sub2_2."\t";
							}
							else {
								print HTMLSNP "</td>";
								print TABSNP "\t";
							}
							print TABSNP "\n";
							print HTMLSNP "</tr>\n";
							if ($alAlt2 ne $alAlt) { # 	(2 alleles)	P1 [A/G] P2 [G]
								$nbDifferent ++ ;
								$alleleCommun ++ ;
								$nbPolyploid2 ++ ;
								$taille++;
							}
							else { #					(3 alleles) P1 [A/G] P2 [C]
								$nbDifferent ++ ;
								$alleleDifferent ++ ;
								$nbPolyploid2 ++ ;
								$taille++;
							}
						}
					}
					# [3] [7] P1 = 0/1 ; P2 = 0/0 (2 alleles) P1 [A/G] P2 [A]
					if ((($GT_poly =~ /^0.1$/) || ($GT_poly =~ /^0.1$/)) && (($GT_poly2 =~ /^0.0$/) || ($GT_poly2 eq ""))) {
						if (($SG1> $value_filter_p1) && ($SG2> $value_filter_p1)) {
							print HTMLSNP "<td style=\"border-left:3px solid black\">".$c."</td>";
							print HTMLSNP "<td class=\"border-left:3px solid black\">".$alRef."</td>";
							print HTMLSNP "<td>".$code_snp."</td>";
							print HTMLSNP "<td>".$code_snp2."</td>";
							print HTMLSNP "<td>".$FDP."/".$DP_P;
							print TABSNP $s . "\t";
							print TABSNP $c . "\t";
							print TABSNP $alRef . "\t";
							print TABSNP $code_snp . "\t";
							print TABSNP $code_snp2 . "\t";
							print TABSNP $FDP."/".$DP_P;
							if (($DP_P) != 0) {
								print HTMLSNP "<br><img src=\"".$REPimages."r1.png\" height=5 width=".sprintf("%.0f", $FDP/$DP_P*100).">";
								print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=5 width=".sprintf("%.0f", 100-$FDP/$DP_P*100).">";
								print HTMLSNP "<br>".$sub1_1." - ".$sub1_2."</td>";
								print TABSNP "," . $sub1_1." - ".$sub1_2."\t";
							}
							else {
								print HTMLSNP "</td>";
								print TABSNP "\t";
							}
							print TABSNP $FDP2."/".$DP_P2;
							print HTMLSNP "<td>".$FDP2."/".$DP_P2;
							if (($DP_P2) != 0) {
								print HTMLSNP "<br><img src=\"".$REPimages."r1.png\" height=5 width=".sprintf("%.0f", $FDP2/$DP_P2*100).">";
								print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=5 width=".sprintf("%.0f", 100-$FDP2/$DP_P2*100).">";
								print HTMLSNP "<br>".$sub2_1." - ".$sub2_2."</td>";
								print TABSNP "," . $sub2_1." - ".$sub2_2."\t";
							}
							else {
								print HTMLSNP "</td>";
								print TABSNP "\t";
							}
							print TABSNP "\n";
							print HTMLSNP "</tr>\n";

							$nbDifferent ++ ;
							$alleleCommun ++ ;
							$nbPolyploid1 ++ ;
							$taille++;
						}
					}
					if ((($GT_poly2 =~ /^0.1$/) || ($GT_poly2 =~ /^0.1$/)) && (($GT_poly =~ /^0.0$/) || ($GT_poly eq ""))) {
						if (($SG3> $value_filter_p2) && ($SG4> $value_filter_p2)) {
							print HTMLSNP "<td style=\"border-left:3px solid black\">".$c."</td>";
							print HTMLSNP "<td class=\"border-left:3px solid black\">".$alRef."</td>";
							print HTMLSNP "<td>".$code_snp."</td>";
							print HTMLSNP "<td>".$code_snp2."</td>";
							print HTMLSNP "<td>".$FDP."/".$DP_P;
							print TABSNP $s . "\t";
							print TABSNP $c . "\t";
							print TABSNP $alRef . "\t";
							print TABSNP $code_snp . "\t";
							print TABSNP $code_snp2 . "\t";
							print TABSNP $FDP."/".$DP_P;
							if (($DP_P) != 0) {
								print HTMLSNP "<br><img src=\"".$REPimages."r1.png\" height=5 width=".sprintf("%.0f", $FDP/$DP_P*100).">";
								print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=5 width=".sprintf("%.0f", 100-$FDP/$DP_P*100).">";
								print HTMLSNP "<br>".$sub1_1." - ".$sub1_2."</td>";
								print TABSNP "," . $sub1_1." - ".$sub1_2."\t";
							}
							else {
								print HTMLSNP "</td>";
								print TABSNP "\t";
							}
							print TABSNP $FDP2."/".$DP_P2;
							print HTMLSNP "<td>".$FDP2."/".$DP_P2;
							if (($DP_P2) != 0) {
								print HTMLSNP "<br><img src=\"".$REPimages."r1.png\" height=5 width=".sprintf("%.0f", $FDP2/$DP_P2*100).">";
								print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=5 width=".sprintf("%.0f", 100-$FDP2/$DP_P2*100).">";
								print HTMLSNP "<br>".$sub2_1." - ".$sub2_2."</td>";
								print TABSNP "," . $sub2_1." - ".$sub2_2."\t";
							}
							else {
								print HTMLSNP "</td>";
								print TABSNP "\t";
							}
							print TABSNP "\n";
							print HTMLSNP "</tr>\n";
							############
							# HERE P2  #
							############
							$nbDifferent ++ ;
							$alleleCommun ++ ;
							$nbPolyploid2 ++ ;
							$taille++;
						}
					}
					# [6] [8] P1 = 1/1 ; P2 = 0/0 
					if (($GT_poly =~ /^1.1$/) && (($GT_poly2 =~ /^0.0$/) || ($GT_poly2 eq ""))) {
						print HTMLSNP "<td style=\"border-left:3px solid black\">".$c."</td>";
						print HTMLSNP "<td class=\"border-left:3px solid black\">".$alRef."</td>";
						print HTMLSNP "<td>".$code_snp."</td>";
						print HTMLSNP "<td>".$code_snp2."</td>";
						print HTMLSNP "<td>".$FDP."/".$DP_P;
						print TABSNP $s . "\t";
						print TABSNP $c . "\t";
						print TABSNP $alRef . "\t";
						print TABSNP $code_snp . "\t";
						print TABSNP $code_snp2 . "\t";
						print TABSNP $FDP."/".$DP_P;
						if (($DP_P) != 0) {
							print HTMLSNP "<br><img src=\"".$REPimages."r1.png\" height=5 width=".sprintf("%.0f", $FDP/$DP_P*100).">";
							print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=5 width=".sprintf("%.0f", 100-$FDP/$DP_P*100).">";
							print HTMLSNP "<br>".$sub1_1." - ".$sub1_2."</td>";
							print TABSNP "," . $sub1_1." - ".$sub1_2."\t";
						}
						else {
							print HTMLSNP "</td>";
							print TABSNP "\t";
						}
						print TABSNP $FDP2."/".$DP_P2;
						print HTMLSNP "<td>".$FDP2."/".$DP_P2;
						if (($DP_P2) != 0) {
							print HTMLSNP "<br><img src=\"".$REPimages."r1.png\" height=5 width=".sprintf("%.0f", $FDP2/$DP_P2*100).">";
							print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=5 width=".sprintf("%.0f", 100-$FDP2/$DP_P2*100).">";
							print HTMLSNP "<br>".$sub2_1." - ".$sub2_2."</td>";
							print TABSNP "," . $sub2_1." - ".$sub2_2."\t";
						}
						else {
							print HTMLSNP "</td>";
							print TABSNP "\t";
						}
						print TABSNP "\n";
						print HTMLSNP "</tr>\n";
							$nbDifferent ++ ;
							$alleleCommun ++ ;
							$nbPolyploid1 ++ ;
							$taille++;
					}
					if (($GT_poly2 =~ /^1.1$/) && (($GT_poly =~ /^0.0$/) || ($GT_poly eq ""))) {
						print HTMLSNP "<td style=\"border-left:3px solid black\">".$c."</td>";
						print HTMLSNP "<td class=\"border-left:3px solid black\">".$alRef."</td>";
						print HTMLSNP "<td>".$code_snp."</td>";
						print HTMLSNP "<td>".$code_snp2."</td>";
						print HTMLSNP "<td>".$FDP."/".$DP_P;
						print TABSNP $s . "\t";
						print TABSNP $c . "\t";
						print TABSNP $alRef . "\t";
						print TABSNP $code_snp . "\t";
						print TABSNP $code_snp2 . "\t";
						print TABSNP $FDP."/".$DP_P;
						if (($DP_P) != 0) {
							print HTMLSNP "<br><img src=\"".$REPimages."r1.png\" height=5 width=".sprintf("%.0f", $FDP/$DP_P*100).">";
							print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=5 width=".sprintf("%.0f", 100-$FDP/$DP_P*100).">";
							print HTMLSNP "<br>".$sub1_1." - ".$sub1_2."</td>";
							print TABSNP "," . $sub1_1." - ".$sub1_2."\t";
						}
						else {
							print HTMLSNP "</td>";
							print TABSNP "\t";
						}
						print TABSNP $FDP2."/".$DP_P2;
						print HTMLSNP "<td>".$FDP2."/".$DP_P2;
						if (($DP_P2) != 0) {
							print HTMLSNP "<br><img src=\"".$REPimages."r1.png\" height=5 width=".sprintf("%.0f", $FDP2/$DP_P2*100).">";
							print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=5 width=".sprintf("%.0f", 100-$FDP2/$DP_P2*100).">";
							print HTMLSNP "<br>".$sub2_1." - ".$sub2_2."</td>";
							print TABSNP "," . $sub2_1." - ".$sub2_2."\t";
						}
						else {
							print HTMLSNP "</td>";
							print TABSNP "\t";
						}
						print TABSNP "\n";
						print HTMLSNP "</tr>\n";
							$nbDifferent ++ ;
							$alleleCommun ++ ;
							$nbPolyploid2 ++ ;
							$taille++;
					}
				#}
			}
			
			
			#print TABSNP $s."\t".$c."\t".$alRef."\t".$code_snp."\t".$code_snp2."\t".$FDP."\t".$DP_P."\t".$FDP2."\t".$DP_P2;
			
			$ligneOK = 1 ;
		}
	
	
		if (($nbCommuns + $nbCommunHomo + $nbDifferent + $nbHomoDiff + $alleleCommun + $alleleDifferent + $alleleCommunH + $nbPolyploid1 + $nbPolyploid2) > 0 ) {
		
			if ($ligneInter == 0) {	
				print HTMLCOUNT "<td class=\"ted2\" style=\"border-right:3px solid black\">".$s."</td>";
			}
			else {
				print HTMLCOUNT "<td class=\"ted\" style=\"border-right:3px solid black\">".$s."</td>";
			}
			print TABCOUNT $s."\t";	
				
			#######################################
			if ($ligneInter == 0) { $ligneInter = 1 ; }
			else { $ligneInter = 0 ; }
			#######################################
			
			# Calcul des intervalles #
			##########################
			$taille_totale = 0 ;
			my $ref = $intervalle2{$s};
			my %hash = %$ref;
			
			foreach my $interval(keys(%hash)){
				my @pos = split(/-/,$interval);
				$taille_inter = $pos[1]-$pos[0]+1 ;
				$taille_totale = $taille_totale + $taille_inter;
			}
			$total1 = $case5 + $case1 + $case2 + $casePolyplother;
			$total2 = $case5 + $case3ou4 + $caseDiplother;
			
			# SYNTHESIS
			
			print HTMLCOUNT "<td>".$taille_totale."</td><td style=\"border-left:3px solid black\">".$taille. "</td></td>";	
			print HTMLCOUNT "<td style=\"border-left:3px solid black\">";
			print HTMLCOUNT $nbCommuns."</td><td>".$nbCommunHomo."</td><td style=\"border-left:3px solid black\">".$nbDifferent."</td><td style=\"border-left:3px solid black\">";
			print HTMLCOUNT $nbHomoDiff."</td><td>".$alleleCommun."</td><td>".$alleleDifferent."</td><td>".$alleleCommunH."</td>"; 
			print HTMLCOUNT "<td style=\"border-left:3px solid black\">".$nbPolyploid1."</td><td>".$nbPolyploid2."</td>";
			print TABCOUNT $taille_totale."\t".$taille."\t";
			print TABCOUNT $nbCommuns."\t".$nbCommunHomo."\t".$nbDifferent."\t".$nbHomoDiff."\t".$alleleCommun."\t".$alleleDifferent."\t".$alleleCommunH."\t".$nbPolyploid1."\t".$nbPolyploid2."\t";

			$nbTotGenesAna ++ ;
		
			print HTMLCOUNT "</tr>";
			print TABCOUNT "\n";
			
			$totalSize = $totalSize + $taille_totale ;
			$totalSNP = $totalSNP + $taille ;
			$totalNbPolyploid1 = $totalNbPolyploid1 + $nbPolyploid1 ;			# SNP heterozygosity for P1
			$totalNbPolyploid2 = $totalNbPolyploid2 + $nbPolyploid2 ;			# SNP heterozygosity for P2
			$totalNbCommuns = $totalNbCommuns + $nbCommuns ;					# SNP heterozygosity [P1] = [P2]
			$totalNbCommunsHomo = $totalNbCommunsHomo + $nbCommunHomo ;			# SNP homozygosity [P1] = [P2]
			$totalNbDifferent = $totalNbDifferent + $nbDifferent ;				# [P1] ne [P2]
			$totalNbAlleleCommun = $totalNbAlleleCommun + $alleleCommun ;		# Example : P1 = [A/G] ; P2 = [A]
			$totalAlleleDifferent = $totalAlleleDifferent + $alleleDifferent ;	# Example : P1 = [A/G] ; P2 = [C] or [T]
			$totalAlleleCommunH = $totalAlleleCommunH + $alleleCommunH ;		# Example : P1 = [A/G] ; P2 = [A/C]
			$totalNbHomoDiff = $totalNbHomoDiff + $nbHomoDiff ;					# Example : P1 = [A/G] ; P2 = [A/C]
		}
		
		
		
		
	}
	########## MODIF DERNIERE MINUTE ################"
	print HTMLCOUNT "<tr class=\"td3\">\n<td>";

	print HTMLCOUNT $nbTotGenesAna."<td style=\"border-left:3px solid black\">";
	print HTMLCOUNT $totalSize."</td><td style=\"border-left:3px solid black\">";
	print HTMLCOUNT $totalSNP."</td><td style=\"border-left:3px solid black\">";
	print HTMLCOUNT $totalNbCommuns."</td><td>";
	print HTMLCOUNT $totalNbCommunsHomo."</td><td>";
	print HTMLCOUNT $totalNbDifferent."</td><td style=\"border-left:3px solid black\">";
	print HTMLCOUNT $totalNbHomoDiff."</td><td style=\"border-left:3px solid black\">";
	print HTMLCOUNT $totalNbAlleleCommun."</td><td style=\"border:3px solid black\">";
	print HTMLCOUNT $totalAlleleDifferent."</td><td style=\"border:3px solid black\">";
	print HTMLCOUNT $totalAlleleCommunH."</td><td>";
	print HTMLCOUNT $totalNbPolyploid1."</td><td>";
	print HTMLCOUNT $totalNbPolyploid2."</td>";
	print HTMLCOUNT "</tr>";


	print TABCOUNT "$nbTotGenesAna\t$totalSize\t$totalSNP\t$totalNbCommuns\t$totalNbCommunsHomo\t$totalNbDifferent\t$totalNbHomoDiff\t$totalNbAlleleCommun\t$totalAlleleDifferent\t$totalAlleleCommunH\t$totalNbPolyploid1\t$totalNbPolyploid2\t";
	print TABCOUNT "\n";
 
	####################################################				
	print HTMLSNP "</table>\n";
	print HTMLSNP "</html>\n";
	close HTMLSNP ;

	print HTMLCOUNT "</table>\n";
	print HTMLCOUNT "</html>\n";
	close HTMLCOUNT ;  

	close TABSNP;
	close TABCOUNT ;

	# tie @array, 'Tie::File', $SNP_count or die ;
	# $array[82] = "<table class=\"tab2\"><th class=\"th\"  style=\"text-align:left;\">"; 
	# $array[83] = "<br>".$nbTotGenesAna." analysed genes";
	# $array[84] = "<br>".$nbTotGenesVal." with SNP validation";
	# $array[85] = "<br>Analysis performed on ".$totalSize." bp";
	# $array[86] = "<br>".$totalSNP." SNP";
	# $array[87] = "<br><img src=\"".$REPimages."5v.png\" WIDTH=20> : ".$total5." validated SNP";
	# $array[88] = "<br><br><img src=\"".$REPimages."1.png\" WIDTH=20> : ".$total11."";
	# $array[89] = "<br><img src=\"".$REPimages."2.png\" WIDTH=20> : ".$total22."";
	# $array[90] = "<br><img src=\"".$REPimages."3ou4.png\" WIDTH=20> : ".$total3ou4."";
	# $array[91] = "<br>Other SNP types : ".$totalOther."";
	# $array[92] = "<br>Heterozygosity for genome 1 : ".$totalGenome2."";
	# $array[93] = "<br>SNP between parental genomes (diploids) : ".$total512."";
	# $array[94] = "<br>SNP polyploid : ".$total534."";
	# $array[95] = "<th class=\"th\"><img src=\"".$REPimages."arbre.png\" WIDTH=400></th></table>";
}

$time2 = time ;
$tmps = $time2 - $time;
print STDOUT "\n\nTemps execution : ".$tmps."\n";
