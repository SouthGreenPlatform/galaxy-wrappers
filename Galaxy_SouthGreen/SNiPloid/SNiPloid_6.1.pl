#!/usr/local/bin/perl

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
$genome1Name = "" ;
$genome2Name = "" ;
#___________________________________
# VCF files
#-----------------------------------
$VCFpolyploid = "" ;
$VCFgenome1 = "" ;
$VCFgenome2 = "" ;
$merged_VCF = "" ; # Polyploid + Genome1 + Genome 2
#___________________________________
# Depth of Coverage File
#-----------------------------------
$DOCpolyploid = "" ;
$DOCgenome1 = "" ;
$DOCgenome2 = "" ;
$merged_DOC = "" ; # Polyploid + Genome1 + Genome 2
#___________________________________
# Depth for each sample
#-----------------------------------
$depthPolyploid = 0 ;
$depthGenome1 = 0 ;
$depthGenome2 = 0 ;
#___________________________________
# Output Files
#-----------------------------------
$SNP_csv = "SNP_tab.xls";
$SNP_html = "SNP_view.html";
$SNP_count = "SNP_synthesis_tab.html";
$SNP_count_csv = "SNP_synthesis_tab.xls";
#___________________________________
# Other parameters
#-----------------------------------
$enableLowQuality = 0 ;  #default value for enable quality SNP = only PASS SNP are considered
$ref = 0 ; # default parameter = extern
$merged = 1 ; # default value

$REPimages = "http://marmadais.cirad.fr:7070/static/images/";

=POD
Add option for "Heterozygosity"
 ° Enable "heterozygosity" for genome 1 (reference intern) - not necessary...
 ° Enable "heterozygosity" for genome 1 and genome 2 (reference extern)
=cut

GetOptions (
	"pn=s"	=>	\$polyploidName,
	"gn1=s"	=>	\$genome1Name,  
	"gn2=s"	=>	\$genome2Name,
	"vp=s"	=>	\$VCFpolyploid,
	"vg1=s"	=>	\$VCFgenome1,
	"vg2=s"	=>	\$VCFgenome2,
	"vm=s"	=>	\$merged_VCF,
	"cpp=s"	=>	\$DOCpolyploid,
	"cg1=s"	=>	\$DOCgenome1,
	"cg2=s"	=>	\$DOCgenome2,
	"cm=s"	=>	\$merged_DOC,
	"dp=i"	=>	\$depthPolyploid,
	"dg1=i"	=>	\$depthGenome1,
	"dg2=i"	=>	\$depthGenome2,
	"oc=s"	=>	\$SNP_csv,
	"oh=s"	=>	\$SNP_html,
	"ohs=s"	=>	\$SNP_count,
	"ocs=s"	=>	\$SNP_count_csv,
	"elq=i"	=>	\$enableLowQuality,
	"ref=i"	=>	\$ref,
# h = i	= >	\ $ heterozygosity ,
);


# Validation - Samples names

print STDOUT "\nPname : ".$polyploidName ;

if ($polyploidName eq "") {	die ("*** /!\\ ERROR : Missing name for polyploid - You have to specify a name for the polyploid species [--pn \"polyploid_name\"] $!") ;}
if ($genome1Name eq "") {	die ("*** /!\\ ERROR : Missing name for genome 1 - You have to specify a name for the genome 1 species") ;}
if ($genome2Name eq "") {	die ("*** /!\\ ERROR : Missing name for genome 2 - You have to specify a name for the genome 2 species") ;}
# Validation - depth
if ($depthPolyploid == 0) {			die ("*** /!\\ ERROR : Missing depth information for polyploid");}
if ($depthGenome1 == 0) {			die ("*** /!\\ ERROR : Missing depth information for genome 1");}
if ($ref == 0 && $depthGenome2 == 0) {die ("*** /!\\ ERROR : Missing depth information for genome 2");}

%intervalle1 ;
%intervalle2 ;
%snp = () ;
my %snp_final ;
$nbTotGenes = 0 ;
$nbTotGenesVal = 0 ;
$nbTotGenesAna = 0 ;


if (($merged_VCF eq "") && ($merged_DOC eq "")) {
	$merged = 0 ; #if no merged files
}

$time = time ;

=POD PART 1 : CREATING COMMON INTERVALS 
=cut
if ($merged == 1) {
	&Intervall_part1($merged_DOC) ;
}	
else {
	&Intervall_part1($DOCpolyploid) ;
	&Intervall_part2($DOCgenome1,$depthGenome1) ;
	if ($ref == 0) { #genome2 => no parental genome as reference
		&Intervall_part2($DOCgenome2,$depthGenome2) ;
	}
}
=POD PART 2 and 3 : CREATING SNP TAB AND OUTPUTS
=cut
if ($merged == 1) {
	&VCF_merged_Analysis($merged_VCF);
	
	if ($ref == 1) { # Reference = one of two parental genomes
	print STDOUT "\n$merged_VCF";
		&int_output ;
	}
	else { 
		&ext_output ;
	}
}
else {
	&VCF_Analysis($VCFpolyploid);
	if ($ref == 1) { # Reference = one of two parental genomes
		&VCF_Analysis($VCFgenome1);
		&intro_output ;
		&int_output ;
	}
	else { # Extern Reference
		&VCF_Analysis($VCFgenome1);
		&VCF_Analysis($VCFgenome2);
		&intro_output ;
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
	#Locus	Total_Depth	Average_Depth_sample	Depth_for_BD54	Depth_for_cat	Depth_for_eug
	#Locus	Total_Depth	Average_Depth_sample	Depth_for_Eug
	if ($DOC[0] =~ /^Locus\tTotal_Depth\tAverage_Depth_sample(.+)$/) {		
		
			@order = split(/\t/,$1);
			@order = split(/\tDepth_for_/,$1);
			
			# print STDOUT "\n1:$order[1]";
			# print STDOUT "\n2:$order[2]";
			# print STDOUT "\n3:$order[3]";
			switch ($order[1]) {
				case "$polyploidName" { $indicePolyploid1 = 3 }
				case "$genome1Name" { $indiceGenome1 = 3 }
				if ($VCFgenome2 ne ""){ # Genome 2 present
					case "$genome2Name" { $indiceGenome2 = 3 }
				}
			}
			switch ($order[2]) {
				case "$polyploidName" { $indicePolyploid1 = 4 }
				case "$genome1Name" { $indiceGenome1 = 4 }
				if ($VCFgenome2 ne ""){
					case "$genome2Name" { $indiceGenome2 = 4 }
				}
			}
			if ($order[3] ne "") {
				switch ($order[3]) {
					case "$polyploidName" { $indicePolyploid1 = 5 }
					case "$genome1Name" { $indiceGenome1 = 5 }
					if ($VCFgenome2 ne ""){
						case "$genome2Name" { $indiceGenome2 = 5 }
					}
				}
			}
		}
		
	foreach $line(@DOC) {
		if ($line ne $DOC[0]) {
			@ligne = split(/\t/ , $line);
			@position = split(/:/ , $ligne[0]);
			$name_gene = $position[0] ;
			
			if ($merged == 0) { # 1st File	 - Polyploid
				$depthcov = $ligne[1] ;
				if ($name_gene){
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
						$rec = 0 ;
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
	my %tab;
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
			if ($interval eq "SH3.74_45-4E.13") {
				#print STDOUT "\n1  ".$tab{$interval}{$i};
				#print STDOUT "\n2  $interval - $i";
			}
			$final = 2 ;
			$rec = 0 ;
			($debut,$fin) = split(/-/,$intervall);
			for ($i=$debut; $i <=$fin; $i++) {
				if ($interval eq "SH3.74_45-4E.13") {
					#print STDOUT "\n1  ".$tab{$interval}{$i};
					#print STDOUT "\n2  $interval - $i";
				}
				#print STDOUT "\n".$tab{$interval}{$i};
				#print STDOUT "\n$interval - $i";
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
						#print STDOUT "\n$val_fin";
						$intervalle2{$name_gene}{$val_fin} = "ok" ;
					}
					$rec = 0 ;
				}
			}
			if ($final == 0) {
				$val_fin = $val_deb.$fin ;
				#print STDOUT "\n$val_fin";
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
	
		foreach my $intervall(sort (keys(%intervalls))){ 
			#if ($intervalle2{$s}{$snip} ne "LQ"){
				# print STDOUT "\n$interval \t".$intervall ;
			#}
		}
	}	
	
	
	return (%intervalle2) ;
}

sub VCF_Analysis {
	%snp_final = () ;
	#%snp = () ;
	my(@args) = @_;
		open (TABSNP, "$args[0]") or die ("ERROR : file $args[0] don't exists");
		@VCF = <TABSNP> ;
	close TABSNP ;
	foreach $line(@VCF){
		if ($line !~ /^#/){
			@infos_line = split(/\t/,$line) ;
			$gene = $infos_line[0];
			$position = $infos_line[1];
			$ref_allele = $infos_line[3];
			$alt_allele = $infos_line[4];
			$snp_code = "[$ref_allele/$alt_allele]";
			$quality_of_snp = $infos_line[6];
			$depth_recuperation = $infos_line[7];
			$alleles = $infos_line[9];
		
			($GT,$AD,$FDP,$GQ,$PL) = split(":",$alleles);
			($AC,$AF,$AN,$DP,$DS,$Dels,$HRun,$HaplotypeScore,$MQ,$MQ0,$QD,$SB,$sumGLbyD) = split(";",$depth_recuperation);
			
			($sub1,$sub2) = split(",",$AD);
			$somme = $sub1 + $sub2 ;
			#print STDOUT "\n".$gene." ".$position." ".$sub1." ".$sub2;
			if ($somme == 0 ) {
				print STDOUT "ERROR : Cannot calculate ratio for ".$gene." [pos:".$position."]\n\"".$line."\"";
				die ("ERROR : Cannot calculate ratio for ".$gene." [pos:".$position."]\n\"".$line."\"");
			}
			else {
				$ratio = ($sub1/($sub1+$sub2))*100;
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
			if ($enableLowQuality eq "true") {
				
				if ($test_inside_interval == 1 ){ #
					if ($args[0] eq $VCFpolyploid) { # Polyploid
						
						$snp{$gene}{$position} = $snp_code."\t".$ratio."\t".$GT."\t".$DP[1]."-".$FDP ;
						#print STDOUT "\nPOLY : ".$snp{$gene}{$position};
					}
					else { 
						if ($args[0] eq $VCFgenome1) { # genome1
							if (exists $snp{$gene}{$position}) { # if polyploid SNP
								$snp{$gene}{$position} = $snp{$gene}{$position}."\t".$snp_code."\t".$GT ;
							}
							else { # if no polyploid SNP, key is empty
								$snp{$gene}{$position} = $ref_allele."\t\t\t\t".$snp_code."\t".$GT ;
							}
						}
						else { #genome2
							if ($args[0] eq $VCFgenome2) {
								if (exists $snp{$gene}{$position}) { # if polyploid SNP
									$snp{$gene}{$position} = $snp{$gene}{$position}."\t".$snp_code."\t".$GT;
								}
								else { # if no polyploid SNP and no genome1, key is empty
									$snp{$gene}{$position} = $ref_allele."\t\t\t\t".$ref_allele."\t\t".$snp_code."\t".$GT ;
								}
							}
						}
					}
				}
			}
			# ONLY PASS SNP CONSIDERED
			else {
				# MULTIALLELIX
				if (($test_inside_interval == 1 ) && ($quality_of_snp eq "PASS") && ($snp{$gene}{$position} ne "LQ")){ #
					
				#print STDOUT "\n$VCFgenome1";
					if ($AC eq "AC=2"){
						if ($DP[1] == ($sub1+$sub2)) {
						}
						else {
							$percent = (($DP[1]-($sub1+$sub2))/($DP[1]))*100 ;
							print STDOUT "\n".$gene."\t".$position."\t".$ref_allele."\t".$alt_allele."\t".$percent."\t".$AC."\t".$GT;
						}	
					}
				
					if ($args[0] eq $VCFpolyploid) { # Polyploid
						$snp{$gene}{$position} = $snp_code."\t".$ratio."\t".$GT."\t".$DP[1]."-".$FDP ;
					}
					else { 
						if ($args[0] eq $VCFgenome1) { # genome1
							#print STDOUT "\nG1";
							if (exists $snp{$gene}{$position}) { # if polyploid SNP
								$snp{$gene}{$position} = $snp{$gene}{$position}."\t".$snp_code."\t".$GT ;
								#print STDOUT "\nEXIST";
							}
							else { # if no polyploid SNP, key is empty
								$snp{$gene}{$position} = $ref_allele."\t\t\t\t".$snp_code."\t".$GT ;
								#print STDOUT "\n-------";
							}
						}
						else { #genome2
							if ($args[0] eq $VCFgenome2) {
								if (exists $snp{$gene}{$position}) { # if polyploid SNP
									$snp{$gene}{$position} = $snp{$gene}{$position}."\t".$snp_code."\t".$GT;
								}
								else { # if no polyploid SNP and no genome1, key is empty
									$snp{$gene}{$position} = $ref_allele."\t\t\t\t".$ref_allele."\t\t".$snp_code."\t".$GT ;
								}
							}
						}
					}
				}
				else {
					if ($quality_of_snp ne "PASS") {  
						$snp{$gene}{$position} = "LQ";
					}
				}
			}
			################################################################################################################################	
		}
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
	
	foreach my $s(sort(keys(%snp_final))){
		my $ref = $snp_final{$s};
		my %hash = %$ref;
		foreach my $snip(keys(%hash)){
			#print STDOUT "\n$s - $snip - $snp_final{$s}{$snip}";
		}
	}
	
	return (%snp_final) ;
}
sub VCF_merged_Analysis {
	
	my(@args) = @_;
	open (TABSNP, "$args[0]") or die ("Pbm a l'ouverture du fichier : $VCF");
		@VCF = <TABSNP> ;
	close TABSNP ;
	
	foreach $line(@VCF){
			#Determine genomes order
			if ($line =~ /^#CHROM.+FORMAT\t(.+)$/) {
				@order = split(/\t/,$1);
				switch ($order[0]) {
					case "$polyploidName" { $indicePolyploid1 = 9 }
					case "$genome1Name" { $indiceGenome1 = 9 }
					if ($VCFgenome2 ne "") { # Genome 2 present
						case "$genome2Name" { $indiceGenome2 = 9 }
					}
				}
				switch ($order[1]) {
					case "$polyploidName" { $indicePolyploid1 = 10 }
					case "$genome1Name" { $indiceGenome1 = 10 }
					if ($VCFgenome2 ne ""){ # Genome 2 present
						case "$genome2Name" { $indiceGenome2 = 10 }
					}
				}
				if ($order[2] ne "") {
					switch ($order[2]) {
						case "$polyploidName" { $indicePolyploid1 = 11 }
						case "$genome1Name" { $indiceGenome1 = 11 }
						if ($VCFgenome2 ne ""){ # Genome 2 present
							case "$genome2Name" { $indiceGenome2 = 11 }
						}
					}
				}
			}
		
		if ($line !~ /^#/){
			@infos_line = split(/\t/,$line) ;
			$gene = $infos_line[0];
			$position = $infos_line[1];
			$ref_allele = $infos_line[3];
			$alt_allele = $infos_line[4];
			$snp_code = "[$ref_allele/$alt_allele]";
			$quality_of_snp = $infos_line[6];
			$depth_recuperation = $infos_line[7];
			$allelesP = $infos_line[$indicePolyploid1];
			$allelesG1 = $infos_line[$indiceGenome1];
			$allelesG2 = $infos_line[$indiceGenome2];
			($AC,$AF,$AN,$DP,$DS,$Dels,$HRun,$HaplotypeScore,$MQ,$MQ0,$QD,$SB,$sumGLbyD) = split(";",$depth_recuperation);
			# POLYPLOID
			if (($allelesP !~ /^\.\/\.$/)&&($allelesP !~ /^0\/0/)) {
				($GT_P,$AD_P,$DP_P,$GQ_P,$PL_P) = split(":",$allelesP);
				($sub1_P,$sub2_P) = split(",",$AD_P);
				$somme_P = $sub1_P + $sub2_P ;
				#print "\n$sub1_P - $sub2_P";
				$ratio_P = ($sub1_P/($sub1_P+$sub2_P))*100;
				$ratio_P = sprintf("%.0f", $ratio_P);
			}
			# GENOME 1
			if ($allelesG1 !~ /^\.\/\.$/) {
				($GT_G1,$AD_G1,$DP_G1,$GQ_G1,$PL_G1) = split(":",$allelesG1);
			}
			# GENOME 2
			if ($VCFgenome2 ne "") {
				if ($allelesG2 !~ /^\.\/\.$/) {
					($GT_G2,$AD_G2,$DP_G2,$GQ_G2,$PL_G2) = split(":",$allelesG2);
				}
			}
			else {$GT_G2 = ""} ;
			
			$test_inside_interval = 0 ;
		
			my $ref = $intervalle1{$gene};
			my %hash = %$ref;
			
			foreach my $interval(keys(%hash)){
				my @pos = split(/-/,$interval);
				if ($position >= $pos[0] && $position <= $pos[1] ){
					$test_inside_interval = 1 ;
				}
			}
			if (($allelesG1 !~ /^\.\/\.$/)&&($allelesP !~ /^\.\/\.$/)) {
				if ($enableLowQuality eq "true") {
					if ($test_inside_interval == 1 ){ #
						@DP = split ("=",$DP) ;
						$snp_final{$gene}{$position} = $snp_code."\t".$ratio_P."\t".$GT_P."\t".$DP_P."\t".$snp_code."\t".$GT_G1."\t".$snp_code."\t".$GT_G2 ;
					}
				}
				else {
					if (($test_inside_interval == 1 )&& ($quality_of_snp eq "PASS")){ #
						@DP = split ("=",$DP) ;
						$snp_final{$gene}{$position} = $snp_code."\t".$ratio_P."\t".$GT_P."\t".$DP_P."\t".$snp_code."\t".$GT_G1."\t".$snp_code."\t".$GT_G2 ;
					}
				}
			}
		}
	}
	return (%snp_final);
}

sub intro_output {
	#print STDOUT "ici ?";
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
print HTMLCOUNT " border-width:0px; \n";
print HTMLCOUNT " border-style:solid; \n";
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

print HTMLCOUNT ".auto-style1 {";
print HTMLCOUNT "	font-weight: normal;";
print HTMLCOUNT "	font-size: x-small;";
print HTMLCOUNT "}";

print HTMLCOUNT "</style>\n";

###################################################################################################################################################################################

print HTMLSNP "</head>\n";
print HTMLSNP "<center><img src=\"".$REPimages."SNiPloid7.png\" WIDTH=250></center>";
print HTMLSNP "<center><img src=\"".$REPimages."arbre.png\" WIDTH=450></center>";
print HTMLSNP "<p>\n";
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print HTMLCOUNT "</head>\n";
print HTMLCOUNT "<center><img src=\"".$REPimages."SNiPloid7.png\" WIDTH=250></center>";
print HTMLCOUNT "<p>\n";
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print HTMLSNP "<body>\n";
print HTMLSNP "<center><h3><font face=\"calibri\">Result of SNP comparison of a Polyploid and its Parental Genomes (Genome 1 and Genome 2 as reference)</font></h3></center>";
print HTMLSNP "<p>\n";
# COLUMNS - HTMLSNP SNP VIEW
print HTMLSNP "<table border=\"1\" border cellpadding=\"5\" style=\"text-align:center;\"> \n";
print HTMLSNP "<tr>\n";
print HTMLSNP "<th>Gene</th>"; 																					# (1) Gene
print HTMLSNP "<th>Position</th>"; 																				# (2) Position
print HTMLSNP "<th>Polyploid<br><span class=\"auto-style1\">".$polyploidName."</span></th>"; 					# (3) Polyploid
print HTMLSNP "<th>Genome 1<br><span class=\"auto-style1\">".$genome1Name."</span></th>"; 						# (4) Genome 1
print HTMLSNP "<th>Genome 2<br><span class=\"auto-style1\">".$genome2Name."</span></th>"; 						# (5) Genome 2
print HTMLSNP "<th>Validation</th>"; 																			# (6) Validation
print HTMLSNP "<th>Ratio (%)<br><span class=\"auto-style1\">".$genome2Name." : ".$genome1Name."</span></th>";	# (7) Ratio
print HTMLSNP "<th>Filtered<br>Depth</th>"; 																		# (8) Filtered Depth
print HTMLSNP "<th>Total<br>Depth</th>"; 																			# (9) Total Depth
print HTMLSNP "</tr>\n";
# Entête fichier SNP VIEW TAB
print TABSNP "Gene\t";																							# (1) Gene
print TABSNP "Position\t";																						# (2) Position
print TABSNP "Polyploid : ".$polyploidName."\t";																# (3) Polyploid
print TABSNP "Genome 1 : ".$genome1Name."\t";																	# (4) Genome 1
print TABSNP "Genome 2 : ".$genome2Name."\t";																	# (5) Genome 2
print TABSNP "Validation\t";																					# (6) Validation
print TABSNP "Ratio (%) ".$genome2Name."\t";																	# (7) Ratio Genome 1
print TABSNP "Ratio (%) ".$genome1Name."\t";																	# (8) Ratio Genome 2
print TABSNP "Filtered Depth\t";																				# (9) Filtered Depth
print TABSNP "Total Depth\n";																					# (10) Total Depth

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print HTMLCOUNT "<body>\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";
print HTMLCOUNT "<center><h3><font face=\"calibri\">Synthesis of the analysis</font></h3></center>";
print HTMLCOUNT "<p>\n";
# COLUMNS - HTML SYNTHESIS
print HTMLCOUNT "<table border=\"1\" border cellpadding=\"5\" style=\"text-align:center;\"> \n";
print HTMLCOUNT "<tr>\n";
print HTMLCOUNT "<th>Gene</th>";																							# (1) Gene
print HTMLCOUNT "<th>Interval Size<br>Analysed (pb)</th>";																	# (2) Interval Size Analyzed
print HTMLCOUNT "<th>nb Positions<br>with SNP</th>";																		# (3) nB Positions With SNP
print HTMLCOUNT "<th><img src=\"".$REPimages."1.png\" height=30></th>";														# (4) [1]
print HTMLCOUNT "<th><img src=\"".$REPimages."2.png\" height=30></th>";														# (5) [2]
print HTMLCOUNT "<th><img src=\"".$REPimages."3ou4.png\" height=30></th>";													# (6) [3 or 4]
print HTMLCOUNT "<th>Other</th>";																							# (7) [other]
print HTMLCOUNT "<th>SNP<br>Heterozygosity<br>for Genome 1<br><span class=\"auto-style1\">".$genome1Name."</span></th>";	# (8) Heterozygosity For Genome 1	
print HTMLCOUNT "<th>SNP Diploids</th>";																					# (9) SNP Diploids
print HTMLCOUNT "<th>SNP Polyploid</th>";																					# (10) SNP Polyploid
print HTMLCOUNT "<th><img src=\"".$REPimages."5.png\" height=30></th>";														# (11) [5]
print HTMLCOUNT "<th>Ratio (%)<br><span class=\"auto-style1\">".$genome2Name." : ".$genome1Name."</span></th>";				# (12) Ratio %
print HTMLCOUNT "</tr>\n";
# Entête HTMLSNP Synthesis
print TABCOUNT "Gene\t";																									# (1) Gene
print TABCOUNT "Interval Size Analysed (pb)\t";																				# (2) Interval Size Analysed (pb)
print TABCOUNT "nb SNP positions\t";																						# (3) nb SNP positions
print TABCOUNT "1\t";																										# (4) [1]
print TABCOUNT "2\t";																										# (5) [2]
print TABCOUNT "3 or 4\t";																									# (6) [3 or 4]
print TABCOUNT "other\t";																									# (7) [other]
print TABCOUNT "SNP Heterozygosity Genome 1\t";																				# (8) SNP Heterozygosity Genome 1
print TABCOUNT "SNP Diploids\t";																							# (9) SNP Diploids
print TABCOUNT "SNP Polyploid\t";																							# (10) SNP Polyploid
print TABCOUNT "5\t";																										# (11) [5]
print TABCOUNT "Ratio (%) ".$genome2Name."\t";																				# (12) Ratio (%) Genome 2
print TABCOUNT "Ratio (%) ".$genome1Name."\n";																				# (13) Ratio (%) Genome 1
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
		$caseOther = 0; 
		$casePolyplother = 0 ;	# SNP chez le polyploide dans un cas "Other"
		$caseDiplother = 0 ; 	# SNP chez les diploides dans un cas "Other"
	
	
	
		#Moyenne Ponderee
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
		
			### Recuperation des informations ###
			
			#print STDOUT "\n\n\n".$snp_final{$s}{$c} ;
			($code_snp,$ratio,$GT_poly,$DP_P,$code_G1,$GT_G1) = split(/\t/,$snp_final{$s}{$c});
			($DP_P, $FDP) = split(/-/, $DP_P);
			print STDOUT "\n($code_snp:$GT_poly) - ($code_G1:$GT_G1)" ;
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
				print HTMLSNP "<td></td><td></td><td></td><td></td>";
			}	
			print TABSNP $s."\t".$c."\t".$code_snp."\t".$code_G1."\t".$alRef."\t";		
	
	
	
	
			if (($GT_poly =~ /^0.1$/)||($GT_poly =~ /^1.0$/)) { #	SNP POLYPLOID - [0/1] [0|1] [1|0]
				# Moyenne du Ratio -----------------------------------
				$ratio = sprintf("%.0f", $ratio);
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
							print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5.png\" height=30></td>";
						}
						else {
							print HTMLSNP "<td class=\"ted\">".$ratio.":".$ratio2."<br>";							# RATIO %
							print HTMLSNP "<img src=\"".$REPimages."r1.png\" height=10 width=".$ratio.">";			# IMG RATIO 1
							print HTMLSNP "<img src=\"".$REPimages."r2.png\" height=10 width=".$ratio2.">";			# IMG RATIO 2
							print HTMLSNP "<td class=\"ted\">".$FDP."</td>";
							print HTMLSNP "<td class=\"ted\">".$DP_P."</td>";
							print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5.png\" height=30></td>";
						}
						print TABSNP "OK\t".$ratio."\t".$ratio2."\t".$FDP."\t".$DP_P."\t5";
						$case5 ++ ;
					}
					else { # Other 0.1 - 1.1 (O GA A)					# Other [SNP DIPLO + SNP POLY]
						print HTMLSNP "<td class=\"final\">other</td>";
						print TABSNP "\t\t\t\t\tother";
						$caseOther ++ ;
						$casePolyplother ++ ;
						$caseDiplother ++ ;
					}
				}
				else {
					if (($GT_G1 =~ /^0.0$/)||($GT_G1 =~ /^$/)){			# 3 ou 4
						print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."3ou4.png\" height=30></td>";
						print TABSNP "\t\t\t\t\t3or4";
						$case3ou4 ++ ;
					}
					else { #0/1											# heterozygosity G1
						print HTMLSNP "<td class=\"final\">Heterozygosity<br>Genome 1</td>";
						print TABSNP "\t\t\t\t\theterozygosity for genome 1";
						$nbGenome2 ++ ;
						$casePolyplother ++ ;
					}
				}
			}
				
			if (($GT_poly =~ /^1.1$/)) { #	POLYPLOID NE REFERENCE - [1/1]
				if (($GT_G1 !~ /^0.0$/) && ($GT_G1 !~ /^1.1$/)){ # SNP Genome1 intra [0/1] [0|1] [1|0]
					print HTMLSNP "<td class=\"final\">other</td>";			
					print TABSNP "\t\t\t\t\tother";
					$nbGenome2 ++ ;
				}
				else { # Pas de SNP Genome1 [0/0] [0|0] [1/1] [1|1] SNP entre genome1 et genome2
					if ($GT_G1 =~ /^0.0$/){									# Other [NOTHING]
						print HTMLSNP "<td class=\"final\">other</td>";
						print TABSNP "\t\t\t\t\tother";
						$caseOther ++ ;
					}
					else {
						if ($code_G1 eq $alAltP) {							# 2
							print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."2.png\" height=30></td>";
							print TABSNP "\t\t\t\t\t2";
							$case2 ++ ;
						}
						else {											# Other [SNP DIPLO]
							print HTMLSNP "<td class=\"final\">other</td>";
							print TABSNP "\t\t\t\t\tother";
							$caseOther ++ ;
							$caseDiplother ++ ;
						}
					}
				}
			}
			
			if (($GT_poly =~ /^0.0$/)||($GT_poly =~ /^$/)) { #POLYPLOID == REFERENCE - [0|0]
				####################################
				if (($GT_G1 !~ /^0.0$/) && ($GT_G1 !~ /^1.1$/)){ # SNP Genome1 intra [0/1] [0|1] [1|0]
					print HTMLSNP "<td class=\"final\">other</td>";
					print TABSNP "\t\t\t\t\tother";
					$nbGenome2 ++ ;
				}
				else { # Pas de SNP Genome1 [0/0] [0|0] [1/1] [1|1] SNP entre genome1 et genome2
					if ($GT_G1 =~ /^1.1$/){
						print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."1.png\" height=30></td>";
						print TABSNP "\t\t\t\t\t1";
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
		$total1 = $case5 + $case1 + $case2 + $casePolyplother;
		$total2 = $case5 + $case3ou4 + $caseDiplother;
		
		# SYNTHESIS
		print HTMLCOUNT "<td>".$taille_totale."</td><td style=\"border-left:3px solid black\">".$taille."</td></td>";	
		print HTMLCOUNT "<td style=\"border-left:3px solid black\">";
		print HTMLCOUNT $case1."</td><td>".$case2."</td><td>".$case3ou4."</td><td>";
		print HTMLCOUNT $caseOther."</td><td>".$nbGenome2."</td><td style=\"border-left:3px solid black\">".$total1."</td>";
		print HTMLCOUNT "<td style=\"border-left:3px solid black\">".$total2."</td><td>".$case5."</td>";
		print TABCOUNT $taille_totale."\t".$taille."\t";
		print TABCOUNT $case1."\t".$case2."\t".$case3ou4."\t".$caseOther."\t".$nbGenome2."\t".$total1."\t".$total2."\t".$case5."\t";

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
	print HTMLCOUNT $totalOther."</td><td>";
	print HTMLCOUNT $totalGenome2."</td><td style=\"border-left:3px solid black\">";
	print HTMLCOUNT $total512."</td><td style=\"border:3px solid black\">";
	print HTMLCOUNT $total534."</td><td style=\"border:3px solid black\">";
	print HTMLCOUNT $total5."</td><td>";
	print HTMLCOUNT $nbTotGenesVal."</td>";

	print HTMLCOUNT "</tr>";


	print TABCOUNT "$nbTotGenesAna\t$totalSize\t$totalSNP\t$total11\t$total22\t$total3ou4\t$totalOther\t$totalGenome2\t$total512\t$total534\t$total5\t";
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
	$array[82] = "<table class=\"tab2\"><th class=\"th\"  style=\"text-align:left;\">"; 
	$array[83] = "<br>".$nbTotGenesAna." analysed genes";
	$array[84] = "<br>".$nbTotGenesVal." with SNP validation";
	$array[85] = "<br>Analysis performed on ".$totalSize." bp";
	$array[86] = "<br>".$totalSNP." SNP";
	$array[87] = "<br><img src=\"".$REPimages."5.png\" WIDTH=20> : ".$total5." validated SNP";
	$array[88] = "<br><br><img src=\"".$REPimages."1.png\" WIDTH=20> : ".$total11."";
	$array[89] = "<br><img src=\"".$REPimages."2.png\" WIDTH=20> : ".$total22."";
	$array[90] = "<br><img src=\"".$REPimages."3ou4.png\" WIDTH=20> : ".$total3ou4."";
	$array[91] = "<br>Other SNP types : ".$totalOther."";
	$array[92] = "<br>Heterozygosity for genome 1 : ".$totalGenome2."";
	$array[93] = "<br>SNP between parental genomes (diploids) : ".$total512."";
	$array[94] = "<br>SNP polyploid : ".$total534."";
	$array[95] = "<th class=\"th\"><img src=\"".$REPimages."arbre.png\" WIDTH=400></th></table>";

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
		
		print STDOUT "\n $c $GT_poly $GT_G1 $GT_G2";
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
							print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5.png\" height=30></td>";
							print TABSNP "$alRef";
							print TABSNP "OK";
							print TABSNP "[$ratio/$ratio2;]";
							print TABSNP "$DP_P;";
							print TABSNP "<td class=\"final\"><img src=\"".$REPimages."5.png\" height=30></td>";
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
							print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5.png\" height=30></td>";
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
			print STDOUT "GTPOLY:$GT_poly@";
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
							print HTMLSNP "<td class=\"final\"><img src=\"".$REPimages."5.png\" height=30></td>";
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



$time2 = time ;
$tmps = $time2 - $time;
print STDOUT "\n\nTemps execution : ".$tmps."\n";