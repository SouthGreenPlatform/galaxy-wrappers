#!/usr/bin/perl

# ------------------------------------------------------------------
# Copyright (C) 2004 INRA <eugene@ossau.toulouse.inra.fr>
#
# This program is open source; you can redistribute it and/or modify
# it under the terms of the Artistic License (see LICENSE file).
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#
# You should have received a copy of Artistic License along with
# this program; if not, please see http://www.opensource.org
#
# $Id: egn_evalpred.pl,v 1.1 2005/03/17 13:39:52 cros Exp $
# ------------------------------------------------------------------
# File:     egn_evalpred.pl
# Contents: evaluation of a gene annotation
# ------------------------------------------------------------------


$setoffset=0;
$offset=0;
$borneg=0;
$borned=0;
$sortie=1; # 0-> short, 2-> detaillee
$name2="";
$ntlevel=0;
$NOMduPROG = $0;
@CMD="$NOMduPROG "."@ARGV";

sub description{
  die <<EOF

############################################################################ 
# Ce prog evalue la qualite de predictions d'EuGene sur des seq nucleiques,
# gere les predictions reverses et plsrs genes ou predictions par sequence. 
#
# arguments :
#
#     - 1 : fichier coordonnees reelles des exons
# -les sequences annotees sont dans le meme ordre que les predictions
# -un espace devant chaque coordonnee exon
# -un signe moins devant la coord si le sens est reverse
# -un gene par ligne
#  (retour chariot apres le dernier exon d'un gene)
# -une ligne vide apres chaque sequence (derniere comprise)
#  (ligne vide = seulement retour chariot "\n")
# exemple format coord:
#
#  13 56 1875 1945 6211 6574
#
#  125 653 2789 3052
#  -3238 -3267 -3289 -3302
#
#  76 331
#  ...
#
#     - 2 : coordonnees predites.
# sortie standard d'EuGeneAS 
# ex: EuGeneAS *.fasta > out, ou EuGeneAS `cat seqfiles.list` > out
# exemple format predictions:
#
#  Init  +   1   13   13   +1   +1   0     14   1.0
#  Term  +  267 1282 1016  +2   +2   266   1283  1.
#  Init  + 1394 1437 ...
#  Intr ...
#  Term ...
#
#  Init  + 1111 1120  88   +1   +3  1111  1120 1.0
#  Sngl  - 1890 1945  ...
# ...
# depuis mai 2002: tolere les lignes "Utr5" et "Utr3"
#
#     - 3, 4, 5, 6 : arguments facultatifs (ordre indifferent).
# -oN
# N=offset, qui definit des bornes de part et d'autre des 
# coordonnees des genes reels au dela desquelles les predictions
# ne seront pas prises en compte.(cf evaluation avec araset, ou offset=300)
# -pX
# affichage de la sortie, X=s pour short (utilise pour l optimisation),
# et X=l pour long (equivaut a sortie detaillee) (par defaut:intermediaire)
# -fX
# sert pour sortie detaillee, le nom des sequences proviendra non pas du
# fichier predictions (comme par defaut, preleve sur la 1ere colonne),
# mais du fichier X (un nom par ligne, pouvant etre l ID, le chemin,...)
# -nt
# option ajoutant les valeurs sensibilite/specificite au niveau nucleotide,
# calculees sur tte la seq, ou la region contenant: gene(s) + offset
# attention! codee dans l'urgence, la fonction coute bcp de temps de calcul!
#
# ex: $NOMduPROG araclean.exons.coord eugene.outpred -o300 -pl -fID.list
#
############################################################################

EOF
}

############################################################################ 
sub usage{
  die <<EOF

usage : 
$NOMduPROG exons_coord_file pred_file [-o{offset}] [-p{s|l}] [-f{IDfile}]
    -o  : offset (integer), only predictions in the region 
          "real_gene_frontiers +/- offset" are evaluated.
          DEFAULT=INF -> all predictions are considerate
    -ps : short output (used for optimisation scripts)
    -pl : long output  (used for expertised analyses)
          (default = intermediaire)
    -f  : name (or ID, path...) file (one name/line), used if -pl is active
    -nt : compute accuracy at nucleotide level - warning: time expensive
"-" arguments are optionnal, must be adjacent to the "-", order no matters
ex : $NOMduPROG seq.exons.coord eugene.outpred -o300 -pl -fseq.ID.list

EOF
}

############################################################################ 
# Fonction qui lit une sequence fasta et renvoie sa longueur 
sub seqlength{
  my$filename = $_[0];
  my$len=0;
  open(FASTASEQ,"$filename") || die "Can't open seq file $filename in seqlength\n";
  if (!(<FASTASEQ> =~ /^\s?>/ )) { 
    die "Fasta format is needed for length calcul for nt level in seq $filename\n";
  }
  while (<FASTASEQ>){
    chomp;
    $len += length($_);
  }
  return $len;
}

############################################################################ 
# Fonction d'evaluation au niveau gene qui,
# pour une sequence donnee, regarde si chaque
# gene reel a une prediction parfaitement egale
sub evaluation_niveau_gene{
  my@COORD2D=@{@_[0]}; # 1er argument
  my@PRED2D=@{@_[1]};  # 2eme argument

  # nbre de genes reels/predits
  $nGR+=($#COORD2D + 1);
  $nGP+=($#PRED2D + 1);

  # Vrais positifs:
  # pour chaque gene reel de la sequence
  foreach $pointeur_gene_reel (@COORD2D){
    my@genereelcoord = @{$pointeur_gene_reel};
    foreach $pointeur_gene_pred (@PRED2D){
      if("@genereelcoord" eq "@{$pointeur_gene_pred}"){
	$VPg++;
      }
    }
  }
}

############################################################################ 
## Fonction d'evaluation au niveau exon
sub evaluation_niveau_exon{
  my@E=@{@_[0]};  # 1er argument
  my@P=@{@_[1]};  # 2eme argument
  my$i=0;;
  my$j=0;;

  # nbre d'exons reels et predits:
  $nER+=($#E+1)/2;
  $nEP+=($#P+1)/2;

  # Vrais positifs:
  for($i=1;$i<=$#E;$i+=2){ # toutes les fins d'exons reels
    for($j=1;$j<=$#P;$j+=2){ # toutes les fins d'exons predits
      # exon correctement predit:
      $VPe+=(($E[$i-1]==$P[$j-1])&&($E[$i]==$P[$j]));
    }
  }
}

############################################################################ 
## Fonction determinant si la position i est predite dans un etat codant
sub ispredcoding{
  my$pos=@_[0]; # 1er arg: position
  my@P=@{@_[1]};# 2eme arg: tab. des frontieres d'exons predits
  my$i=0;
  for($i=0;$i<$#P;$i+=2) {
    if ($P[$i]>0) {
      if ( ($P[$i] <= $pos) && ($P[$i+1] >= $pos) ) {
	return 1;
      }
    }
    else {
      if ( ($P[$i+1] <= (-$pos)) && ($P[$i] >= (-$pos)) ) {
	return 1;
      }
    }
  }
  return 0;
}

############################################################################ 
## Fonction determinant s'il y a un chgmt d'etat predit dans un intervalle
sub ispredchanged{
  my$beg=@_[0]; # 1er arg
  my$end=@_[1];
  my@P=@{@_[1]};
  my$i=0;
  foreach $i (@P) {
    if ( (abs($i) >= $beg) || (abs($i) <= $end) ) {
      return 1;
    }
  }
  return 0;
}

############################################################################ 
## Fonction d'evaluation au niveau nucleotide
sub evaluation_niveau_nucleotide{
  my@E=@{@_[0]};  # 1er argument
  my@P=@{@_[1]};  # 2eme argument
  my$i=0;;
  my$j=0;;
  my$precedent=$borneg; # cas de l'offset gauche

#  print "evaluation niveau nt gau=$borneg dro=$borned\n";
  for( $i=0 ; $i<$#E ; $i+=2 ) { # pour tous les debuts d'exons reels
    for( $j=$precedent ; $j< abs($E[$i]) ; $j++ ){  # chq nt a partir de la derniere frontiere C->NCodant
      $FPn += (&ispredcoding($j,\@P));    # (en amont de l'exon reel, donc dans de l'intron)
      if (&ispredcoding($j,\@P)) { $nNP ++ }
      
#      print"i=$i, pos=$j, Ei= $E[$i], Ei+1= $E[$i+1], vpn= $VPn, fpn= $FPn, fnn= $FNn, nNP= $nNP\n";
    }
    for ($j = abs($E[$i]); $j <= abs($E[$i+1]); $j++) { # chq nt de l'exon reel
      $nNR ++;   # nbre total de nucleotides codants
      if (&ispredcoding($j,\@P)) {      # predit codant aussi, OK
	$VPn++;
	$nNP ++ ;
#      print"i=$i, pos=$j, Ei= $E[$i], Ei+1= $E[$i+1], vpn= $VPn, fpn= $FPn, fnn= $FNn, nNP= $nNP\n";
      }
      else { # loupe!
#	print" non codant\n";
	$FNn++;
#      print"i=$i, pos=$j, Ei= $E[$i], Ei+1= $E[$i+1], vpn= $VPn, fpn= $FPn, fnn= $FNn, nNP= $nNP\n";
      }
    }
    $precedent= abs($E[$i+1])+1;
  }
  for ($j=$precedent; $j <= $borned; $j++) { # cas de l'offset droit
    $FPn += (&ispredcoding($j,\@P));
    if (&ispredcoding($j,\@P)) { $nNP ++ }
#      print"i=$i, pos=$j, Ei= $E[$i], Ei+1= $E[$i+1], vpn= $VPn, fpn= $FPn, fnn= $FNn, nNP= $nNP\n";
  }
#  for($i=0;$i<$#P;$i+=2) { # pour tous les debuts d'exons predits ## ATTENTION!! cas de l'offset chevauchant!
#    $nNP += (abs($P[$i+1]) - abs($P[$i])) +1;   # nbr total de nt predits codants
#      print"i=$i, pos=$j, Ei= $E[$i], Ei+1= $E[$i+1], vpn= $VPn, fpn= $FPn, fnn= $FNn, nNP= $nNP\n";
#  }
}

############################################################################ 
# Renvoie le nombre d'exons identiques entre 2 genes
sub nombre_exons_identiques{
  my@gene1=@{@_[0]};
  my@gene2=@{@_[1]};
  my$i=0;
  my$j=0;
  my$neID=0;
  for($i=1;$i<=$#gene1;$i+=2){ # toutes les fins d'exons gene1
    for($j=1;$j<=$#gene2;$j+=2){ # toutes les fins d'exons gene2
      # exon identique:
      $neID+=(($gene1[$i-1]==$gene2[$j-1])&&($gene1[$i]==$gene2[$j]));
    }
  }
  return $neID;
}

############################################################################ 
# Fonction nouvelle d'evaluation tous niveaux
# et d'affichage des resultats
sub evaluation_et_affichage{
  my@COORD=@{@_[0]}; # 1er argument
  my@PRED=@{@_[1]};  # 2eme argument
  my@E=@{@_[2]};
  my@P=@{@_[3]};
  my$seq=@_[4];

  # nbre de genes reels/predits
  my$ngrtot=($#COORD + 1);
  my$ngptot=($#PRED + 1);
  $nGR+=$ngrtot;
  $nGP+=$ngptot;

  # nbre d'exons reels et predits:
  my$nertot=($#E+1)/2;
  my$neptot=($#P+1)/2;
  $nEP+=$neptot;
  $nER+=$nertot;

  my$i=0;
  my$j=0;
  my$c="";
  my$p="";
  my$npredparfaites=0;
  my$nexonstrouves=0;
  print("\n -SEQUENCE ANALYSEE: $seq\n");

  my$ngr=0;
  my$ngp=0;
  @PRINT=();

  # pour chaque gene reel de la sequence
  foreach $c (@COORD){
    $ngr++;
    my@genereel = @{$c};
    my$ne=($#genereel+1)/2;
    # affichage coordonnees reelles du gene:
    print("gene reel  $ngr/$ngrtot:");
    for($i=0;$i<$#genereel;$i+=2){
      print(" $genereel[$i]-$genereel[$i+1]");
    }
    print(' ... ');
    $predit=0;
    $ngp=0;
    foreach $p (@PRED){
      # on cherche parmi les genes predits
      $ngp++;
      if("@genereel" eq "@{$p}"){
	# Vrais positifs genes:
	$VPg++;
	$npredparfaites++;
	$predit=2;
	print("OK! (detecte en pred n°$ngp)");
      }
      else{
	$neid=&nombre_exons_identiques(\@genereel,$p);
	if ($neid>0){
	  print("rate! ($neid exons trouve(s) sur $ne en pred n°$ngp)");
	  $predit=1;
	}
      }
    }
    ($predit<2) ? push(@PRINT,"fng") : push(@PRINT,"vpg");
    ($predit>=1)? print("\n") : print("tous les exons rates!\n");
  }

  # affichage coordonnees des genes predits:
  $ngp=0;
  foreach $p (@PRED){
    $ngp++;
    my@genepred=@{$p};
    my$ne=($#genepred+1)/2;
    print("prediction $ngp/$ngptot:");
    for($i=0;$i<$#genepred;$i+=2){
      print(" $genepred[$i]-$genepred[$i+1]");
    }
    print(" ... ");
    $predok=0;
    $numerogene=0;
    foreach $p (@COORD){
      $numerogene++;
      if("@genepred" eq "@{$p}"){
	$predok=2;
	print("OK! (parfaite pour gene n°$numerogene)");
      }
      else{
	$neid=&nombre_exons_identiques(\@genepred,$p);
	if ($neid>0){
	  print("imparfaite (cf gene n°$numerogene)");
	  $predok=1;
	}
      }
    }
    if ($predok<2) { push(@PRINT,"fpg") } # vpg deja mis
    ($predok==0)? print("faux positif!\n") : print("\n");
  }

  # EXONS:
  $tmpe=&nombre_exons_identiques(\@E,\@P);
  $VPe+= $tmpe;
  $nexonstrouves+= $tmpe;

  for($i=1;$i<=$#P;$i+=2){ # toutes les fins d'exons P
    my$flag=0;
    for($j=1;$j<=$#E;$j+=2){ # toutes les fins d'exons E
      # exon identique:
      $flag+=(($E[$i-1]==$P[$j-1])&&($E[$i]==$P[$j]));
    }
    ($flag==0) ? push(@PRINT,"fpe") : push(@PRINT,"vpe");
  }
  for($i=1;$i<=$#E;$i+=2){ # toutes les fins d'exons E
    my$flag=0;
    for($j=1;$j<=$#P;$j+=2){ # toutes les fins d'exons P
      # exon identique:
      $flag+=(($E[$i-1]==$P[$j-1])&&($E[$i]==$P[$j]));
    }
    if ($flag==0) { push(@PRINT,"fne") }
  }

  print("Code  pour $seq: @PRINT\n");

  print("Total pour $seq: ");
#  print("Pour $ngrtot genes reels, $ngptot prediction(s), dont $npredparfaites parfaite(s) \n");
#  print("Pour $nertot exons reels, $neptot prediction(s), dont $nexonstrouves parfaite(s)\n\n");
#  print("$npredparfaites genes bien predits sur $ngrtot (avec $ngptot predictions)\n");
#  print("$nexonstrouves exons bien predits sur $nertot (avec $neptot predictions)\n");
  print("$ngrtot gene(s), $ngptot predit(s), $npredparfaites trouve(s), $nertot exon(s), $neptot predit(s), $nexonstrouves trouve(s)\n");
}

########    "LEGENDE" ou guide des symboles ###################
##
##  FP/FN=Faux positifs/negatifs,
##  VP/VN=Vrais positifs/negatifs,
##  suivi de g=gene, e=exon, ou n=nucleotide.
##
##  tableau E =coordonnees des exons reels pour le gene courant
##  tableau P =coordonnees des exons predits pour le gene courant
##
##  nG/nE/nN = nombre de genes, exons et nucleotides,
##  suivi de R=reellement codants, ou P=predits comme codants.
##  
################################################################

$FPg=0;$FNg=0;$VNg=0;$VPg=0;$nGR=0;$nGP=0;
$FPe=0;$FNe=0;$VNe=0;$VPe=0;$nER=0;$nEP=0;
$FPn=0;$FNn=0;$VNn=0;$VPn=0;$nNR=0;$nNP=0;$nNtot=0;
$totalseqlen=0;$AC=0;$CC=0;
@coord=();@pred=();
@E=();
@P=();
@COORD2D=();
@PRED2D=();
$nseqreelles=1;
$nseqpred=1;

######################################################################
####################          ARGMT            #######################
######################################################################
# lecture des arguments (pas joli, mais pas besoin de use Getopt::Std)
if ($#ARGV != 1) {
  if ($ARGV[0] eq '-h') { description() }
  if ($#ARGV < 1) { usage() }
  for (my$i=2 ; $i<=$#ARGV ; $i++) {
    if ($ARGV[$i] =~ /^\-(.)(.+)/) {
      if ($1 eq 'o') {
	if (!($2 =~ /^(\d+)$/)) { usage() }
	$setoffset=1;
	$offset= $1; 
      }
      elsif ($1 eq 'f') {
	$name2="1";
	open(LS,"$2") || die "Can't open seqID list $2\n";
      }
      elsif ($1 eq 'p') {
	if (($2 eq 's') || ($2 eq 'l')) {
	  $sortie= ( ($2 eq 's') ? 0 : 2 );
	}
	else { usage() }
      }
      elsif ($1 eq 'n') {
	if ($2 eq 't'){
	  $ntlevel=1;
	}
	elsif ($2 eq 'tplus') {
	  if ($name2) {
	    $ntlevel=2;
	  }
	  else {
	    die "special ntplus option requires first the list of fasta sequences (-f)\n";
	  }
	}
	else { usage() }
      }
      else { usage() }
    }
    else { usage() }
  }
}

open(COORD,"$ARGV[0]") || die "Can't open fich coord reelles $ARGV[0]\n";
open(PRED,"$ARGV[1]") || die "Can't open fich coord predites $ARGV[1]\n";

$nlCOORD=0;$nlPRED=0;
$numgene=0;$j=0;
$newseqreelle=1;

if ($sortie!=0) {
  print("\nEVALUATION DES PREDICTIONS D'EUGENE (@CMD)\n");
}

while($lCOORD=<COORD>) {
  chomp $lCOORD;
  $nlCOORD++;
  if (($name2) && ($newseqreelle==1)) {
    chomp($name2=<LS>);
    $newseqreelle=0;
  }
  # Si on est dans une sequence
  if($lCOORD =~ /(\s+\-?([0-9]+))+/) {
    $numgene++;
    @TMP=split(/\s+/,$lCOORD);
#TMP!! on pourrait verifier si vide avant de shifter (+ de souplesse de format coord)
    shift @TMP;  # case vide
    foreach $coord (@TMP){
      push (@E,$coord); #stocke ensemble des exons de la sequence
    }
    @COORD2D[$numgene-1]=([@TMP]);   #stocke ensemble des genes de la sequence
  }
  else{
    if($lCOORD == ""){ # Fin de la sequence
      $newseqreelle=1;

      # calcul des bornes droite et gauche
      $borneg= ( (abs($COORD2D[0][0]) - $offset) > 0 ) ? abs($COORD2D[0][0]) - $offset : 1 ;
      $borned= abs($COORD2D[$#COORD2D][$#{$#COORD2D}]) + $offset;

      ## LECTURE des predictions d'EuGene
      $j=0;
      $gene_en_cours=0;
      foreach $pointeur2D (@PRED2D){
	@{$pointeur2D}=(); }
      @PRED2D=();
      @P=();
      @TMP=();

      chomp($lPRED=<PRED>);$nlPRED++;
      if($lPRED eq ""){ # pas de gene trouve
	$j++;
	$name= ( ($name2) ? $name2 : " -NO_NAME- " );
      }

      while($lPRED ne ""){ # tant qu'il y a un exon trouve
	if ($lPRED =~ /^([^\s]+)\s+([a-zA-Z]+)\s+([\+\-])\s+([0-9]+)\s+([0-9]+)\s+[0-9]+\s+[^\s]+\s+[^\s]+\s+([0-9]+)\s+([0-9]+)\s+/) {
	  if ( ($setoffset==0) ||
	       ((abs($5) >= $borneg) && (abs($4) <= $borned))) {
	    $name= ( ($name2) ? $name2 : $1);
	    $type=$2;
	    $sens=$3;
	    $beg=$4;
	    $end=$5;
#	    ($name= $1) =~ s/\.\d+\.\d+\.\d+//;
	    $gene_en_cours=1;
	    # Si format ou fichier mauvais:
	    if (!( ($type eq "Init")||($type eq "Intr")||($type eq "Term")||($type eq "Sngl"))){ die" Pb(1) format fichier predictions $ARGV[1] ligne $nlPRED\n"};

	    # Stockage des exons
	    $signe=( ("$sens" eq '+')? "" : '-');
	    push(@P,"$signe$beg");
	    push(@P,"$signe$end");
	    push(@TMP,"$signe$beg");
	    push(@TMP,"$signe$end");

	    ## fin d'un gene
	    if ( (($type eq "Term")&&("$sens" eq '+'))||(($type eq "Init")&&("$sens" eq '-'))||($type eq "Sngl")){
	      $j++;
	      $gene_en_cours=0;
	      @PRED2D[$j-1]=[@TMP];
	      @TMP=();
	    }
	  }
	}
	elsif($lPRED !~ /\s+Utr\d\s/) {die"Pb(2) format fichier predictions $ARGV[1] ligne $nlPRED\n"};
	chomp($lPRED=<PRED>);
	$nlPRED++;
      }
      if ($gene_en_cours==1) {
	# Cas ou la derniere prediction "continue" apres la sequence
	# (gene non termine par un "+Term" ou "-Init")
	$j++;
	@PRED2D[$j-1]=[@TMP];
	@TMP=();
	$gene_en_cours=0;
      }
      # Les tableaux de coordonees et predictions sont prets
      # Si l'offset n'a pas ete defini (ttes les pred ont ete prises), on prend les valeurs extremes (pour niveau nt)
      if ($setoffset==0) {
	$borneg= (abs($PRED2D[0][0]) < $borneg) ? abs($PRED2D[0][0]) : $borneg ;
	$borned= (abs($PRED2D[$#PRED2D][$#{$#PRED2D}]) > $borned) ? abs($PRED2D[$#PRED2D][$#{$#PRED2D}]) : $borned; 
      }
      else { # par contre, si offset est defini, on prend l'intervalle mini
	if (abs($PRED2D[0][0]) < abs($COORD2D[0][0])) {
	  $borneg= (abs($PRED2D[0][0]) < $borneg) ? $borneg : abs($PRED2D[0][0]);
	}
	else { 
	  $borneg= abs($COORD2D[0][0]);
	}
	$borned= (abs($PRED2D[$#PRED2D][$#{$#PRED2D}]) > $borned) ? $borned : abs($PRED2D[$#PRED2D][$#{$#PRED2D}]);
      }

      if ($sortie == 2) {
	if ($ntlevel>0) { &evaluation_niveau_nucleotide(\@E,\@P); }
	&evaluation_et_affichage(\@COORD2D,\@PRED2D,\@E,\@P,$name);
      }
      else {
	&evaluation_niveau_gene(\@COORD2D,\@PRED2D);
	&evaluation_niveau_exon(\@E,\@P);
	if ($ntlevel>0) { &evaluation_niveau_nucleotide(\@E,\@P);}
      }
      $numgene=0;
      if ($ntlevel==2) {
	$totalseqlen += &seqlength($name);
      }

      # liberation memoire tableaux:
      @E=();
      @P=();
      foreach $pointeur2D (@COORD2D){
	@{$pointeur2D}=(); }
      @COORD2D=();
      foreach $pointeur2D (@PRED2D){
	@{$pointeur2D}=(); }
      @PRED2D=();
    }
    else{
      die"Pb(3) format fichier coordonnees $ARGV[0] ligne $nlCOORD\n";
    }
  }
}

if (($nGR==0)||($nER==0)){die"Pb(4) : nbre de seq soumises (ou d'exons) nul!\n"}

$SENSg=($VPg*100)/$nGR;
$SENSe=($VPe*100)/$nER;
if ($ntlevel>0) {$SENSn=($VPn*100)/$nNR;}

if (($nGP==0)||($nEP==0)){
  $SPECg=100;
  $SPECe=100;
  $SPECn=100;
}
else{
  $SPECg=($VPg*100)/$nGP;
  $SPECe=($VPe*100)/$nEP;
  if ($ntlevel>0) { $SPECn=($VPn*100)/$nNP; }
}

if ($ntlevel>0) { 
  if ($SENSn != (100*$VPn)/($VPn+$FNn)) {
    die "internal error computing nt sens:\n VPn=$VPn , nNR=$nNR , sensn= $SENSn\n VPn=$VPn , FNn=$FNn\n";
  }
  if ($SPECn != (100*$VPn)/($VPn+$FPn)) {
    die "internal error computing nt spec:\n$VPn*100/$nNP != 100*($VPn/($VPn+$FPn))\n VPn=$VPn , nNP=$nNP , specn= $SPECn\n VPn=$VPn , FPn=$FPn\n";
  }
}
#print"VPn=$VPn , nNP=$nNP , specn= $SPECn\n VPn=$VPn , FPn=$FPn\n";
if ($sortie==0) {
  print"$SENSg $SPECg $SENSe $SPECe\n";
}
else {
  print "\n>TOTAL (@CMD)\n";
  print "$VPg genes bien detectes sur $nGR avec $nGP predictions\n";
  print "$VPe exons bien detectes sur $nER avec $nEP predictions\n";
  if ($ntlevel>0) {  print "$VPn nt bien detectes sur $nNR avec $nNP predictions\n";}
#print "$VPn nt codants bien detectes sur $nNR ($FNn rates);$nNP predits dont $FPn faux pos\n";
#print "\n";
#print "SENSIBILITE GENES : $SENSg\n";
#print "SPECIFICITE GENES : $SPECg\n";
#print "SENSIBILITE EXONS : $SENSe\n";
#print "SPECIFICITE EXONS : $SPECe\n";
#print "SENSIBILITE NUCL  : $SENSn\n";
#print "SPECIFICITE NUCL  : $SPECn\n";
  print"SNG: $SENSg SPG: $SPECg SNE: $SENSe SPE: $SPECe";
  if ($ntlevel>0) {  print " SNN: $SENSn SPN: $SPECn";}
}
if ($ntlevel == 2) {
  $VNn= $totalseqlen -$FPn -$FNn -$VPn;
  $AC= 0.5*( $VPn/($VPn+$FNn) + $VPn/($VPn+$FPn) + $VNn/($VNn+$FPn) + $VNn/($VNn+$FNn) ) -1;
  $CC= ( ($VPn*$VNn)-($FNn*$FPn) ) / ( sqrt( ($VPn+$FNn)*($VNn+$FPn)*($VPn+$FPn)*($VNn+$FNn) ) );
  print " AC: $AC CC: $CC";
}
print"\n";
if($name2) {close LS}
close COORD;
close PRED;
