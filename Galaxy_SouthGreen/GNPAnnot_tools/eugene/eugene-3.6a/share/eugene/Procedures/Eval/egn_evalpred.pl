#!/usr/bin/perl -I../lib

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
# $Id: egn_evalpred.pl,v 1.4 2010-01-25 08:18:20 sallet Exp $
# ------------------------------------------------------------------
# File:     egn_evalpred.pl
# Contents: evaluation of a gene annotation
# ------------------------------------------------------------------

BEGIN
{
    my ($dir,$file) = $0 =~ /(.+)\/(.+)/;
    unshift (@INC,"$dir/../lib");
}


use ParamParser;
use strict;



# Variables globales
my $NOMduPROG = $0;
my @CMD       = "$NOMduPROG "."@ARGV";

# Nombre de genes reels et predits
my $nGR = 0;
my $nGP = 0;
# genes FP/FN/VN/VP
my $FPg = 0;
my $FNg = 0;
my $VNg = 0;
my $VPg = 0;
# Nombre d'exons reels et predits
my $nER = 0;my $nEP = 0;
# exons VP
my $VPe = 0;
# FP/VP/FN au niveau nucleotides
my $FPn = 0;
my $VPn = 0; 
my $FNn = 0;

# Nombre de nucleotides predits et reels
my $nNP=0; my $nNR=0;

sub description
{
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
sub usage
{
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
sub seqlength
{
  my ($filename) = @_;  

  my $len = 0;
  open(FASTASEQ,"$filename") || die "Can't open seq file $filename in seqlength\n";
  # Check first character of the file is >
  if (!( <FASTASEQ> =~ /^\s?>/ ) ) 
  { 
    die "$filename has to be a fasta file to allow length sequence computing.\n";
  }

  while ( <FASTASEQ> )
  {
    chomp;
    $len += length($_);
  }

  return $len;
}


=head2 procedure evaluation_niveau_gene

 Title        : evaluation_niveau_gene
 Usage        : -
 Function     : Fonction d'evaluation au niveau gene qui pour une sequence donnee, 
                regarde si chaque gene reel a une prediction parfaitement egale
                Met a jour le nombre de gene reel et predit ($nGR, $n,GP) et le nombre
                de gene Vrai positif ($VPg)
 Args         : $ra_coord2d: reference sur un tableau de genes reels
                $ra_pred2d : reference sur un tableau de genes predits
 Globals      : none

=cut
sub evaluation_niveau_gene
{
  my ($ra_coord2d, $ra_pred2d) = @_;

  # nbre de genes reels/predits
  my $nGeneReels    = scalar(@$ra_coord2d);
  my $nGenesPredits = scalar(@$ra_pred2d);
  $nGR += $nGeneReels;
  $nGP += $nGenesPredits;
  
  # Si la structure des deux genes est identique, incrementer le nb de VP
  foreach my $pointeur_gene_reel (@$ra_coord2d)
  {
    my @genereelcoord = @{$pointeur_gene_reel};
    foreach my $pointeur_gene_pred (@$ra_pred2d)
    {
	$VPg++ if("@genereelcoord" eq "@{$pointeur_gene_pred}");
    }
  }

  return;
}


=head2 procedure evaluation_niveau_exon

 Title        : evaluation_niveau_exon
 Usage        : -
 Function     : Fonction d'evaluation au niveau exon qui pour une sequence donnee, 
                regarde le nombre d exon correctement predit
                Met a jour le nombre d'exons reel et predit ($nER, $nEP) et le nombre
                d'exon Vrai positif ($VPe)
 Args         : $ra_e: reference sur un tableau de position d exon reel
                $ra_p : reference sur un tableau de positions d exon predit
 Globals      : none

=cut
sub evaluation_niveau_exon
{
  my ($ra_e, $ra_p) = @_;

  # nbre d'exons reels et predits:
  my $nExonsReels   = scalar(@$ra_e);
  my $nExonsPredits = scalar(@$ra_p);
  $nER += $nExonsReels/2;
  $nEP += $nExonsPredits/2;

  # Incrementation du nombre de vrai positif
  # toutes les fins d'exons reels
  for(my $i = 1; $i < $nExonsReels; $i += 2 )
  { 
    # toutes les fins d'exons predits
    for(my $j = 1; $j < $nExonsPredits; $j += 2 )
    { 
      # exon correctement predit:
      $VPe += (($ra_e->[$i-1] == $ra_p->[$j-1])&&
               ($ra_e->[$i]   == $ra_p->[$j]  ) );
    }
  }

  return;
}

=head2 function ispredcoding

 Title        : ispredcoding
 Usage        : -
 Function     : determine si la position pos est predite dans un etat codant
 Args         : $pos  : $position
                $ra_p : reference sur un tableau des positions des exons predits
 Return       : 1 si la position est predit comme codante
 Globals      : none

=cut
sub ispredcoding
{
  my ($pos, $ra_p) = @_;

  my $p_size = scalar(@$ra_p); # Taille du tableau d'exons
  
  for (my $i = 0; $i < ($p_size-1); $i+=2) 
  {
    # Cas ou pos est inclu dans un exon 
    if ($ra_p->[$i] > 0) 
    {
	return 1 if ( ($ra_p->[$i] <= $pos) && ($ra_p->[$i+1] >= $pos) );
    }
    else 
    {
	return 1 if ( ($ra_p->[$i+1] <= (-$pos)) && ($ra_p->[$i] >= (-$pos)) );
    }
  }

  return 0;
}


=head2 procedure evaluation_niveau_nucleotide

 Title        : evaluation_niveau_nucleotide
 Usage        : -
 Prerequisite : -
 Function     : Evaluation tous niveaux et affichage des resultats
 Args         : $ra_e    : reference sur un tableau des positions des exons reels
                $ra_p    : reference sur un tableau des positions des exons predits
                $borneg  : position 5' de la region a evaluer
                $borned  : position 3' de la region a evaluer
 Globals      : none

=cut
sub evaluation_niveau_nucleotide
{
  my ($ra_e, $ra_p, $borneg, $borned) = @_;
 
  my $precedent = $borneg; # cas de l'offset gauche
  my $nb_exons  = scalar(@$ra_e);

  # pour tous les debuts d'exons reels
  for( my  $i = 0; $i < $nb_exons-1; $i+=2 ) 
  {
    # chq nt a partir de la derniere frontiere C->NCodant
    for( my $j = $precedent; $j < abs($ra_e->[$i]); $j++ )
    {
      # En amont de l'exon reel, donc dans de l'intron
      $FPn += (&ispredcoding  ($j, $ra_p) );
      $nNP++ if (&ispredcoding($j, $ra_p));
    }

    # Pour chaque nt de l'exon reel
    for (my $j = abs($ra_e->[$i]); $j <= abs($ra_e->[$i+1]); $j++) 
    {
      $nNR++;   # nbre total de nucleotides codants
      if ( &ispredcoding($j, $ra_p) ) 
      {    
	# predit codant
	$VPn++;
	$nNP++;
      }
      else 
      { # loupe!
	$FNn++;
      }
    }
    $precedent = abs($ra_e->[$i+1])+1;
  }

  # cas de l'offset droit
  for (my $j = $precedent; $j <= $borned; $j++) 
  {
    $FPn += (   &ispredcoding($j,$ra_p) );
    $nNP++ if ( &ispredcoding($j,$ra_p) );
  }

#  for($i=0;$i<$#P;$i+=2) { # pour tous les debuts d'exons predits ## ATTENTION!! cas de l'offset chevauchant!
#    $nNP += (abs($P[$i+1]) - abs($P[$i])) +1;   # nbr total de nt predits codants
#      print"i=$i, pos=$j, Ei= $E[$i], Ei+1= $E[$i+1], vpn= $VPn, fpn= $FPn, fnn= $FNn, nNP= $nNP\n";
#  }

  return;
}

############################################################################ 
# Renvoie le nombre d'exons identiques entre 2 genes
=head2 function evaluation_niveau_nucleotide

 Title        : nombre_exons_identiques
 Usage        : -
 Prerequisite : -
 Function     : Renvoie le nombre d'exons identiques entre 2 genes
 Args         : $ra_gene1  : reference sur un tableau de position d exon du gene1 
                $ra_gene2  : reference sur un tableau de position d exon du gene2
 Return       : nombre d exons identique entre les deux genes
 Globals      : none

=cut
sub nombre_exons_identiques
{
  my ($ra_gene1, $ra_gene2) = @_; 

  my $nb_gene1 = scalar(@$ra_gene1);
  my $nb_gene2 = scalar(@$ra_gene2);
  
  # Nb d'exons identiques
  my $neID = 0;

  # toutes les fins d'exons gene1
  for(my $i=1; $i < $nb_gene1; $i+=2)
  { # toutes les fins d'exons gene2
      for(my $j=1; $j < $nb_gene2; $j+=2)
      {
	  # nombre d exon identique
	  $neID += 1 if ( ($ra_gene1->[$i-1] == $ra_gene2->[$j-1]) && 
			  ($ra_gene1->[$i]   == $ra_gene2->[$j]  ) );
      }
  } 

  return $neID;
}

=head2 procedure evaluation_et_affichage

 Title        : evaluation_et_affichage
 Usage        : -
 Prerequisite : -
 Function     : Evaluation tous niveaux et affichage des resultats
 Args         : $ra_coord: reference sur un tableau de genes reels
                $ra_pred : reference sur un tableau de genes predits
                $ra_e    : reference sur un tableau des positions des exons reels
                $ra_p    : reference sur un tableau des positions des exons predits
                $seq     : nom de la sequence 
 Globals      : none

=cut

sub evaluation_et_affichage
{
  my ($ra_coord, $ra_pred, $ra_e, $ra_p, $seq) = @_;

  # met a jour le nombre de gene reels et predits
  my $ngrtot = scalar(@$ra_coord);
  my $ngptot = scalar(@$ra_pred);
  $nGR += $ngrtot;
  $nGP += $ngptot;
  #  met a jour le nombre d'exons reels et predits
  my $nertot = scalar(@$ra_e)/2;
  my $neptot = scalar(@$ra_p)/2;
  $nEP += $neptot;
  $nER += $nertot;

  my$npredparfaites = 0;
  my$nexonstrouves  = 0;

  print("\n -SEQUENCE ANALYSEE: $seq\n");

  my$ngr    = 0; # numero du gene reel etudie
  my$ngp    = 0; # numero du gene predit etudie
  my $neid;      # Var utilisee pour connaitre le nb d exon identique entre 2 genes    
  my @PRINT = ();

  # Affichage des informations sur chaque gene reel
  foreach my $c (@$ra_coord)
  {
    $ngr++;
    my @genereel = @{$c};
    my $ne        = ($#genereel+1)/2; # Nombre d'exons du gene
 
    print("gene reel  $ngr/$ngrtot:"); 
    # affichage des positions des exons reels du gene
    for ( my $i=0; $i < $#genereel; $i+=2 )
    {
      print(" $genereel[$i]-$genereel[$i+1]");
    }
    print(' ... ');

    my $predit = 0;
    $ngp       = 0; # Numero de la prediction etudiee
    # Recherche si le gene a bien ete predit 
    foreach my $p (@$ra_pred)
    {
      $ngp++; 
      # Cas ou on trouve exactement le gene
      if ("@genereel" eq "@{$p}")
      {
	$VPg++;
	$npredparfaites++;
	$predit = 2;
	print("OK! (detecte en pred n $ngp)");
      }
      else
      {
	$neid = &nombre_exons_identiques(\@genereel, $p);
	# Cas ou la structure exon n'a pas ete totalement trouvee
	if ($neid > 0)
	{
	  print("rate! ($neid exons trouve(s) sur $ne en pred n $ngp)");
	  $predit = 1;
	}
      }
    }

    ($predit <  2) ? push(@PRINT,"fng") : push(@PRINT,"vpg");
    ($predit >= 1) ? print("\n") : print("tous les exons rates!\n");
  }

  # Affichage des informations sur chaque gene predit
  $ngp = 0;
  foreach my $p (@$ra_pred)
  {
    $ngp++;
    my @genepred = @{$p};
    my $ne       = ($#genepred+1)/2;  # Nombre d'exons de la prediction
    print("prediction $ngp/$ngptot:");
    # affichage des positions des exons predits
    for ( my $i=0; $i < $#genepred; $i+=2 )
    {
      print(" $genepred[$i]-$genepred[$i+1]");
    }
    print(" ... ");

    my $predok     = 0;
    my $numerogene = 0;
    # Recherche si la prediction est correcte
    foreach my $c (@$ra_coord)
    {
      $numerogene++;
      # Cas ou la prediction du gene est correct
      if("@genepred" eq "@{$c}")
      {
	$predok = 2;
	print("OK! (parfaite pour gene n $numerogene)");
      }
      else
      {
	$neid = &nombre_exons_identiques(\@genepred, $c);
        # Cas ou la structure exon n'est pas completement correcte
	if ( $neid > 0 )
	{
	  print("imparfaite (cf gene n $numerogene)");
	  $predok=1;
	}
      }
    }
    if ( $predok < 2 ) { push(@PRINT,"fpg") } # vpg deja mis
    ( $predok == 0 )? print("faux positif!\n") : print("\n");
  }

  # Nombre global d exons identiques
  my $tmpe        = &nombre_exons_identiques($ra_e, $ra_p);
  $VPe           += $tmpe;
  $nexonstrouves += $tmpe;

  # 2 boucles pour construire le code de la sequence 
  # (cad pour chaque exon dire si vp, fp ou fn)

  # toutes les fins d'exons predits
  for (my $i = 1; $i < scalar(@$ra_p); $i+=2)
  {
    my $flag = 0;
    # pour toutes les fins d'exons reels
    for (my $j = 1; $j < scalar(@$ra_e); $j += 2 )
    {
      # exon identique:
      $flag += 1 if ( ( $ra_e->[$i-1] == $ra_p->[$j-1] ) && 
		      ( $ra_e->[$i]   == $ra_p->[$j]   ) );
    }
    ( $flag == 0 ) ? push(@PRINT, "fpe") : push(@PRINT, "vpe"); 
  }

  # toutes les fins d'exons reels
  for (my $i = 1; $i < scalar(@$ra_e); $i+=2 )
  { 
    my $flag = 0; 
    # pour toutes les fins d'exons predits
    for(my $j = 1; $j < scalar(@$ra_p); $j+=2 )
    { 
      # exon identique:
      $flag += 1 if ( ($ra_e->[$i-1] == $ra_p->[$j-1]) && 
		      ($ra_e->[$i]   == $ra_p->[$j]  ) );
    }
    if ( $flag == 0) { push(@PRINT,"fne") }
  }

  print("Code  pour $seq: @PRINT\n");
  print("Total pour $seq: ");
  print("$ngrtot gene(s), $ngptot predit(s), $npredparfaites trouve(s), $nertot exon(s), $neptot predit(s), $nexonstrouves trouve(s)\n");

  return;
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



MAIN: 
{

my $setoffset = 0;
my $offset    = 0;
my $sortie    = 1; # 0-> short, 2-> detaillee
my $name2     = "";
my $ntlevel   = 0;

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
	  $sortie = ( ($2 eq 's') ? 0 : 2 );
	}
	else { usage() }
      }
      elsif ($1 eq 'n') {
	if ($2 eq 't'){
	  $ntlevel = 1;
	}
	elsif ($2 eq 'tplus') {
	  if ($name2) 
	  {
	    $ntlevel = 2;
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

my $borneg      = 0;
my $borned      = 0;
my $totalseqlen = 0;
my @E       = (); # Tableau des exons
my @P       = ();
my @COORD2D = (); # Contient les coordonnees des exons des genes reels
my @PRED2D  = ();

# Numero de lignes dans les fichiers d'entree
my $nlCOORD = 0; 
my $nlPRED  = 0;

my $numgene      = 0; # Numero du gene courant
my $j            = 0; # Index pour sauvegarder les predictions (pas clair)
my $newseqreelle = 1; # Utiliser pour savoir si on est dans une nouvelle seq
my $gene_en_cours;

my @TMP;

if ($sortie != 0) 
{
  print("\nEVALUATION DES PREDICTIONS D'EUGENE (@CMD)\n");
}

# For each line of the coordinates file 
# MtC00289_GC#MtC00289_GC -1001 -1160 -1265 -1857
while(my $lCOORD = <COORD>) 
{
  chomp $lCOORD;
  $nlCOORD++;

  # Cas ou on lit le nom de la sequence dans un fichier (option -fNomFichier)
  if ( ( $name2 ) && ($newseqreelle == 1 ) ) 
  {
    chomp($name2 = <LS> );
    $newseqreelle = 0;
  }

  # Si on est dans une sequence
  if ($lCOORD =~ /(\s+\-?([0-9]+))+/) 
  {
    $numgene++;
    @TMP = split(/\s+/,$lCOORD);  # [MtC00289_GC#MtC00289_GC, -1001, -1160, -1265, -1857]
    #TMP!! on pourrait verifier si vide avant de shifter (+ de souplesse de format coord)
    shift @TMP;  # case vide
    foreach my $coord (@TMP)
    {
      push (@E, $coord); #stocke ensemble des exons de la sequence
    }
    
    $COORD2D[$numgene-1] = ([@TMP]); #stocke ensemble des genes de la sequence
  }
  else
  {
    if ($lCOORD == "")
    { # Fin de la sequence
      $newseqreelle = 1;

    # calcul des bornes droite & gauche de la region contenant des genes
      #  (pos 5' du 1er exon du 1er gene, pos 3' du dernier exon du dernier gene)
      $borneg = ( (abs($COORD2D[0][0])-$offset) > 0 ) ? abs($COORD2D[0][0])-$offset : 1;
      my $s_tab  = $#COORD2D;
      my $s_cell = $#{$COORD2D[$s_tab]};
      $borned    = abs($COORD2D[$s_tab][$s_cell]) + $offset;

      ## LECTURE des predictions d'EuGene
      $j             = 0;
      $gene_en_cours = 0;
      
      @PRED2D = ();
      @P      = ();
      @TMP    = ();
      my $name;

      my $lPRED = <PRED>; # Lecture de la ligne suivante dans le fichier de prediction
      chomp($lPRED);
      $nlPRED++;

      if($lPRED eq "")
      { # pas de gene trouve
	$j++;
	$name = ( ($name2) ? $name2 : " -NO_NAME- " );
      }

      while($lPRED ne "")
      { 
        # tant qu'il y a un exon trouve
	# Id  Term  -   1001  1288  288  -1  -3    1289    1000     0.0
	if ($lPRED =~ /^([^\s]+)\s+([a-zA-Z]+)\s+([\+\-])\s+([0-9]+)\s+([0-9]+)\s+[0-9]+\s+[^\s]+\s+[^\s]+\s+([0-9]+)\s+([0-9]+)\s+/) 
	{
	  my $beg = $4;  # 1001, start position
	  my $end = $5;  # 1288, stop position

	  if ( ($setoffset == 0) || 
               ((abs($end) >= $borneg) && (abs($beg) <= $borned)) ) 
	  {
	    $name          = ( ($name2) ? $name2 : $1);
	    my $type       = $2; # 'Term', type of the feature
	    my $sens       = $3; # '-', sign
	    $gene_en_cours = 1;

	    # Si format ou fichier mauvais:
	    if (!( ($type eq "Init")||($type eq "Intr")||($type eq "Term")||($type eq "Sngl")))
	    {
		die " Pb(1) format fichier predictions $ARGV[1] ligne $nlPRED type $type\n";
	    }

	    # Stockage des exons predits
	    my $signe = ( ("$sens" eq '+') ? "" : '-');
	    push(@P,   "$signe$beg");
	    push(@P,   "$signe$end");
	    push(@TMP, "$signe$beg");
	    push(@TMP, "$signe$end");

	    ## fin d'un gene
	    if ( (($type eq "Term")&&("$sens" eq '+')) ||
                 (($type eq "Init")&&("$sens" eq '-')) ||
                  ($type eq "Sngl"))
	    {
	      $j++;
	      $gene_en_cours = 0;
	      @PRED2D[$j-1]  = [@TMP]; # Save the exon prediction
	      @TMP           = ();
	    }
	  }
	}

	elsif($lPRED !~ /\s+Utr\d\s/ && $lPRED !~ /\sncRNA\s/) 
	{
	    die "Pb(2) format fichier predictions $ARGV[1] ligne $nlPRED\n";
	}

	chomp($lPRED = <PRED>); # read next prediction line
	$nlPRED++;
      }


      if ($gene_en_cours == 1) 
      {
	# Cas ou la derniere prediction "continue" apres la sequence
	# (gene non termine par un "+Term" ou "-Init")
	$j++;
	@PRED2D[$j-1]  = [@TMP];
	@TMP           = ();
	$gene_en_cours = 0;
      }

      $s_tab  = $#COORD2D;
      $s_cell = $#{$COORD2D[$s_tab]};

      # Les tableaux de coordonnees et predictions sont prets
      # Si l'offset n'a pas ete defini (ttes les pred ont ete prises), on prend les valeurs extremes (pour niveau nt)
      if ( $setoffset == 0 ) 
      {
	  #print "\nOOOO: Offset 0\n";
	$borneg = (abs($PRED2D[0][0]) < $borneg) ? abs($PRED2D[0][0]) : $borneg ;
	$borned = (abs($PRED2D[$s_tab][$s_cell]) > $borned) ? abs($PRED2D[$s_tab][$s_cell]) : $borned; 
      }
      else 
      { 
	  #print "\nOOOO: Offset 1\n";
        
        # par contre, si offset est defini, on prend l'intervalle mini
	if (abs($PRED2D[0][0]) < abs($COORD2D[0][0])) 
	{
	  $borneg = (abs($PRED2D[0][0]) < $borneg) ? $borneg : abs($PRED2D[0][0]);
	}
	else 
	{ 
	  $borneg = abs($COORD2D[0][0]);
	}
	$s_tab  = $#COORD2D;
	$s_cell = $#{$COORD2D[$s_tab]};
	$borned = (abs($PRED2D[$s_tab][$s_cell]) > $borned) ? $borned : abs($PRED2D[$s_tab][$s_cell]);
      }

      
      if ($sortie == 2) 
      {
	if ($ntlevel > 0) 
	{ 
	    &evaluation_niveau_nucleotide(\@E,\@P, $borneg, $borned); 
	}
	&evaluation_et_affichage(\@COORD2D, \@PRED2D, \@E, \@P, $name);
      }
      else 
      {
	&evaluation_niveau_gene(\@COORD2D,\@PRED2D);
	&evaluation_niveau_exon(\@E,\@P);
	if ( $ntlevel > 0 ) 
	{ 
	  &evaluation_niveau_nucleotide(\@E,\@P, $borneg, $borned);
	}
      }

      $numgene = 0;
      if ($ntlevel == 2) 
      {
	  $totalseqlen += &seqlength($name);
      }

      # liberation memoire tableaux:
      @E       = ();
      @P       = ();
      @COORD2D = ();
      @PRED2D  = ();

    }
    else
    {
      die"Pb(3) format fichier coordonnees $ARGV[0] ligne $nlCOORD\n";
    }
  }
}

if ( ($nGR == 0) || ($nER == 0) )
{
  die "Pb(4) : nbre de seq soumises (ou d'exons) nul!\n";
}

# Calcul de sensibilite
my $SENSg = ($VPg * 100)/$nGR;
my $SENSe = ($VPe * 100)/$nER;
my $SENSn;  
if ( $ntlevel > 0 ) 
{
    $SENSn = ($VPn*100)/$nNR;
}

# Calcul de la specificite
my $SPECg, my $SPECe, my $SPECn;
if ( ($nGP == 0) || ($nEP == 0) )
{
  $SPECg = 100;
  $SPECe = 100;
  $SPECn = 100;
}
else
{
  $SPECg = ($VPg*100)/$nGP;
  $SPECe = ($VPe*100)/$nEP;
  if ($ntlevel > 0) 
  { 
      $SPECn = ($VPn*100)/$nNP; 
  }
}

# Verification de la coherence des calculs de sensibilite et de specificite
if ( $ntlevel > 0 ) 
{ 
  if ($SENSn != (100*$VPn)/($VPn+$FNn)) 
  {
    die "Internal error computing nt sens:\n VPn=$VPn , nNR=$nNR , sensn= $SENSn\n VPn=$VPn , FNn=$FNn\n";
  }
  if ($SPECn != (100*$VPn)/($VPn+$FPn))
  {
    die "internal error computing nt spec:\n$VPn*100/$nNP != 100*($VPn/($VPn+$FPn))\n VPn=$VPn , nNP=$nNP , specn= $SPECn\n VPn=$VPn , FPn=$FPn\n";
  }
}

#print"VPn=$VPn , nNP=$nNP , specn= $SPECn\n VPn=$VPn , FPn=$FPn\n";
if ( $sortie == 0 ) 
{
    print "$SENSg $SPECg $SENSe $SPECe\n";
}
else
{
  print "\n>TOTAL (@CMD)\n";
  print "$VPg genes bien detectes sur $nGR avec $nGP predictions\n";
  print "$VPe exons bien detectes sur $nER avec $nEP predictions\n";
  if ( $ntlevel > 0) 
  {  
      print "$VPn nt bien detectes sur $nNR avec $nNP predictions\n";
  }
  
#print "$VPn nt codants bien detectes sur $nNR ($FNn rates);$nNP predits dont $FPn faux pos\n";
#print "\n";
#print "SENSIBILITE GENES : $SENSg\n";
#print "SPECIFICITE GENES : $SPECg\n";
#print "SENSIBILITE EXONS : $SENSe\n";
#print "SPECIFICITE EXONS : $SPECe\n";
#print "SENSIBILITE NUCL  : $SENSn\n";
#print "SPECIFICITE NUCL  : $SPECn\n";

  print"SNG: $SENSg SPG: $SPECg SNE: $SENSe SPE: $SPECe";
  if ( $ntlevel > 0 )
  {  
      print " SNN: $SENSn SPN: $SPECn";
  }
}

if ($ntlevel == 2) 
{
  my $VNn = $totalseqlen -$FPn -$FNn -$VPn;
  my $AC  = 0.5*( $VPn/($VPn+$FNn) + $VPn/($VPn+$FPn) + $VNn/($VNn+$FPn) + $VNn/($VNn+$FNn) ) -1;
  my $CC  = ( ($VPn*$VNn)-($FNn*$FPn) ) / ( sqrt( ($VPn+$FNn)*($VNn+$FPn)*($VPn+$FPn)*($VNn+$FNn) ) );
  print " AC: $AC CC: $CC";
}
print"\n";

if($name2) {close LS}
close COORD;
close PRED;

}



