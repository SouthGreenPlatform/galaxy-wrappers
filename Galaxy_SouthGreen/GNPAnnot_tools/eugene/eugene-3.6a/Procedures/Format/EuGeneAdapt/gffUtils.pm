package gffUtils;

use strict;
use warnings;
use GFF;

#=============================================================================#
#=             Copyright (c) 2003 by INRA. All rights reserved.              =#
#=                 Redistribution is not permitted without                   =#
#=                 the express written permission of INRA.                   =#
#=                     Mail : tschiex@toulouse.inra.fr                       =#
#=---------------------------------------------------------------------------=#
#= File         : gffUtils.pl                                                =#
#= Description  : Package for eugeneAdapt.pl                                 =#
#= Authors      : P.Bardou, S.Foissac, M.J.Cros, A.Moisan, T.Schiex          =#
#= History      : version 1.0 (30 oct. 2003)				     =#
#= Improvements :                                                            =#
#=============================================================================#


my $gff;


##---------------------------------------------------------------------------##
sub gffGetCoord {
  my $inputGff = shift;
  my @geneCoord;

  $gff = new GFF::GeneFeatureSet;
  gffRead   ($inputGff);
  gffVerify ();

  my @lineList = $gff->order_gf();
  my $strand = $lineList[0]->strand();
  if ($strand eq "+") { $strand = ""; }
  my $i = 1;

  if ($lineList[0]->feature() eq "UTR5" ||
      $lineList[0]->feature() eq "UTR3") {
    push(@geneCoord, $strand.$lineList[0]->start());
  }
  else { push(@geneCoord, 0); $i--; }

  for (; $i<$gff->count()-1; $i++) {
    if ($lineList[$i]->feature() ne "Intron") {
      push(@geneCoord, $strand.$lineList[$i]->start());
      push(@geneCoord, $strand.$lineList[$i]->end());
    }
  }

  if ($lineList[-1]->feature() eq "UTR5" ||
      $lineList[-1]->feature() eq "UTR3") {
    push(@geneCoord, $strand.$lineList[-1]->end());
  }
  else {
    push(@geneCoord, $strand.$lineList[-1]->start());
    push(@geneCoord, $strand.$lineList[-1]->end());
    push(@geneCoord, 0);
  }

  return @geneCoord;
}

##---------------------------------------------------------------------------##
sub gffRead {
  my $inputGff = shift;

  open(GFFHANDLE,$inputGff)||die "ERROR: Could not open GFF file $inputGff\n";
  $gff->read(\*GFFHANDLE);
  close(GFFHANDLE);
}

##---------------------------------------------------------------------------##
sub gffWrite {
  my ($id, $outGffFile, $outInfoFile, $fastaFile, @tab) = @_;
  my $st;                   # pour le champ strand
  my $u1;                   # utr1
  my $u2;                   # utr2
  my $nbExon = (@tab-2)/2;
  my $ne1;                  # nameExon1
  my $ne2;                  # nameExon2
  my $rank;                 # pour le champ rank (initialisation)

  open(TFA,$fastaFile)||die"ERROR: Could not open new fasta file $fastaFile\n";
  my @fasta = <TFA>;
  my $seq   = "";
  close TFA;
  shift(@fasta);
  my $seqLen = 0;
  foreach my $line (@fasta) {
    chomp $line;
    $seqLen += length($line);
    $seq    .= $line;
  }

  # Calcul des effectifs
  my @nuc = computeNuc($seq, @tab);

  # Calcul des frames
  my @frame = computeFrame($seqLen, @tab);
  my $cmpFrame = 0;
  my $cmpNuc   = 8;

  open(GFFHDLE, ">$outGffFile") ||
    die"ERROR: Could not open GFF file $outGffFile\n";
  open(INFOHDLE,">$outInfoFile")||
    die"ERROR: Could not open Inof file $outInfoFile\n";
  if ($tab[1] > 0)  {
    $u1 = "UTR5";
    $u2 = "UTR3";
    $st = "+";
    $rank = 1;
    if ($nbExon == 1) { $ne1 = "E.Sngl"; }
    else              { $ne1 = "E.Init"; $ne2 = "E.Term"; }
  }
  else {
    $u1 = "UTR3";
    $u2 = "UTR5";
    $st = "-";
    $rank = $nbExon;
    if ($nbExon == 1) { $ne1 = "E.Sngl"; }
    else              { $ne1 = "E.Term"; $ne2 = "E.Init"; }

    for (my $i=0; $i<@tab; $i++) {
      $tab[$i] = $tab[$i]*-1;
    }
  }

  my $nbTab = "";
  if (length($id) > 7)  { $nbTab =  "\t"; }
  if (length($id) > 14) { $nbTab .= "\t"; }
  print GFFHDLE "#S.Name$nbTab\tSource\tFeature\tStart\tEnd\t".
    "Score\tStrand\tFrame\tA\tT\tC\tG\tOthers\tGC%\tLength\t%3\tRank\n";

  if ($tab[0] != 0) {
    print GFFHDLE "$id\teugAdap\t$u1\t".abs($tab[0])."\t".($tab[1]-1).
      "\t.\t$st\t.\t#$nuc[0]\t$nuc[1]\t$nuc[2]\t$nuc[3]\t$nuc[4]\t$nuc[5]".
	"\t$nuc[6]\t$nuc[7]\tU0\n";
  }

  print GFFHDLE "$id\teugAdap\t$ne1\t".$tab[1]."\t".$tab[2].
    "\t.\t$st\t".$frame[$cmpFrame]."\t#";
  for (my $c=0; $c<8; $c++) { print GFFHDLE $nuc[$cmpNuc++]."\t"; }
  print GFFHDLE "E$rank";

  for (my $i=3; $i+1<$#tab-1; $i+=2) {
    print GFFHDLE "\n$id\teugAdap\tIntron\t".($tab[$i-1]+1)."\t".($tab[$i]-1).
      "\t.\t$st\t.\t#";
    for (my $c=0; $c<8; $c++) { print GFFHDLE $nuc[$cmpNuc++]."\t"; }
    if ($st eq "-") { $rank--; }
    print GFFHDLE "I$rank";
    print GFFHDLE "\n$id\teugAdap\tE.Intr\t".$tab[$i]."\t".$tab[$i+1].
      "\t.\t$st\t".$frame[++$cmpFrame]."\t#";
    for (my $c=0; $c<8; $c++) { print GFFHDLE $nuc[$cmpNuc++]."\t"; }
    if ($st eq "+") { $rank++; }
    print GFFHDLE "E$rank";
  }

  if ($nbExon != 1) {
    print GFFHDLE "\n$id\teugAdap\tIntron\t".($tab[-4]+1)."\t".($tab[-3]-1).
      "\t.\t$st\t.\t#";
    for (my $c=0; $c<8; $c++) { print GFFHDLE $nuc[$cmpNuc++]."\t"; }
    if ($st eq "-") { $rank--; }
    print GFFHDLE "I$rank";
    print GFFHDLE "\n$id\teugAdap\t$ne2\t".$tab[-3]."\t".$tab[-2].
      "\t.\t$st\t".$frame[++$cmpFrame]."\t#";
    for (my $c=0; $c<8; $c++) { print GFFHDLE $nuc[$cmpNuc++]."\t"; }
    if ($st eq "+") { print GFFHDLE "E".++$rank; }
    else            { print GFFHDLE "E$rank";    }
  }

  if ($tab[-1] != 0) {
    print GFFHDLE "\n$id\teugAdap\t$u2\t".($tab[-2]+1)."\t".abs($tab[-1]).
      "\t.\t$st\t.\t#";
    for (my $c=0; $c<8; $c++) { print GFFHDLE $nuc[$cmpNuc++]."\t"; }
    print GFFHDLE "U0";
  }
  else { $cmpNuc += 8; }

  print INFOHDLE "ID\t\t$id";
  print INFOHDLE "\nNbExon  \t".$nbExon;
  print INFOHDLE "\nLenExon \t".$nuc[$cmpNuc++];
  print INFOHDLE "\nGCCDS   \t".$nuc[$cmpNuc++];
  print INFOHDLE "\nGCExonMin\t".$nuc[$cmpNuc++];
  print INFOHDLE "\nGCExonMax\t".$nuc[$cmpNuc++];
  print INFOHDLE "\nNbIntron\t".($nbExon-1);
  print INFOHDLE "\nLenIntron\t".$nuc[$cmpNuc++];
  print INFOHDLE "\nGCIntron\t".$nuc[$cmpNuc++];
  print INFOHDLE "\nGCUTR   \t".$nuc[$cmpNuc++];
  print INFOHDLE "\nLenGene \t".$nuc[$cmpNuc++];
  print INFOHDLE "\nLencDNA \t".$nuc[$cmpNuc++];
  close(GFFHDLE);
  close(INFOHDLE);
}

##---------------------------------------------------------------------------##
sub gffVerify {
  my $nbUtr5 = 0;
  my $nbUtr3 = 0;
  my $nbInit = 0;
  my $nbIntr = 0;
  my $nbTerm = 0;
  my $nbSngl = 0;

  my $feature_filter =
    sub
      {
	my $self = shift;
	my $s = $self->feature();
	if ($s eq "UTR5"   || $s eq "UTR3" ||
	    $s eq "E.Init" || $s eq "E.Intr" ||
	    $s eq "E.Term" || $s eq "E.Sngl" ||
	    $s eq "Intron") {
	  if ($s eq "UTR5")   { $nbUtr5++; }
	  if ($s eq "UTR3")   { $nbUtr3++; }
	  if ($s eq "E.Init") { $nbInit++; }
	  if ($s eq "E.Intr") { $nbIntr++; }
	  if ($s eq "E.Sngl") { $nbSngl++; }
	  if ($s eq "E.Term") { $nbTerm++; }
	  return 1;
	}
	else {
	  print "\nERROR in GFF file.\n".
	    "WARNING :\n".
	      "   - Feature must be UTR5, UTR3, E.Init,".
		" E.Intr, E.Term, E.Sngl or Intron.\n".
		  "   - Fields must be tab-separated\n";
	  exit;
	}
      };

  my $tempgff = new GFF::GeneFeatureSet;
  $tempgff = $gff->filter($feature_filter);

  if($nbSngl == 0 && $nbInit==0 && $nbTerm == 0 && $nbIntr == 0) {
    print "\nERROR in GFF file. No exons.\n";
    exit;
  }
  if(($nbUtr3 >1 || $nbUtr5 >1 || $nbSngl >1 || $nbInit >1 || $nbTerm >1) ||
     ($nbSngl == 1 && ($nbInit == 1 || $nbTerm == 1))){
    print "\nERROR in GFF file.".
      " This file must contain only ONE complete gene.\n";
    exit;
  }
  if(($nbInit == 1 && $nbTerm == 0) ||
     ($nbInit == 0 && $nbTerm == 1) ||
     ($nbInit == 0 && $nbTerm == 0 && $nbIntr > 0)) {
    print "\nERROR in GFF file. Only one COMPLETE gene in GFF file.\n";
    exit;
  }
}

##---------------------------------------------------------------------------##
sub computeFrame {
  my ($seqLen, @tab) = @_;

  my @frame  = ();
  my $sens  = ( ($tab[1] > 0) ? "" : "-" );

  if ($sens eq "") {
    push(@frame, ($tab[1]-1) % 3);
    for (my $i=3; $i<$#tab; $i+=2) {
      my $intronFrame = ($tab[$i] - 1 - $tab[$i-1]) % 3;
      push(@frame, (($frame[-1] + $intronFrame) % 3));
    }
  }
  else {
    push(@frame, ($seqLen-abs($tab[-2])) % 3);
    for (my $i=$#tab-2; $i>1; $i-=2) {
      my $intronFrame = (abs($tab[$i]) - 1 - abs($tab[$i-1])) % 3;
      unshift(@frame, (($frame[0] + $intronFrame) % 3));
    }
  }
  return @frame;
}

##---------------------------------------------------------------------------##
sub computeNuc {
  my ($seq, @tab) = @_;

  my $seqLen = length($seq);
  my $seqTmp;
  my @nuc = ();
  my %c;
  my $gc;
  my $nbGC;
  my $nbGCExs  = 0;  # nombre de GC dans les Exons
  my $nbGCInts = 0;  # nombre de GC dans les Introns
  my $nbGCUTRs = 0;  # nombre de GC dans les UTRs
  my $lenUTRs  = 0;  # longueur des UTRs
  my $lenExs   = 0;  # longueur de tous les exons
  my $exsGC;         # GC% de tous les exons
  my $exMinGC = 100; # GC% min des exons
  my $exMaxGC = 0;   # GC% max des exons
  my $lenInts = 0;   # longueur de tous les introns
  my $intsGC;        # GC% de tous les introns

  for (my $i=0; $i<$#tab; $i++) {
    my $begin;
    my $len;
    $c{"a"} = $c{"A"} = $c{"c"} = $c{"C"} = 0;
    $c{"g"} = $c{"G"} = $c{"t"} = $c{"T"} = 0;

    if ($i == 0)
      { $begin = abs($tab[$i]) - 1; }
    else
      { $begin = abs($tab[$i]) - 1 + (($i+1)%2); }

    if ($i == 0  ||  $i == $#tab-1)
      { $len = abs($tab[$i+1]) - abs($tab[$i]) + 1 - (($i+1)%2); }
    else
      { $len = abs($tab[$i+1]) - abs($tab[$i]) + 1 - (($i+1)%2*2); }

    $seqTmp = substr($seq, $begin, $len) ;
    $seqTmp =~s /[^aAtTcCgG]//g;

    my @l = split(//,$seqTmp);
    foreach my $l (@l) {
      $c{$l} += 1;
    }

    push(@nuc, $c{"a"} + $c{"A"});
    push(@nuc, $c{"t"} + $c{"T"});
    push(@nuc, $c{"c"} + $c{"C"});
    push(@nuc, $c{"g"} + $c{"G"});
    push(@nuc, $len - ($nuc[-1] + $nuc[-2] + $nuc[-3] + $nuc[-4]));
    $nbGC = $nuc[-2] + $nuc[-3];
    $gc   = sprintf ("%.2f", $nbGC / $len * 100);
    push(@nuc, $gc);
    push(@nuc, $len);
    push(@nuc, $len%3);

    if ($i != 0 && $i != $#tab-1) {
      if ($i%2 != 0) {
	$lenExs  += $len;
	$nbGCExs += $nbGC;
	if ($gc < $exMinGC) { $exMinGC = $gc; }
	if ($gc > $exMaxGC) { $exMaxGC = $gc; }
      }
      if ($i%2 == 0) {
	$lenInts  += $len;
	$nbGCInts += $nbGC;
      }
    }
    elsif ($nbGC != 0) {
      $lenUTRs  += $len;
      $nbGCUTRs += $nbGC;
    }
  }

  push(@nuc, $lenExs);
  push(@nuc, sprintf ("%.2f", ($nbGCExs / $lenExs * 100)));
  push(@nuc, $exMinGC);
  push(@nuc, $exMaxGC);
  push(@nuc, $lenInts);
  if ($lenInts != 0)
    {
      push(@nuc, sprintf ("%.2f", $nbGCInts / $lenInts * 100)); }
  else
    {
      push(@nuc, 0);
    }
  if ($tab[0] + $tab[-1] != 0)
    {
      push(@nuc, sprintf ("%.2f", $nbGCUTRs / $lenUTRs * 100));
    }
  else
    {
      push(@nuc, 0);
    }
  push(@nuc, $lenExs + $lenInts + $lenUTRs);
  push(@nuc, $lenExs + $lenUTRs);

  return @nuc;
}

1;
