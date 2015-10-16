#!/usr/bin/perl
# tout petit script perso pour parser les sorties de
# blastall -p tblastx
# et les rendre EuGene-compatibles
# accepte aussi les 
# blastall -p blastn !!
# et les blastx
# filtre les hits "cognates" (meme nom et memes coordonnees du hit)

($#ARGV==0) || die "Usage : <nom d'un fichier sortie blastall -p tblastx ou blastn a parser>\n";

$blastfile=$ARGV[0];
open(BLASTFILE,"$blastfile") || die "pb ouverture fichier blast $blastfile!\n";

$programme=<BLASTFILE>;
if ($programme =~ /BLASTN/) {$programme="blastn";}
elsif ($programme =~ /^BLASTX/) {$programme="blastx";}
else {$programme="tblastx";}

$nbrehits=0;
$querydeb= $queryfin = $sbjctdeb = $sbjctfin = $querylgr = $hitlgr = 0;
$score= $bits = $E = $queryframe = $sbjctframe = 0;
$idnt = $idpc =0; #identite (nbre de nt et pourcentage)
$brin= "";
$SEQ="";

while ($l=<BLASTFILE>){
  if ($l =~ /^Query=\s([^\n\s]+)\s?/) {
    $nomquery = $1;
    $l=<BLASTFILE>;
    if ($l=~/\((\d+)\sletters\)/){
      $querylgr=$1;
    }
  }
  elsif ($l =~ /^Database:\s([a-zA-Z0-9]+)\s?/) {
    $nomdb = $1;
  }
  elsif ($l =~ /^>([^\n\s]+)/){
    $nomsbjctsuivant=$1;
  }
  elsif ($l =~ /Identities\s=\s(\d+)\/(\d+)\s\((\d+)/){
    $idnt=$1;
    $idpc=$3;
    $hitlgr=$2;
  }
  # CHAQUE HIT (et cas du premier et dernier hit)
  elsif ( ($l =~ /^\sScore\s+=\s+(.+)\s+bits\s\(([0-9]+)\)\,\sExpect.*\s=\s([^\n]+)/) ||
	  ($l =~ /^\s\sDatabase:\s/)) {
#    $nbrehits= ( ($nbrehits==0) ? 1 : (print"$querydeb $queryfin $score $E $queryframe $nomsbjct".'_'."$sbjctdeb".'_'."$sbjctfin $sbjctdeb $sbjctfin\n")  ) ;
    if ($nbrehits==0){
      $nbrehits++;
    }
    else {
      if ( ("$nomsbjct" ne "$nomquery") ||
	   ($querydeb != $sbjctdeb)     ||
	   ($queryfin != $sbjctfin) ) {
	## (filtre du cognate)
	if ($programme eq "tblastx"){
	  print "$querydeb $queryfin $bits $E $queryframe $nomsbjct".'_'."$sbjctdeb".'_'."$sbjctfin $sbjctdeb $sbjctfin $SEQ\n";
	}
	elsif ($programme eq "blastn") {
	  for ($i=-3;$i<=3;$i++){
	    $brin = ( ($i>0) ? "+" : "" ) ;
	    if ($i!=0){
	      print "$querydeb $queryfin $bits $E $brin$i $nomsbjct".'_'."$sbjctdeb".'_'."$sbjctfin".'_'."$i $sbjctdeb $sbjctfin\n"
	    }
	  }
	}
	elsif ($programme eq "blastx"){
	  print "$querydeb $queryfin $bits $E $queryframe $nomsbjct $sbjctdeb $sbjctfin\n";
	}
      }
    }
    $nomsbjct=$nomsbjctsuivant;
    $score=$1;
    $bits= $2;
    $E=$3;
    $querydeb= 0;
    $queryfin =0;
    $sbjctdeb= 0;
    $sbjctfin =0;
    $queryframe = $sbjctframe = 0;
    $SEQ="";
  }
  elsif ($l =~ /^\sFrame\s=\s([^\s]+)\s\/\s(.+)$/){ # tblastx
    $queryframe = $1;
    $sbjctframe = $2;
  }
  elsif ($l =~ /^\sFrame\s=\s([^\s]+)$/){ # blastx
    $queryframe = $1;
  }
# brin non informatif pour blastn
#  elsif ($l =~ /^\sStrand\s=\s([a-zA-Z]+)/){
#    $brin= ( ("$1" eq "Plus") ? "+"  : "-" )  ;
#  }
  elsif ($l =~ /^Query:\s([0-9]+).+\s([0-9]+)/){
    $querydeb= ( ($querydeb==0) ? $1 : $querydeb );
    $queryfin= $2;
  }
  elsif ($l =~ /^Sbjct:\s([0-9]+)\s+(.+)\s+([0-9]+)/){
    $sbjctdeb= ( ($sbjctdeb==0) ? $1 : $sbjctdeb );
    $SEQ= "$SEQ"."$2";
    $sbjctfin= $3;
  }
}
close BLASTFILE;
