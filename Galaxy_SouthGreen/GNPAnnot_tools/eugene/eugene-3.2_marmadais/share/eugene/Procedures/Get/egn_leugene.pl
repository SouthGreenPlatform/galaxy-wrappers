#!/usr/bin/perl -w

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
# $Id: egn_leugene.pl,v 1.1 2005/03/17 13:41:23 cros Exp $
# ------------------------------------------------------------------
# File:     egn_leugene.pl
# Contents: if needed, lauch local NetStart, NetGene2, SplicePredictor 
# and eugene on sequences files
# ------------------------------------------------------------------

use HTTP::Request::Common qw(POST);
use LWP::UserAgent;

=head1 Description

 This program takes as parameter a nucleic acid sequence file in FASTA
 format or a folder containing such files (*.tfa *.fasta). Extra parameters
 are consider as eugene additional parameters.
 For each sequence, launch local NetStart, NetGene2, SplicePredictor
 and if only the EUGENEDIR environment is defined (ending with '/'), 
 lauch eugene with HTML output.

 Beware, if you want eugene to take into account EST or Blast
 results, the files have to be present (correctly named).
 Note that used paths correspond to the sources organization (therefor
 need to be updated fot others configurations).

=cut

#------------------------------------------------------------
sub usage( $ )
  {
    printf STDERR "%s\n";
    system("pod2text $0");
    exit(-1);
  }

#------------------------------------------------------------
sub getsites( $ )
  {
    local($sequence, $rsequence, $l);
    local($flag);
    local($i, $n, $r, $file,$DoIt,$Olap,$StepLength);
    local($MaxLength);

    $sequence = &fasta2flat($_[0]);
    $l = length($sequence);
    $rsequence = &reverse(&complement($sequence));


    # NetStart ###############################################
    $Olap = 100;
    $MaxLength = 20000;
    $StepLength =  $MaxLength - 2*$Olap;

    # NetStart on direct strand 
    $DoIt = 0;
    printf STDERR "NetStart [2*%d request(s)]:",($l-2*$Olap+$StepLength-1)/$StepLength;
    $Last = 0;

    if (!( -e "$_[0].starts") || (-s "$_[0].starts" == 0))    {$DoIt = 1; }

    if ($DoIt)  { open(STARTS,">$_[0].starts") || die "Can't create $_[0].starts file\n";}
    else { printf STDERR "[F on disk] ";}

    # split input file into MaxLength kb parts because of netstart limitation
    if ($DoIt) {
      for ($i = 0; $i + 2*$Olap < $l; $i += $StepLength) {

	# create temp file with the sequence portion
	&flat2fasta("/Annotation/SunSoft/netstart-1.0/tmp/$$.tfa", 'tmp', substr($sequence,$i,$MaxLength), 60);
	
	printf STDERR "Calling netstart on [%d - %d] / $l of $_[0]\n", $i, $i + $MaxLength;
	open(FIN,"cd /Annotation/SunSoft/netstart-1.0; ./netstart -at ./tmp/$$.tfa |") || die "Can't launch NetStart on $[0]";
	$flag = 0;
	while (<FIN>)
	  {
	    chomp;
	    if ($flag == 1 && $_ eq "") { $flag = 0;}
	    if ($flag == 1) {
	      if (($n,$r) = (/^\s+(\d+)(.+)$/)) {
		if ((($n >= $Olap) || ($i == 0)) && (($n < $MaxLength-$Olap) || $Last  == 1))
		  { printf STARTS "%7d%s\n", $n + $i, $r; }
	      }
	    }
	    if (/----------------------/) { $flag = 1; }
	  }
	close(FIN);
	unlink("/Annotation/SunSoft/netstart-1.0/tmp/$$.tfa");
      }
      close(STARTS);
    }

    # NetStart on reverse strand
    $DoIt=0;

    if (!( -e "$_[0].startsR") || (-s "$_[0].startsR" == 0))  {$DoIt = 1;}

    if ($DoIt) { open(RSTARTS,">$_[0].startsR") || die "Can't create $_[0].startsR file\n";}
    else { printf STDERR "[R on disk] ";}

    if ($DoIt) {
      # split input file into 50 kb parts because of netstart limitation
      for ($i = 0; $i + 2*$Olap < $l; $i += $StepLength) {

	# create temp file with the sequence portion
	&flat2fasta("/Annotation/SunSoft/netstart-1.0/tmp/$$.tfa", 'tmp', substr($rsequence,$i,$MaxLength), 60);
	printf STDERR "Calling netstart on [%d - %d] / $l of $_[0] (rev)\n", $i, $i + $MaxLength;
	open(FIN,"cd /Annotation/SunSoft/netstart-1.0; ./netstart -at ./tmp/$$.tfa |") || die "Can't launch NetStart on $_[0]";
	$flag = 0;
	while (<FIN>)
	  {
	    chomp;
	    if ($flag == 1 && $_ eq "") { $flag = 0;}
	    if ($flag == 1) {
	      if (($n,$r) = (/^\s+(\d+)(.+)$/)) {
		if ((($n >= $Olap) || ($i == 0)) && (($n < $MaxLength-$Olap) || $Last  == 1))
		  { printf RSTARTS "%7d%s\n", $n + $i, $r; }
	      }
	    }
	    if (/----------------------/) { $flag = 1; }
	  }
	close(FIN);
	unlink("/Annotation/SunSoft/netstart-1.0/tmp/$$.tfa");
      }
      close(RSTARTS);
    }

    # NetGene2 ##################################################################

    $DoIt = 0;
    if (!( -e "$_[0].splices") || (-s "$_[0].splices" == 0))    {$DoIt = 1; }
    if (!( -e "$_[0].splicesR") || (-s "$_[0].splicesR" == 0))    {$DoIt = 1; }

    printf STDERR "\nNetGene2: ";

    if ($DoIt) {
      &flat2fasta("$$.tfa", "$$", $sequence, 60);
      system("netgene2 -a -e -p -r $$ -s at $$.tfa");
      ($? == 0) || die "fail to run netgene2";
      (-e "$$.score.txt" && -e "$$.comp.score.txt") || die "no result file";
      rename("$$.score.txt", "$_[0].splices");
      rename("$$.comp.score.txt", "$_[0].splicesR");
      unlink("$$.results");
      unlink("$$.tfa");
    }
    else {print STDERR "[on disk]";}

    # Splice Predictor ######################################################
    # direct strand
    $DoIt = 0;
    if (!( -e "$_[0].spliceP") || (-s "$_[0].spliceP" == 0))    {$DoIt = 1;}
    printf STDERR "\nSplicePredictor: ";

    if ($DoIt) {
      print STDERR "Issuing request to SplicePredictor (direct strand)...\n";
      open(FORWARD,">$_[0].spliceP");
      open(FIN,"SplicePredictorLL  -c 1 -s Arabidopsis -p 5 -l $_[0] |") || die "Can't start Splice Predictor";
      while(<FIN>)
	{
	print FORWARD;
	}
      close(FIN);
      close(FORWARD);
    }
    else {print STDERR "[F on disk]";}

    # reverse strand
    $DoIt = 0;
    if (!( -e "$_[0].splicePR") || (-s "$_[0].splicePR" == 0))    {$DoIt = 1; }

    if ($DoIt) {
      print STDERR "Issuing request to SplicePredictor (reverse strand)...\n";
      open(REVERSE,">$_[0].splicePR");
      open(FIN,"SplicePredictorLL -r -c 1 -s Arabidopsis -p 5 -l $_[0] |") || die "Can't start Splice Predictor";
      while(<FIN>)
	{
	  print REVERSE;
	}
      close(FIN);
      close(REVERSE);
    }
    else {print STDERR "[R on disk]";}
  }

#------------------------------------------------------------
# reverse a flat sequence
sub reverse ( $ )
  {
    my ($seq) = @_;

    return join('', reverse (split '', $seq));
  }

#------------------------------------------------------------
# complement a flat sequence
sub complement( $ )
  {
    my ($seq) = @_;

    $seq = uc($seq);
    $seq =~ tr/('A','T','G','C')/('T','A','C','G')/;
    return $seq;
  }

#------------------------------------------------------------
# read a FASTA file and return the first flat sequence if no AC is given
# or the one that have the given AC (empty if none match)
sub fasta2flat
  {
    my($fastaF, $AC);
    my ($ac, $seq, $inseq) = ('', '', 0);

    if (scalar @_ > 0)
      {
	$fastaF = $_[0];
	if (scalar @_ == 1) {$AC = '';} 
	else {$AC = $_[1];}
      }
    else { die "usage: fasta2flat file_name [AC]";}

    open(FF,$fastaF) || die "Can't open $fastaF";
    while (<FF>)
      {
	$_ =~ s/[\n\r\f]//g;

	if (/^>/)
	  {
	    if ($inseq == 1) {last;} # end of the sequence
	    ($ac) = (/^>([^\s]+)/); # retrieve current accession number
	    if ($AC eq '' || $AC eq $ac) # it's the sequence we are looking for
	      {
		$inseq = 1;
		next;
	      }
	  }

	if ($inseq == 1)
	  {
	    $seq .= $_;
	    next;
	  }
      }
    close(FF);
    return $seq;
  }

#------------------------------------------------------------
# Append a sequence to a FASTA file
sub flat2fasta ( $ $ $ $ )
  {
    my($fname, $AC, $seq, $line_length) = @_;
    my($i);

    open(FF,">>$fname") || die "can't write into $fname";
    print FF ">$AC\n";
    for ($i=0; $i < length($seq); $i += $line_length)
      {
	printf FF "%s\n", substr($seq,$i,$line_length);
      }
    close(FF);
  }

#------------------------------------------------------------
# Call preprocessing for a sequence and launch eugene
sub eugene( $ @ )
  {
    local($file, @params) = @_;
    local($cmd);

    print STDERR "\nprocessing $file\n\n";
    &getsites($file);
    if (exists($ENV{'EUGENEDIR'}))
      {
	$cmd = sprintf("$ENV{'EUGENEDIR'}eugene -ph %s -g $file > $file.html 2> $file.trace", join(' ',@params));
	system($cmd);
      }
  }

#------------------------------------------------------------
# Main procedure
($#ARGV != -1) || &usage('usage:');

if (!exists($ENV{'EUGENEDIR'}))
  {
    system('echo "########################################################"');
    system('echo "WARNING: Only preprocessing of sequence(s) will be realized,"');
    system('echo "eugene will not be launched because"');
    system('echo "the EUGENEDIR environment variable is not set."');
    system('echo "########################################################"');
  }

system('echo "Started on "`date`');

$param1 = shift @ARGV;

if (-f "$param1") # just launch  on one file
  {
    &eugene($param1, @ARGV);
  }
elsif (-d "$param1") # launch on all fasta files in this folder
  {
    $param1 =~ s/\/$//;
    @L = glob("$param1/*.tfa");
    push @L, glob("$param1/*.fasta");
    (scalar @L > 0) || die "$param1 contains no fasta files (*.tfa *.fasta)";
    foreach $f (@L)
      {
	&eugene($f, @ARGV);
      }
  }
else
  {
    &usage("$param1 is neither a file, nor a folder");
  }

system('echo ""');system('echo ""');
system('echo "Finished on "`date`');
if (!exists($ENV{'EUGENEDIR'}))
  {
    system('echo "########################################################"');
    system('echo "WARNING: Only preprocessing of sequence(s) was realized,"');
    system('echo "eugene was not launched because"');
    system('echo "the EUGENEDIR environment variable is not set."');
    system('echo "########################################################"');
  }
