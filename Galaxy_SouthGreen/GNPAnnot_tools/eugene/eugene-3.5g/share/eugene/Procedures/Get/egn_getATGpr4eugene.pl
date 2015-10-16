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
# $Id: egn_getATGpr4eugene.pl,v 1.2 2005-05-13 15:19:31 bardou Exp $
# ------------------------------------------------------------------
# File:     egn_getATGpr4eugene.pl
# Contents: launch of ATGpr on sequence files
# ------------------------------------------------------------------

use HTTP::Request::Common qw(POST);
use HTTP::Request::Common qw(GET);
use LWP::UserAgent;

=head1 Description

 This program takes as parameter a nucleic acid sequence file in FASTA
 format  or  a folder containing  such files   (*.tfa *.fasta)  If the
 sequence is longer than 20000, it split it into  segment of 20000 bp,
 launch  ATGpr on  each segment, combine the resulting files to obtain
 the same results as if it was  sumitted once.

 Initial version: T. Schiex (INRA, Toulouse)
 Enhancements by: P. Dehais (INRA / Univ. of Gand, Belgium)

=cut


#------------------------------------------------------------

sub usage( $ )
  {
    printf STDERR "%s\n";
    system("pod2text $0");
    exit(-1);
  }

#------------------------------------------------------------
# just a copy of getsites.pl from EuGene dist

sub getsites( $ )
  {
    local($sequence, $rsequence, $seq, $l);
    local($ua, $req, $res, $line,$DoIt,$DoItR);
    local($flag, @content, $zob1, $zob2);
    local($MaxLength,$Last,$LastLength,$Olap);
    local($i, $p, $n, $r, @result, $file);

    $sequence = &fasta2flat($_[0]);
    $l = length($sequence);
    $rsequence = &reverse(&complement($sequence));

    $MaxLength = 20000;
    $Olap = 100;
    $StepLength =  $MaxLength - 2*$Olap;
    $DoIt = 0;
    $DoItR = 0;

    printf STDERR "ATGpr [2*%d request(s)]: ",
      ($l-2*$Olap+$StepLength-1)/$StepLength;
    $Last = 0;

    if (!( -e "$_[0].atgpr")  || (-s "$_[0].atgpr" == 0))   {$DoIt = 1; }
    if (!( -e "$_[0].atgprR") || (-s "$_[0].atgprR" == 0))  {$DoItR = 1;}

    if ($DoIt)  { open(STARTS,">$_[0].atgpr") || die "Can't create $_[0].atgpr file\n";}
    else { printf STDERR "[F on disk] ";}
    if ($DoItR) { open(RSTARTS,">$_[0].atgprR") || die "Can't create $_[0].atgprR file\n";}
    else { printf STDERR "[R on disk] ";}


    for ($i = 0; $i + 2*$Olap < $l; $i += $StepLength) {

      if ($DoIt) {
	# Create a request for ATGpr
	$seq = substr($sequence, $i, $MaxLength);
	
	if ($i+$MaxLength >= $l) {$Last = 1;}
	
	$ua = new LWP::UserAgent;
	$ua->timeout(3600); # 1 hour !!!
	$ua->agent("AgentName/0.1 " . $ua->agent);
	
	$req = POST 'http://www.hri.co.jp/atgpr/cgi-bin/atgpr.cgi',
	  [
	   "number" => 1000000000,
	   "seq" => $seq,
	  ];

	printf STDERR "F%d", 
	  ($i)/$StepLength + 1;
	# Pass request to the user agent and get a response back
	$res = $ua->request($req);
	
	# Check the outcome of the response
	if ($res->is_success) {
	  print STDERR ".";
	  @result = split '\n',$res->content;
	
	  $output = "";
	  foreach $line (@result) {
	    if ($line =~ /<TD/) {
	      if ($line =~ /<.+><.+>(.+)<.+><.+>/) {
		$output .= "$1\t";
	      }
	      else { $output .= "\n"; }
	    }
	  }
	
	  @out = split '\n',$output;
	  @out = sort {(split '\t', $a)[0] <=> (split '\t', $b)[0]} @out;
	  foreach $line (@out) {
	    @spLine = split '\t', $line;
	    if ((($spLine[4] >= $Olap) || ($i == 0)) && (($spLine[4] < $MaxLength-$Olap) || $Last  == 1))
	      {
		$spLine[4] += $i;
		$spLine[5] += $i;
		foreach $field (@spLine) {
		  print STARTS $field."\t";
		}
		print STARTS "\n";
	      }
	  }
	}
	else {
	  print STDERR "Bad luck this time\n";
	}
      }

      if ($DoItR) {
	# Create a request for ATGpr on reverse strand
	$seq = substr($rsequence, $i, $MaxLength);
	
	$ua = new LWP::UserAgent;
	$ua->agent("AgentName/0.1 " . $ua->agent);
	
	$req = POST 'http://www.hri.co.jp/atgpr/cgi-bin/atgpr.cgi',
	  [
	   "number" => 1000000000,
	   "seq" => $seq,
	  ];

	printf STDERR "R%d", 
	  ($i)/$StepLength + 1;
	
	# Pass request to the user agent and get a response back
	$res = $ua->request($req);
	
	# Check the outcome of the response
	if ($res->is_success) {
	  print STDERR ".";
	  @result = split '\n',$res->content;
	  
	  $output = "";
	  foreach $line (@result) {
	    if ($line =~ /<TD/) {
	      if ($line =~ /<.+><.+>(.+)<.+><.+>/) {
		$output .= "$1\t";
	      }
	      else { $output .= "\n"; }
	    }
	  }
	  
	  @out = split '\n',$output;
	  @out = sort {(split '\t', $a)[0] <=> (split '\t', $b)[0]} @out;
	  foreach $line (@out) {
	    @spLine = split '\t', $line;
	    if ((($spLine[4] >= $Olap) || ($i == 0)) && (($spLine[4] < $MaxLength-$Olap) || $Last  == 1))
	      {
		$spLine[4] += $i;
		$spLine[5] += $i;
		foreach $field (@spLine) {
		  print RSTARTS $field."\t";
		}
		print RSTARTS "\n";
	      }
	  }
	}
	else {
	  print STDERR "Bad luck this time\n";
	}
      }
    }
    
    if ($DoIt)  { close STARTS; }
    if ($DoItR) { close(RSTARTS); }
    printf STDERR "done\n";

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

sub eugene( $ @ )
  {
    local($file, @params) = @_;
    local($cmd);

    print STDERR "\nprocessing $file\n\n";
    &getsites($file);
    # JER $cmd = sprintf("./EuGeneAS -p h %s -g $file > $file.html", join(' ',@params));
    # JER system($cmd);
  }
#------------------------------------------------------------
# reas a FASTA file and return the first flat sequence if no AC is given
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

($#ARGV != -1) || &usage('usage:');

if (!exists($ENV{'EUGENEDIR'}))
  {
    $ENV{'EUGENEDIR'} = '/group/biocomp/Soft/src/Eugene/dist';
    $ENV{'PATH'} = "$ENV{'PATH'}:$ENV{'EUGENEDIR'}";
  }

system('echo "started on "`date`');

$param1 = shift @ARGV;

if (-f "$param1") # just launch EuGene on this file
  {
    &eugene($param1, @ARGV);
  }
elsif (-d "$param1") # launch EuGene on all fasta files in this folder
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

system('echo "finished on "`date`');
