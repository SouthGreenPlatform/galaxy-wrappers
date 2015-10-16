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
# $Id: egn_getATGpr4eugene.pl,v 1.1 2005/03/17 13:41:23 cros Exp $
# ------------------------------------------------------------------
# File:     egn_getATGpr4eugene.pl
# Contents: launch of ATGpr on sequence files
# ------------------------------------------------------------------

use HTTP::Request::Common qw(POST);
use HTTP::Request::Common qw(GET);
use LWP::UserAgent;

=head1 Description

 This program takes as parameter a nucleic acid sequence file in FASTA
 format  or  a folder containing  such files   (*.tfa *.fasta), launch
 ATGpr.

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
sub getsites( $ )
  {
    local($file, @params) = @_;
    local($sequence, $rsequence, $seq, $l);
    local($ua, $req, $res, $line,$DoIt,$DoItR);
    local($flag, @content, $zob1, $zob2);
    local($MaxLength,$Last,$LastLength,$Olap);

    print STDERR "\nprocessing $file\n";

    $sequence = &fasta2flat($_[0]);
    #$l = length($sequence);
    $rsequence = &reverse(&complement($sequence));

    $DoIt = 0;
    $DoItR = 0;

    # Request for ATGpr
    print STDERR "ATGpr : ";
    $DoIt = 0;
    if (!( -e "$_[0].atgpr") || (-s "$_[0].atgpr" == 0))    {$DoIt = 1;}
    if (!( -e "$_[0].atgprR") || (-s "$_[0].atgprR" == 0))  {$DoItR = 1;}

    for ($i=0; $i<2; $i++) {
      # $i=0 -> ATGpr FORWARD
      # $i=1 -> ATGpr REVERSE
      if($i==0) {
	print STDERR "F...";
        $Do   = $DoIt;
       	$seq  = $sequence;
	if($Do) {open(FILE,">$_[0].atgpr"); }
      }
      else {
	print STDERR "R...";
	$Do   = $DoItR;
     	$seq  = $rsequence;
	if($Do) { open(FILE,">$_[0].atgprR"); }
      }

      if ($Do) {
	$ua = new LWP::UserAgent;
	$ua->agent("AgentName/0.1 " . $ua->agent);
	
	# Create a request
	$req = POST 'http://www.hri.co.jp/atgpr/cgi-bin/atgpr.cgi',
	  [
	   "number" => 100000000,
	   "seq" => $seq,
	  ];

	# Pass request to the user agent and get a response back
	$res = $ua->request($req);
	
	# Check the outcome of the response
	if ($res->is_success) {
	  if ($i==0) {print STDERR "done  /  "; }
	  else       {print STDERR "done\n";    }
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
	    print FILE $line."\n";
	  }
	  close FILE;
	}
	else { print STDERR "Bad luck this time\n"; }
      }
      else {
	if ($i==0) {print STDERR "[on disk]  /  "; }
	else       {print STDERR "[on disk]\n";    }
      }
    }
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

system('echo "started on "`date`');

$param1 = shift @ARGV;

if (-f "$param1")
  {
    &getsites($param1, @ARGV);
  }
elsif (-d "$param1")
  {
    $param1 =~ s/\/$//;
    @L = glob("$param1/*.tfa");
    push @L, glob("$param1/*.fasta");
    (scalar @L > 0) || die "$param1 contains no fasta files (*.tfa *.fasta)";
    foreach $f (@L)
      {
	&getsites($f, @ARGV);
      }
  }
else
  {
    &usage("$param1 is neither a file, nor a folder");
  }

system('echo "finished on "`date`');
