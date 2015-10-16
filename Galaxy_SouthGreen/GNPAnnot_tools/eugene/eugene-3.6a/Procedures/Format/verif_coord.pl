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
# $Id: verif_coord.pl,v 1.2 2004-09-16 11:56:08 cros Exp $
# ------------------------------------------------------------------
# File:     verif_coord.pl
# Contents: Just a little script for coordinates verification by
#           starts stops and splices extraction from the new genomic
#           sequences. Need a list of fasta files (the genomic seqs)
#           and a coord file (cf. sim2eugene output)
# ------------------------------------------------------------------


use strict;
#use warnings;
use IO::Handle;

my $embosspath = "";
my $cmd_seqret = $embosspath."seqret";
my $cmd_extseq = $embosspath."extractseq";

my $outputF = 'oryclean.verif';
unlink $outputF;

my $begin  = 0;
my $end    = 0;
my $before = 0;   # context nucleotide to extract before and after
my $after  = 0;   # the "GT/AG" consensus
my $options;

($#ARGV==1) || die "Usage: $0 coordfile fastalistfile\n";

my $coordfile = $ARGV[0];
my $list      = $ARGV[1];

open(COORD,"$coordfile") || die "can't open coord file $coordfile\n";
open(SEQ,  "$list")      || die "can't open list file $list\n";

while(my $seq = <SEQ>) {
  chomp $seq;
  print "\n-sequence $seq\n";
  
  while(my $lcoord = <COORD>) {
    if ($lcoord =~ /(\s+\-?\d+)+/) {
      my $negsense = ($lcoord =~ /\-\d+\s/);
      $lcoord =~ s/\-//g;
      my @coord=split(/\s+/,$lcoord);
      shift @coord;               # empty
      
      open(OUTPUT, ">>$outputF")||die "ERROR: Could not open $outputF";
      print OUTPUT "=----------------------------------";
      print OUTPUT "---------------------------=\n";
      print OUTPUT "$seq"."\n";
      
      ###########  START VERIFICATION  ###############
      if ($negsense) {
	print OUTPUT "   Start "."$coord[$#coord]"."\t";
	close OUTPUT;
	$begin= $coord[$#coord]-2;
	$end  = $coord[$#coord];
	$options = "-sbegin $begin -send $end -sreverse -osformat raw";
	system "$cmd_seqret $options -stdout -auto $seq >> $outputF"; 
      }
      else {
	print OUTPUT "   Start "."$coord[0]\t";
	close OUTPUT;
	$begin= $coord[0];
	$end  = $coord[0]+2;
	$options = "-regions $begin-$end -osformat raw";
	system "$cmd_extseq $seq $options -stdout -auto >> $outputF";
      }
      
      if ($lcoord =~ /\s+\-?\d+\s+\-?\d+(\s+\-?\d+)+/) { 
	############  SPLICES SITES ###########
	if ($negsense==0) {
	  for (my $i=1; $i<$#coord-1; $i+=2 ) {
	    # DON
	    open(OUTPUT,">>$outputF") || die "ERROR: Could not open $outputF";
	    print OUTPUT "   Don   ".($coord[$i])."\t";
	    close OUTPUT;
	    $begin= $coord[$i] - ($before-1);
	    $end  = $coord[$i] + ($after +2);
	    $options = "-regions $begin-$end -osformat raw";
	    system "$cmd_extseq $seq $options -stdout -auto >> $outputF";
	    # ACC
	    open(OUTPUT,">>$outputF") || die "ERROR: Could not open $outputF";
	    print OUTPUT "   Acc   ".($coord[$i+1])."\t";
	    close OUTPUT;
	    $begin= $coord[$i+1] - ($before+2);
	    $end  = $coord[$i+1] + ($after -1);
	    $options = "-regions $begin-$end -osformat raw";
	    system "$cmd_extseq $seq $options -stdout -auto >> $outputF";
	  }
	}
	else {
	  for (my $i=$#coord-1; $i>1; $i-=2 ) {
	    # DON
	    open(OUTPUT,">>$outputF") || die "ERROR: Could not open $outputF";
	    print OUTPUT "   Don   ".($coord[$i])."\t";
	    close OUTPUT;
	    $begin= $coord[$i] - ($before+2);
	    $end  = $coord[$i] + ($after -1);
	    $options = "-sbegin $begin -send $end -sreverse -osformat raw";
	    system "$cmd_seqret $options -stdout -auto $seq >> $outputF";
	    # ACC
	    open(OUTPUT,">>$outputF") || die "ERROR: Could not open $outputF";
	    print OUTPUT "   Acc   ".($coord[$i-1])."\t";
	    close OUTPUT;
	    $begin= $coord[$i-1] - ($before-1);
	    $end  = $coord[$i-1] + ($after +2);
	    $options = "-sbegin $begin -send $end -sreverse -osformat raw";
	    system "$cmd_seqret $options -stdout -auto $seq >> $outputF";
	  }
	}
      }
      
      ###########  STOP VERIFICATION  ###############
      open(OUTPUT, ">>$outputF") || die "ERROR: Could not open $outputF";
      if ($negsense) {
	print OUTPUT "   Stop  "."$coord[0]\t";
	close OUTPUT;
	$begin= $coord[0];
	$end  = $coord[0]+2;
	$options = "-sbegin $begin -send $end -sreverse -osformat raw";
	system "$cmd_seqret $options -stdout -auto $seq >> $outputF"; 
      }
      else {
	print OUTPUT "   Stop  "."$coord[$#coord]"."\t";
	close OUTPUT;
	$begin= $coord[$#coord]-2;
	$end  = $coord[$#coord];
	$options = "-regions $begin-$end -osformat raw";
	system "$cmd_extseq $seq $options -stdout -auto >> $outputF";
      }
    }
    elsif ($lcoord == "") {
      last;
    }
  }
}

open(OUTPUT, ">>$outputF") || die "ERROR:";
print OUTPUT "=----------------------------------";
print OUTPUT "---------------------------=\n";
close OUTPUT;
