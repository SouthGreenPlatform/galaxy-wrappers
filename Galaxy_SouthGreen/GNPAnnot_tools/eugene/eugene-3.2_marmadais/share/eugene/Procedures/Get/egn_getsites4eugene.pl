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
# $Id: egn_getsites4eugene.pl,v 1.1 2005/03/17 13:41:23 cros Exp $
# ------------------------------------------------------------------
# File:     egn_getsites4eugene.pl
# Contents: see description below
# ------------------------------------------------------------------

use HTTP::Request::Common qw(POST);
use HTTP::Request::Common qw(GET);
use LWP::UserAgent;

=head1 Description

 This program takes as parameter a nucleic acid sequence file in FASTA
 format  or  a folder containing  such files   (*.tfa *.fasta)  If the
 sequence is longer than 20000, it split it into  segment of 20000 bp,
 launch  NetStart and NetGene2 on  each segment, combine the resulting
 files to obtain the same results as if it was  sumitted once. Call to
 SplicePredictor is not changed because  it can handle huge sequences.
 Launch then EuGeneAS with HTML output Extra parameter are consider as
 EuGene parameter If  you want  EuGeneAS to take  into account  EST or
 Blast results, the files have to be present (correctly named).

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

    printf STDERR "NetStart [2*%d request(s)]: ",
      ($l-2*$Olap+$StepLength-1)/$StepLength;
    $Last = 0;

    if (!( -e "$_[0].starts") || (-s "$_[0].starts" == 0))    {$DoIt = 1; }
    if (!( -e "$_[0].startsR") || (-s "$_[0].startsR" == 0))  {$DoItR = 1;}

    if ($DoIt)  { open(STARTS,">$_[0].starts") || die "Can't create $_[0].starts file\n";}
    else { printf STDERR "[F on disk] ";}
    if ($DoItR) { open(RSTARTS,">$_[0].startsR") || die "Can't create $_[0].startsR file\n";}
    else { printf STDERR "[R on disk] ";}


    for ($i = 0; $i + 2*$Olap < $l; $i += $StepLength) {

      if ($DoIt) {
	# Create a request for NetStart
	$seq = substr($sequence, $i, $MaxLength);
	
	if ($i+$MaxLength >= $l) {$Last = 1;}
	
	$ua = new LWP::UserAgent;
	$ua->timeout(3600); # 1 hour !!!
	$ua->agent("AgentName/0.1 " . $ua->agent);
	
	$req = POST 'http://www.cbs.dtu.dk/cgi-bin/nph-webface',
	  [ 
	   configfile => "/home/genome2/www/Public_html/services/NetStart-1.0/NetStart.cf",
	   type => "-at",
	   SEQPASTE => $seq
	  ];
	
	printf STDERR "F%d", 
	  ($i)/$StepLength + 1;
	# Pass request to the user agent and get a response back
	$res = $ua->request($req);
	
	# Retrieve result notification
	if ($res->is_success) 
	  {
	    @content = split('\n',$res->content);
	    foreach $line (@content) {
	      if (($url) = ($line =~ /location\.replace\("(.+)"\)/)) {
		last;
	      }
	    }
	    ($url =~ /^http/) || die "Can't retrieve job id";
	  }
	else {
	  print STDERR "Can't get NetStart result notification";
	  exit(1);
	}
	
	# retrieve result file
	$req = GET $url;
	do {
	  print STDERR ".";
	  $res = $ua->request($req);
	} while ($res->content =~ /javascript/i);
	
	# Check the outcome of the response
	if ($res->is_success) {
	  @content = split(/\n/, $res->content);
	  $flag = 0;
	  foreach $line (@content) {
	    if ($flag == 1 && $line eq "") { $flag = 0;}
	    if ($flag == 1) {
	      if (($n,$r) = ($line =~ /^\s+(\d+)(.+)$/)) {
		if ((($n >= $Olap) || ($i == 0)) && (($n < $MaxLength-$Olap) || $Last  == 1))
		  { printf STARTS "%7d%s\n", $n + $i, $r; }
	      }
	    }
	    if ($line =~ /----------------------/) { $flag = 1; }
	  }
	} 
	else {
	  print STDERR "Bad luck this time\n";
	}
      }

      if ($DoItR) {
	# Create a request for NetStart on reverse strand
	$seq = substr($rsequence, $i, $MaxLength);
	
	$ua = new LWP::UserAgent;
	$ua->agent("AgentName/0.1 " . $ua->agent);
	
	$req = POST 'http://www.cbs.dtu.dk/cgi-bin/nph-webface',
	  [ 
	   configfile => "/home/genome2/www/Public_html/services/NetStart-1.0/NetStart.cf",
	   type => "-at",
	   SEQPASTE => $seq
	  ];

	printf STDERR "R%d", 
	  ($i)/$StepLength + 1;
	# Pass request to the user agent and get a response back
	$res = $ua->request($req);
	
	# Retrieve result notification
	if ($res->is_success) {
	  @content = split('\n',$res->content);
	  foreach $line (@content) {
	    if (($url) = ($line =~ /location\.replace\("(.+)"\)/)) {
	      last;
	    }
	  }
	  ($url =~ /^http/) || die "Can't retrieve job id";
	}
	else {
	  print STDERR "Can't get NetStart result notification";
	  exit(1);
	}
	
	# retrieve result file
	$req = GET $url;
	do {
	  print STDERR ".";
	  $res = $ua->request($req);
	} while ($res->content =~ /javascript/i);
	
	# Check the outcome of the response
	if ($res->is_success) {
	  @content = split(/\n/, $res->content);
	  $flag = 0;
	  foreach $line (@content) {
	    if ($flag == 1 && $line eq "") { $flag = 0;}
	    if ($flag == 1) {
	      if (($n,$r) = ($line =~ /^\s+(\d+)(.+)$/)) {
		if ((($n >= $Olap) || ($i == 0)) && (($n < $MaxLength-$Olap) || $Last  == 1))
		  { printf RSTARTS "%7d%s\n", $n + $i, $r; }
	      }
	    }
	    if ($line =~ /----------------------/) { $flag = 1; }
	  }
	}
	else { print STDERR "Bad luck this time\n"; }
      }
    }

    if ($DoIt)  { close STARTS; }
    if ($DoItR) { close(RSTARTS); }
    printf STDERR "done\n";

    $ua = new LWP::UserAgent;
    $ua->agent("AgentName/0.1 " . $ua->agent);

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Modification for NetGene2: split the sequence in segments of
    # 20Kbp to avoid server limitations then gather result files
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    $Olap = 200;
    $StepLength =  $MaxLength - 2*$Olap;
    $DoIt = 0;

    if (!( -e "$_[0].splices") || (-s "$_[0].splices" == 0))    {$DoIt = 1; }
    if (!( -e "$_[0].splicesR") || (-s "$_[0].splicesR" == 0))  {$DoIt = 1;}


    if ($DoIt) {
      open(FORWARD, ">$_[0].splices") || die "can't create $_[0].splices";
      open(BACKWARD, ">$_[0].splicesR") || die "can't create $_[0].splicesR";

      print FORWARD ">test $l - Made By netgene2\n";
      print BACKWARD ">test $l - Made By netgene2\n";

      printf STDERR "NetGene2 [%d   request(s)]: ",
	($l-2*$Olap+$StepLength-1)/$StepLength;
      $Last = 0;

      for ($i = 0; $i + 2*$Olap < $l; $i += $StepLength) {
	# cut the sequence in segment
	$seq = substr($sequence, $i, $MaxLength);
	
	if ($i+$MaxLength >= $l) {$Last = 1;}
	
	# Create a request for the current segment
	$req = POST 'http://www.cbs.dtu.dk/htbin/nph-webface',
	  [
	   configfile => "/home/genome2/www/Public_html/services/NetGene2/NetGene2.cf", 
	   thing => 'at',
	   SEQNAME => 'test',
	   SEQ => $seq
	  ];
	
	$p = $i/$StepLength;
	printf STDERR "%d", $p + 1;
	# Pass request to the user agent and get a response back
	$res = $ua->request($req);
	
	# Retrieve result notification un
	if ($res->is_success) 
	  {
	    @content = split('\n',$res->content);
	    foreach $line (@content)
	      {
		if (($url) = ($line =~ /location\.replace\("(.+)"\)/))
		  {
		    last;
		  }
	      }
	    ($url =~ /^http/) || die "Can't retrieve job id";
	  }
	else
	  {
	    print STDERR "Can't get netgene2 result notification";
	    exit(1);
	  }
	
	# retrieve result file
	$req = GET $url;
	do
	  {
	    print STDERR ".";
	    $res = $ua->request($req);
	  } while ($res->content =~ /javascript/i);
	
	
	# Check the outcome of the response
	if ($res->is_success) 
	  {
	    @content = split('\n',$res->content);
	    $zob1 = $zob2 = '';
	    foreach $line (@content)
	      {
		if ($line =~ /http:.+\.comp\.score\.txt\.gz/)
		  {
		    ($zob2) = ($res->content =~ /(http:.+\.comp\.score\.txt\.gz)/);
		  }
		elsif ($line =~ /http:.+\.score\.txt\.gz/)
		  {
		    ($zob1) = ($res->content =~ /(http:.+\.score\.txt\.gz)/);
		  }
	      }
	    ($zob1 ne '' && $zob2 ne '') || 
	      die "Can't retrieve location of NetGene2 result files";

	    # retrieve forward file
	    $req = new HTTP::Request 'GET', $zob1;
	    print STDERR "F";
	    $res = $ua->request($req, "$$.splices.Z");
	    if (!$res->is_success) {
	      print STDERR "Problem in forward\n";
	    }
	    # retrieve backward file
	    $req = new HTTP::Request 'GET', $zob2;
	    print STDERR "R";
	    $res = $ua->request($req, "$$.splicesR_$p.Z");
	    if (!$res->is_success) {
	      print STDERR "Problem in backward\n";
	    }
	    # uncompress result files
	    system("gunzip -f $$.splices.Z $$.splicesR_$p.Z");

	    # reformat the results (on for forward part)
	    open(FIN,"$$.splices") || die "can't open $$.splices";
	    $_ = <FIN>; # skip header line
	    while (<FIN>)
	      {
		chomp;
		if (($n,$r,$m) = (/^\s+(\d+)(\s+\w.+\s+)([-\d]+)\s*$/)) {
		  if ((($n >= $Olap) || ($i == 0)) && (($n < $MaxLength-$Olap) || $Last  == 1))
		    { printf FORWARD "%10d%s%s\n", $n + $i, $r, ($m eq '-') ? '-' : $m + $i; }
		}
	      }
	    close(FIN);
	    print STDERR ".";
	    unlink("$$.splices");
	  }
	else
	  {
	    print STDERR "Bad luck this time\n";
	  }
      }

      close(FORWARD);
      $Offset = 0;
      $Last = 1;
      $LastLength = (($l-2*$Olap+$StepLength) % $StepLength)+2*$Olap;

      # rebuild the backward prediction (put result files in reverse order)
      for ($p = int(($l-2*$Olap+$StepLength-1)/$StepLength)-1, $i = 1; $p >=0; $p--) {
	open(FIN,"$$.splicesR_$p") || die "can't open $$.splicesR_$p";
	$_ = <FIN>; # skip header line
	while (<FIN>)
	  {
	    chomp;
	    if (($n,$r,$m) = (/^\s+(\d+)(\s+\w.+\s+)([-\d]+)\s*$/)) {
	      if ((($n >= $Olap) || ($Last == 1)) && 
		  (($n < $LastLength-$Olap) || ($p  == 0))) {
		printf BACKWARD "%10d%s%s\n", $i, $r, ($m eq '-') ? '-' : $m+$Offset;
		$i++;
	      }
	    }
	  }
	$Last = 0;
	$Offset += ($LastLength-2*$Olap);
	$LastLength = $MaxLength;
	close(FIN);
	print STDERR ".";
	unlink("$$.splicesR_$p");
      }
      close(BACKWARD);
      printf STDERR "done\n";
    }
    else { print STDERR "NetGene2: [on disk]\n";}

    # Request for Splice Predictor
    $DoIt = 0;
    if (!( -e "$_[0].spliceP") || (-s "$_[0].spliceP" == 0))    {$DoIt = 1; }
    if (!( -e "$_[0].splicePR") || (-s "$_[0].splicePR" == 0))  {$DoIt = 1;}

    if ($DoIt) {
      $ua = new LWP::UserAgent;
      $ua->agent("AgentName/0.1 " . $ua->agent);

      # Create a request
      $req = POST 'http://bioinformatics.iastate.edu/cgi-bin/sp.cgi',
	[
	 "ufname" => "schiex_id", 
	 "-gdnap" => $sequence,
	 "-l" => "plain ",
	 "-r" => "both",
	 "-f" => "",
	 "-a" => "",
	 "-b" => "",
	 "-pval" => "not set",
	 "-e" => "",
	 "-q" => "",
	 "-sval" => "not set",
	 "-n" => "not set",
	 "-s" => "Arabidopsis",
	 "-m" => "1",
	 "-c" => "All GU and AG sites", #"100% learning set",
	 "-o" => "position"
	];

      print STDERR "SplicePredictor: ";
      # Pass request to the user agent and get a response back
      $res = $ua->request($req);

      # Check the outcome of the response
      if ($res->is_success)
	{
	  print STDERR "done\n";
	  @result = split '\n',$res->content;
	
	  open(FORWARD,">$_[0].spliceP");
	  open(REVERSE,">$_[0].splicePR");
	
	  $file = FORWARD;
	  foreach $line (@result)
	    {
	      $file = REVERSE if ($line =~ m/reverse strand/);
	      printf $file "$line\n"  if ($line =~ m/^(A|D)\s*(<|>|-)/);
	    }
	  close FORWARD;
	  close REVERSE;
	}
      else
	{
	  print STDERR "Bad luck this time\n";
	}
    }
    else { print STDERR "SplicePredictor: [on disk]\n";}
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
    # JER $cmd = sprintf("./eugene -p h %s -g $file > $file.html", join(' ',@params));
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

if (-f "$param1") # just launch eugene on this file
  {
    &eugene($param1, @ARGV);
  }
elsif (-d "$param1") # launch eugene on all fasta files in this folder
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
