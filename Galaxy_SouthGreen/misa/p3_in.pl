#!/usr/bin/perl 
# Author: Thomas Thiel
# Program name: primer3_in.pl
# Description: creates a PRIMER3 input file based on SSR search results


open (IN,"<$ARGV[0]") || die ("\nError: Couldn't open misa.pl results file (*.misa) !\n\n");

my $filename = $ARGV[0];
$filename =~ s/\.misa//;
open (SRC,"<$filename") || die ("\nError: Couldn't open source file containing original FASTA sequences !\n\n");
open (OUT,">$filename.p3in");

undef $/;
$in = <IN>;
study $in;

$/= ">";

my $count;
while (<SRC>)
  {
  next unless (my ($id,$seq) = /(.*?)\n(.*)/s);

  $seq =~ s/[\d\s>]//g;#remove digits, spaces, line breaks,...
  while ($in =~ /$id\t(\d+)\t\S+\t\S+\t(\d+)\t(\d+)/g)
    {
      my ($ssr_nr,$size,$start) = ($1,$2,$3);
      $count++;
      print OUT "PRIMER_SEQUENCE_ID=$id"."_$ssr_nr\nSEQUENCE=$seq\n";
      print OUT "PRIMER_PRODUCT_SIZE_RANGE=100-280\n";
      print OUT "PRIMER_OPT_SIZE=17\n";
      print OUT "PRIMER_MIN_SIZE=15\n";
      print OUT "PRIMER_MAX_SIZE=23\n";
      print OUT "PRIMER_OPT_TM=56.0\n";
      print OUT "PRIMER_MIN_TM=50.0\n";
      print OUT "PRIMER_MAX_TM=63.0\n";
      print OUT "PRIMER_MIN_GC=35\n";
      print OUT "PRIMER_MAX_DIFF_TM=1.5\n";
      print OUT "PRIMER_SELF_ANY=5\n";
      print OUT "PRIMER_SELF_END=2\n";
      print OUT "PRIMER_MAX_POLY_X=4\n";
      print OUT "PRIMER_WT_COMPL_END=10\n";
      print OUT "PRIMER_PAIR_WT_DIFF_TM=5\n";
      print OUT "PRIMER_LIBERAL_BASE=1\n";
      print OUT "TARGET=",$start-3,",",$size+6,"\n=\n"
    };
  };
#print "\n$count records created.\n";
