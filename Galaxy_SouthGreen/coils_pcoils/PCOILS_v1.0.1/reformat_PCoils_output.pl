#!/usr/bin/perl -w
use strict;

my $line;
my $ignore;
my $full;

unless(defined($ARGV[2]) && defined($ARGV[1])) {die "usage: perl reformat_PCoils_output.pl <file> <ignore> <full>\n";} #ignore: number of residues before classification, full: number of classified residues

$ignore=$ARGV[1];
$full=$ARGV[2];

open(FILE,$ARGV[0]) || die "usage: perl reformat_PCoils_output.pl [file_to_be_reformatted]\n";
while($line=<FILE>)
{
    if($line=~/(\S+)\s+\(/)
    {
	if($ignore<=0 && $ignore>(-$full))
	{
	    print "$1\n";
	}
	$ignore--;
    }
}
close(FILE);
