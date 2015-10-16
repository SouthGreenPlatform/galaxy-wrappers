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
# $Id: egn_compseq.pl,v 1.1 2004-09-30 11:43:41 bardou Exp $
# ------------------------------------------------------------------
# File:     egn_compseq.pl
# Contents: counts composition of dimer/trimer/etc words
# ------------------------------------------------------------------

use Getopt::Long;
use ParamParser;
use warnings;
use strict;
use IO::Handle;

=pod

=head1

=head1 NAME

  egn_compseq - Counts the composition of dimer/trimer/etc words in a sequence

=head1 MANDATORY

  -i, --in    [FILE NAME(S)]
        Fasta sequence file(s).

  -s, --size  [INT]
        Word size to consider (e.g. 2=dimer).

  -o, --out   [FILE NAME]
        Output file [<sequence>.comp]

=head1 EXAMPLES

  egn_compseq.pl -i "*.fasta" -s 10 -o seq.comp

=head1 Copyright (c) 2004 by INRA. All rights reserved.

=head1 Mail : eugene@ossau.toulouse.inra.fr

=head1

=cut

# Parameters
my $rh_param =  New ParamParser('GETOPTLONG',("in=s","size=i","out=s"));

# Print in output file only word with count > X
my $printFilter   = 5;

# Delete keys with value = 1 when keys > X
my $warningMemory = 5000000;

MAIN:
{
    # Parameters
    $rh_param->SetUsage(my $usage=sub { &Usage(); } );
    #$rh_param->AssertFileExists('in');
    $rh_param->AssertDefined('in');
    $rh_param->AssertInteger('size');
    $rh_param->AssertDefined('size');
    $rh_param->AssertDefined('out');

    my($infile) = $rh_param->Get('in');
    my($size)   = $rh_param->Get('size');
    my($outfile) = $rh_param->Get('out');

    my(@a_inFiles) = `ls $infile`;
    my %h_comp = ();
    my $c = 0;

    for (my $i=0; $i<=$#a_inFiles; $i++)
    {
	chomp $a_inFiles[$i];
	print STDERR ($i+1).":".$a_inFiles[$i]."\n";

	my $seq = &GetSeq($a_inFiles[$i]);
	my $j   = 0;
	while ( $j  <= (length($seq)-$size) )
	{
	    my $w = substr($seq, $j, $size);
	    if ($w !~ /N/)
	    {
		$h_comp{$w}++;
		$c++;
	    }
	    if ($c > $warningMemory)
	    {
		print STDERR "\nMemory WARNING : $c keys (>$warningMemory)".
		    " => delete keys with value = 1...";
		foreach my $k (keys %h_comp)
		{
		    if($h_comp{$k} == 1) { delete $h_comp{$k}; }
		}
		$c = 0;
		print STDERR "done\n\n";
	    }
	    $j++;
	}
    }

    open(OUT,">$outfile") or die "FATAL ERROR : can't create >$outfile<\n";
    $c = 0;
    foreach my $k (keys(%h_comp))
    {
	if ($h_comp{$k} > $printFilter)
	{
	    print OUT "$k   $h_comp{$k}\n";
	    $c++;
	}
    }
    print STDERR "\nKeys number (with value > $printFilter) : $c\n";
    close OUT;
}


#######################
#        SUB          #
#######################

=head2 Procedure Usage

 Title          : Print usage
 Usage          : &Usage()
 Prerequisite   : none
 Function       : Print usage on stdout
 Returns        : none
 Args           : none
 Globals        : none

=cut

sub Usage
{
    my @usage = `pod2text $0`;
    foreach my $u ( @usage )
    {
        if ($u =~ /Procedure/) { last; }
        else                { print "$u"; }
    }
}
#### End Usage ####

=head2 Procedure GetSeq

 Title          : Get sequence
 Usage          : &GetSeq()
 Prerequisite   : none
 Returns        : date sequence
 Args           : fasta filename
 Globals        : none

=cut

sub GetSeq
{
    my($seq) = @_;
    my($lign);
    my($seq_seq) = "";
    my($multifasta) = -1;

    open(SEQ,"$seq") or die "FATAL ERROR >$seq<\n";
    while($lign=<SEQ>)
    {
        chomp($lign);
        if ( $lign =~ /^>/ ) { $multifasta++;     }
        else                 { $seq_seq .= $lign; }
        if ( $multifasta )
        {
            print STDERR "FATAL ERROR: this script cannot be used".
                " on a multifasta file ($seq)\n";
            exit;
        }
    }
    close(SEQ);
    $seq_seq =~ s/[^ATGCNX]//gi;
    return "\U$seq_seq\E";
}
#### End GetSeq ####
