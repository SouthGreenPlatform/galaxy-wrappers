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
# $Id: egn_getRepeat4eugene.pl,v 1.1 2005-03-17 13:41:23 cros Exp $
# ------------------------------------------------------------------
# File:     egn_getRepeat4eugene.pl
# Contents: build .ig file for eugene
# ------------------------------------------------------------------

use Getopt::Long;
use ParamParser;
use warnings;
use strict;
use IO::Handle;

=pod

=head1

=head1 NAME

  getRepeat4eugene - Build .ig file for eugene

=head1 MANDATORY

  -s, --seqs     [FILE NAME(S)]
        Fasta sequence file(s).
        The sequence(s) must be reformat by "ref"

  -d, --database [FILE NAME]
        Multiple fasta file (repeats).
        Must be reformat by "ref"

=head1 EXAMPLES

  getRepeat4eugene.pl -s "*.fasta" -d repeat.fasta

=head1 Copyright (c) 2004 by INRA. All rights reserved.

=head1 Mail : eugene@ossau.toulouse.inra.fr

=head1

=cut

my $rh_param =  New ParamParser('GETOPTLONG',("seqs=s","database=s"));

MAIN:
{
    # Parameters
    $rh_param->SetUsage(my $usage=sub { &Usage(); } );
    $rh_param->AssertFileExists('database');
    $rh_param->AssertDefined('seqs');
    $rh_param->AssertDefined('database');
    my($tfafile) = $rh_param->Get('seqs');
    my($dbfile)  = $rh_param->Get('database');

    print STDERR "#\n# Get repeat for eugene\n#\n";

    ## PRESSDB
    print STDERR "pressdb...";
    if (!( -e "$dbfile.csq") || (-s "$dbfile.csq" == 0)  ||
	!( -e "$dbfile.nhd") || (-s "$dbfile.nhd" == 0)  ||
	!( -e "$dbfile.ntb") || (-s "$dbfile.ntb" == 0))
    {
	system "pressdb $dbfile >& /dev/null";
	print STDERR "done\n";
    }
    else { print STDERR "[ON DISK]\n"; }

    my(@a_tfaFiles) = `ls $tfafile`;
    for (my $i=0; $i<=$#a_tfaFiles; $i++)
    {
	chomp $a_tfaFiles[$i];
	my $tfa = $a_tfaFiles[$i];
	print STDERR "- ".($i+1).":$tfa\t";

	## WU
	print STDERR "wu.";
	if (!( -e "$tfa.wu") || (-s "$tfa.wu" == 0))
	{
	    system "/www/meth.DNA/meth.DNA.WU-blast.ok.pl LIB=$dbfile IN=$tfa OUT=$tfa.wu PARAM=\"M=5 N=-5 E=1000 S=100 S2=100\"";
	    print STDERR "done...";
	}
	else { print STDERR "[ON DISK]...";}

	## Filter
	print STDERR "fil.";
	if (!( -e "$tfa.fil") || (-s "$tfa.fil" == 0))
	{
	    system "/www/meth.AA_DNA/meth.AA_DNA.Filter.ok.pl FILE=$tfa.wu IDENTITY=100 PCS=100 LIGHT=YES | grep \\< > $tfa.fil";
	    print STDERR "done...";
	}
	else { print STDERR "[ON DISK]...";}
	
	## Build IG (sort and no overlap)
	print STDERR "ig.";
	if (!( -e "$tfa.ig") || (-s "$tfa.ig" == 0))
	{
	    my @a_hit      = ();
	    my @a_hitmerge = ();

	    open(TMP,">$tfa.filtmp") or die "FATAL ERROR >$tfa.filtmp<\n";
	    @a_hit = `cut -d" " -f2,3 $tfa.fil`;
	    for(my $i=0; $i<=$#a_hit; $i++)
	    {
		my ($d,  $f)  = $a_hit[$i] =~ /^([0-9]+) ([0-9]+)/;
		if ($d > $f)
		{
		    my $tmp = $d;
		    $d = $f;
		    $f = $tmp;
		}
		print TMP "$d $f\n";
	    }
	    @a_hit = ();
	    close TMP;
	
	    @a_hit = `sort -n $tfa.filtmp | uniq`;
	    for(my $i=0; $i<=$#a_hit; $i++)
	    {
		my ($d,  $f)  = $a_hit[$i] =~ /^([0-9]+) ([0-9]+)/;
		if($i == 0) { $a_hitmerge[0] = "$d $f\n"; }
		else
		{
		    my ($dm, $fm) = $a_hitmerge[-1] =~ /^([0-9]+) ([0-9]+)/;
		    if($fm > $d  &&  $f > $fm)
		    {
			$a_hitmerge[-1] = "$dm $f\n";
		    }
		    elsif($d > $fm)
		    {
			push(@a_hitmerge, $a_hit[$i]);
		    }
		}
	    }
	    open(IG,">$tfa.ig") or die "FATAL ERROR >$tfa.ig<\n";
	    if($#a_hitmerge == -1) { print STDERR "-No repeat-."; }
	    else {
		for(my $i=0; $i<=$#a_hitmerge; $i++)
		{
		    print IG $a_hitmerge[$i];
		}
	    }
	    close IG;
	    print STDERR "done\n";
	}
	else { print STDERR "[ON DISK].done\n";}
	unlink ("$tfa.filtmp");
	unlink ("$tfa.wu");
    }
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
        else                   { print "$u"; }
    }
}
#### End Usage ####
