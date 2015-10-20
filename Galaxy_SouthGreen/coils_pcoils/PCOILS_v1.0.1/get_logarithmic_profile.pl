#!/usr/bin/perl -w
use strict;

#creates a coiled coil profile in logarithmic form for the matrix in "test"

my $line;

open(P,$ARGV[0]) || die "usage: perl get_logarithmic_profile.pl <non_log_coils_matrix>\n";

while($line=<P>)
{
    if($line=~/(\w)\s(\d+.\d+)\s(\d+.\d+)\s(\d+.\d+)\s(\d+.\d+)\s(\d+.\d+)\s(\d+.\d+)\s(\d+.\d+)/)
    {
	print $1, " ";
	if($2>0)
	{
	    print log $2, " ";
	}
	else
	{
	    print "-999999.0 ";
	}
	if($3>0)
	{
	    print log $3, " ";
	}
	else
	{
	    print "-999999.0 ";
	}
	if($4>0)
	{
	    print log $4, " ";
	}
	else
	{
	    print "-999999.0 ";
	}
	if($5>0)
	{
	    print log $5, " ";
	}
	else
	{
	    print "-999999.0 ";
	}
	if($6>0)
	{
	    print log $6, " ";
	}
	else
	{
	    print "-999999.0 ";
	}
	if($7>0)
	{
	    print log $7, " ";
	}
	else
	{
	    print "-999999.0 ";
	}
	if($8>0)
	{
	    print log $8, "\n";
	}
	else
	{
	    print "-999999.0\n";
	}
    }
    else
    {
	print $line;
    }
}
