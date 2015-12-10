#!/usr/bin/perl

use strict;

my $in = $ARGV[0];
my $out = $ARGV[1];

open F,$in or die "not possible\n";
open Sav, ">$out" or die "not possible2\n";
my $def;
while (<F>)
{
	my $ligne=$_;

	if ($ligne=~/\<Hit_id\>sp\|.+/)
        {
                $ligne=~s /sp\|/gi\|xa\|sp\|/; # ajout de ca pour que Blast2GO prenne en charge le XML
        }
# blast+ - obliger de rajouter ca pour que le parsing par bioperl se fasse correctement

        if ($ligne=~/\<Iteration_query-ID\>.+\<\/Iteration_query-ID\>/)
        {
        	next;
        }
        if ($ligne=~/\<Iteration_query-def\>(.+)\<\/Iteration_query-def\>/)
        {
        	print Sav "      <Iteration_query-ID>$1</Iteration_query-ID>\n";
        }
# fin commentaire blast
        print Sav "$ligne";
}
close Sav;
close F;
