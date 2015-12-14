#!/usr/bin/perl 
use strict;

my $xml4blast2go = shift;
my $xml = shift;

open F, "$xml" or die "not possible\n";
open Sav, ">$xml4blast2go" or die "not possible2\n";
my $def;

#	<Hit_id>blast|to|go|format|</Hit_id>
#	<Hit_def>sp|P0CH30|RING1_GOSHI E3 ubiquitin-protein ligase RING1 OS=Gossypium hirsutum GN=RING1 PE=2 SV=1</Hit_def>

while (<F>){
    my $ligne=$_; 
    
    if ($ligne=~/\<Hit_id\>.*\<\/Hit_id\>/){
    	next;
    }
    if ($ligne=~/(.*)\<Hit_def\>(.*)\<\/Hit_def\>/){
    	print Sav "$1<Hit_id>$2</Hit_id>\n";
    }
    
	# blast+ - obliger de rajouter ca pour que le parsing par bioperl se fasse correctement 

    if ($ligne=~/\<Iteration_query-ID\>.+\<\/Iteration_query-ID\>/){
	next;
    }
    if ($ligne=~/\<Iteration_query-def\>(.+)\<\/Iteration_query-def\>/){
	print Sav "      <Iteration_query-ID>$1</Iteration_query-ID>\n"; 
    }                                  
    print Sav $ligne;
}

close Sav;
close F;
