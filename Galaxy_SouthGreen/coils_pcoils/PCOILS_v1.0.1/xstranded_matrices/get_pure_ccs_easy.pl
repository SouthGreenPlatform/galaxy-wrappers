#!/usr/bin/perl -w

use strict;

#this program extracts the pure coiled coil sequences from the Socket output files

my @ccarray;
my $ccarray;
my $line;
my $name;
my $chain;
my $sequence;
my $register;
my $neglect;
my $real_register;
my $sequence_line;
my $register_line;
my $rl_length;

open(CC,"cc_data") || die "opening error of cc_data\n";
while($line=<CC>)
{
    chomp $line;
    unshift(@ccarray,$line);
}
close(CC);

#counts valid coiled coils that were treated
my $name_counter=0;


#loop over all coiled coils
foreach $ccarray (@ccarray)
{
    #open coiled coil file in SOCKET format for reading
    open (FILE_POINTER, "../../socket/only_coiled_coils/pdb$ccarray.ent.output") ||
    open (FILE_POINTER, "../../socket/only_coiled_coils_latest/$ccarray.output") || die "error opening only_coiled_coils for $ccarray\n";

    $name=$ccarray;

    #read each line of the coiled coil file
    while($line=<FILE_POINTER>)
    {
	#pattern match for "extent of coiled coil packing: 162-190:F", for example
	if($line=~/extent\sof\scoiled\scoil\spacking:\s+(\-?\d+)\s?\S\s?(\d+)\s?:\s?(\w?)/)
	{
	    undef $chain;
	    $chain=$3;
	    $line=<FILE_POINTER>;
	    if($line=~/(sequence\s+\w+)/)
	    {
		$sequence=$1;
		$line=<FILE_POINTER>;
		if($line=~/(register\s+\w+)/)
		{
		    $register=$1;
		    if($line=~/(register\s+)(\w+)/)
		    {
			$neglect=length($1);
			#print "neglect: $neglect\n";
			$real_register=$2;
		    }
		    if(defined($chain))
		    {
			print ">gi|$name|$chain\n";
		    }
		    else
		    {
			print ">gi|$name|X\n";
		    }
		    if($sequence=~/^.{$neglect}(\S+)/)
		    {
			$sequence_line="sequence:$1";
			$register_line="register:$real_register";
			$rl_length=length($register_line);
			if($sequence_line=~/^(.{$rl_length})/)
			{
			    print "$1\n";
			    print "$register_line\n";
			}
		    }
		}
	    }
	}
    }
}

close(FILE_POINTER);
