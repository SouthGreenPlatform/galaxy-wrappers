#!/usr/bin/perl -w
use strict;

my $hmmline;
my @freq;
my @in;
my @matrix;
my $a;
my $i;
my @p;
my $length;
my $for_print;

# hmm       A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
my @ch2pro=(0, 2, 3, 4, 5, 6, 7, 8, 10,11,12,13,15,16,17,18,19,21,22,24);
#           0  2  3  4  5  6  7  8  10 11 12 13 15 16 17 18 19 21 22 24
# profile   A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z



#run hhmake
system("./hhmake -i temp_files/hhfilter_output.a3m -o temp_files/hhmake_output -pcm 2 -pca 0.5 -pcb 2.5 -cov 20");

#read output file of hhmake
open(HMM,"temp_files/hhmake_output") || die "error opening\n";
while($hmmline=<HMM>)
{
    #this is the start of the part with the interesting values
    if($hmmline=~/\#/)
    {
	#next line gives values that can be used to calculate the background frequencies in the input alignment p(a)
	$hmmline=<HMM>;
	#this should be it...
	if($hmmline=~/^NULL/)
	{
	    #read in values (only one line)
	    @in=split(/\s+/,$hmmline);
	    #get rid of the "NULL" that is the first element of the array
	    shift(@in);
	    #distribute the 20 values over 26 values in a "letter array"
	    for($a=0;$a<20;$a++)
	    {
		$freq[$ch2pro[$a]]=$in[$a];
	    }
	}
    }
    #the following values can be used to calculate the profile from the input alignment
    if($hmmline=~/^\S\s\d+/)
    {
	#read in line by line
	@in=split(/\s+/,$hmmline);
	#get rid of code at the beginning
	shift(@in);
	shift(@in);
	#get rid of number at the end (this is actually the position in the sequence)
	$i=pop(@in);
	#decrement for storage in array
	$i--;
	#put all the values in a matrix
	for($a=0;$a<20;$a++)
	{
	    $matrix[$i][$ch2pro[$a]]=$in[$a];
	}
    }
}
close(HMM);

#length of protein
$length=++$i;

#calculate p(a) (the only purpose of this is being used to calculate p(i,a)
#for($a=0;$a<26;$a++)
#{
#    unless($a==1 ||$a==9 ||$a==14 ||$a==20 ||$a==23 ||$a==25)
#    {
#	$p[$a]=2**($freq[$a]/-1000);
#    }
#}

#calculate p(i,a) and print it to output file
open(PROF,">temp_files/myhmmmake_PCoils_output") || die "error opening\n";
for($i=0;$i<$length;$i++)
{
    for($a=0;$a<26;$a++)
    {

	if($a==1 ||$a==9 ||$a==14 ||$a==20 ||$a==23 ||$a==25 || $matrix[$i][$a] eq "*")
	{
	    print(PROF "0.0000 ");
	}
	else
	{
	    $for_print=2**($matrix[$i][$a]/-1000);
	    printf(PROF "%6.4f ",$for_print);	    
	}
    }
    print(PROF "\n");
}
close(PROF);
