#!/usr/bin/perl

#treats files with sequences in FASTA format
#sequences are put in one line

unless(defined($ARGV[0])) {die "usage: perl manipulate_pdbaa.pl <file_to_be_manipulated>\n"};

#open file to be opened (e.g. pdb database in FASTA format)
open (PDB, $ARGV[0]);

#read one line and print it (this is only done for the first line!)
#the reason for this is simply the handling of newline characters
$line=<PDB>;
print $line;

#for all other lines:
#they are read...
while ($line=<PDB>)
{
    #...non-sequence lines are printed...
    if($line=~/>/)
    {
	print "\n",$line;
    }
    else
    {
	#...sequence lines are put together and printed...
	chomp $line;
	print uc($line);
    }
	
}

#last newline (don't worry about the newline handling!)
print "\n";
