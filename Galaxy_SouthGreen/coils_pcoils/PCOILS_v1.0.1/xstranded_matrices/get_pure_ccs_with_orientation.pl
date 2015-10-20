#!/usr/bin/perl

#this program extracts the pure coiled coil sequences from the Socket output files
#it works like "get_pure_ccs_with_orientation.pl", but it also gives the orientation and number of strands


open(TWO,">temp_files/2stranded");
open(THREE,">temp_files/3stranded");
open(FOUR,">temp_files/4stranded");

$aa2_counter=0;
$aa3_counter=0;
$aa4_counter=0;


#put all the file names of only_coiled_coil directory in @name
opendir (ONLY_CC, "../../socket/only_coiled_coils") || die "Cannot opendir only_coiled_coils: $!";
@filename = readdir(ONLY_CC);
closedir(ONLY_CC);

#counts valid coiled coils that were treated
$name_counter=0;

#loop over all coiled coils
foreach $filename (@filename)
{
    #open coiled coil file in SOCKET format for reading
    open (FILE_POINTER, "../../socket/only_coiled_coils/$filename");

    #extract the PDB code from the filename
    if($filename=~/pdb(\w{4})/)
    {
	#the PDB code is put in the variable $name
	$name=$1;
    }

    #read each line of the coiled coil file
    while($line=<FILE_POINTER>)
    {
	if($line=~/^$name/)
	{
	    if($line=~/\((parallel)/)
	    {
		$orientation=$1;
	    }
	    if($line=~/(antiparallel)/)
	    {
		$orientation=$1;
	    }
	    if($line=~/(\d\-stranded)/)
	    {
		$strand_number=$1;
	    }
	}
	#pattern match for "extent of coiled coil packing: 162-190:F", for example
	if($line=~/extent\sof\scoiled\scoil\spacking:\s+(\-?\d+)\s?\S\s?(\d+)\s?:\s?(\w?)/)
	{
	    #print "extent of coiled coil packing checked\n";
	    #in this example "162" would be the first position of the coiled coil...
	    $beginning=$1;
	    #...and "190" the last position...
	    $end=$2;
	    #...while "F" is the name of the chain
	    $chain=$3;
	    #the first and last positions will not be used except for calculating the length of the cc region

	    #print $name, ": ", $beginning, "-", $end, ":", $chain, "\n";

	    #now the file with the complete sequence of the chain is opened for reading
	    open (SOCKET_CC, "../../fasta/$name.$chain.fasta")
		|| open (SOCKET_CC, "../../fasta/$name.X.fasta");

	    #the first line is the ID - we don't care for it!
	    $get_socket_cc=<SOCKET_CC>;
	    #the second line is the sequence
	    $get_socket_cc=<SOCKET_CC>;
	    $length_get_socket_cc=length $get_socket_cc;
	    $length_get_socket_cc--;
	    #print "NEW get_socket_cc\n";
	    #the sequence is stored
	    #print $get_socket_cc;
	    close (SOCKET_CC);
	}

	#pattern match for "sequence RIRRERNKMAAAKSRNRRRELTDTLQAETDQLEDEKSALQTEIANLLKEKE", for example
	#this is not the complete sequence as in $get_socket_cc, but only the helical part of it
	#it is therefore longer than the coiled coil region
	if($line=~/sequence\s+(\w+)/)
	{
	    #print "sequence checked\n";
	    #the sequence is stored
	    $socket_sequence=$1;

	    #print "socket_sequence: ", $socket_sequence, "\n";

	    #let's see what residues are missing at the beginning
	    #if we only consider the helical part of the sequence
	    if($get_socket_cc=~/(\w*)$socket_sequence/)
	    {
		#how many residues are non-helical? (how many of them do we have?)
		$ignore_for_helix=length $1;
		#print "ignore: ", $ignore_for_helix, "\n";
		#print "cc_length: ", $end-$beginning+1, "\n";
	    }
	    
	}

	#pattern match for something like "register                    "
	#and
	#"abcdefgabcdefgabcdefgabcdefga"
	if($line=~/(register\s+)(\w+)/)
	{
	    undef %hash;
	    #print "register checked\n";
	    #how much later than the helical part does the coiled coil region start?
	    #the obtained value will be 9 too high, because "register " is considered as well...
	    #...and it has 9 characters
	    $ignore_for_reg=length $1;
	    #print "ignore_for_reg: ", length $1, "\n";

	    #the starting point of the coiled coil region is the sum of
	    #the non-helical residues and the only-helical-but-not-in-register residues
	    #of course, one has to subtract 9 from this sum
	    $starting_point=$ignore_for_helix+$ignore_for_reg-9;
	    #print "starting_point changed\n";

	    #the register itself is stored in an array
	    @register = split //, $2;

	    

	    #how many residues are considered?
	    #initialisation
	    $residue_counter=0;
	    
	    #increment the name counter
	    $name_counter++;

	    #the sequence is put in an array
	    undef(@complete_sequence);	    
	    @complete_sequence = split //, $get_socket_cc;

	    #for loop taking into account different starting positions for the coiled coil
	    $limit=$end-$beginning+$starting_point;
	    #for($i=$starting_point;$i<=$limit;$i++)
	    #print $name, "\n";

	    #we would like to consider only the coiled coils that are specified in the file "strict_pdb_ccs"
	    open(STRICT, "get_pure_ccs_easy.output") || die "opening error\n";
	    while($strictline=<STRICT>)
	    {
		if($orientation eq "parallel" && $strand_number eq "2-stranded")
		{
		    $which_file="TWO";
		}
		if($orientation eq "parallel" && $strand_number eq "3-stranded")
		{
		    $which_file="THREE";
		}
		if($orientation eq "parallel" && $strand_number eq "4-stranded")
		{
		    $which_file="FOUR";
		}

		#if the pdb code is specified in the "strict" file, the coiled coil fragment is printed out
		if($strictline=~/$name\S$chain/i)
		{
		    unless(defined($hash{$name.$chain}))
		    {
			printf($which_file ">gi|$name|$chain\n");
			#print "$orientation $strand_number\n";
			printf ($which_file "register:");
			for($i=0;$i<=$end-$beginning;$i++)
			{
			    printf($which_file "$register[$i]");
			    if($which_file eq "TWO") {$aa2_counter++;}
			    if($which_file eq "THREE") {$aa3_counter++;}
			    if($which_file eq "FOUR") {$aa4_counter++;}
			}
			printf($which_file "\n");	
			printf($which_file "sequence:");
			for($i=$starting_point;$i<=$limit;$i++)
			{
			    printf($which_file "$complete_sequence[$i]");
			}
			printf($which_file "\n");
		    }
		    $hash{$name.$chain}="used";
		}
	    }
	    close(STRICT);

	}
	    
	
    }

    close (FILE_POINTER);
}

printf(STDERR "number of two-stranded amino acids: $aa2_counter\n");
printf(STDERR "number of three-stranded amino acids: $aa3_counter\n");
printf(STDERR "number of four-stranded amino acids: $aa4_counter\n");

close(TWO);
close(THREE);
close(FOUR);
