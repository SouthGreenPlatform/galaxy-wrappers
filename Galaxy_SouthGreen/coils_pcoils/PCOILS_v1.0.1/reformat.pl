#! /usr/bin/perl -w
#
use strict;

my $numres=100;             # number of residues per line
my $DESCLEN=1000;           # maximum number of characters in nameline
my $ARGC=scalar(@ARGV);
if ($ARGC<2) {
    print("
Read a multiple alignment in one format and write it in another format
Usage: reformat.pl [informat] [outformat] infile outfile [options] 
  or   reformat.pl [informat] [outformat] 'fileglob' .ext [options] 

Available input formats:
   fas:     aligned fasta; lower and upper case equivalent, '.' and '-' equivalent
   a2m:     aligned fasta; inserts: lower case, matches: upper case, deletes: '-', gaps aligned to inserts: '.'
   a3m:     like a2m, but gaps aligned to inserts MAY be omitted
   sto:     Stockholm format; sequences in several blocks with sequence name at beginning of line (hmmer output)
   clu:     Clustal format; sequences in several blocks with sequence name at beginning of line
   psi:     psiblast-readable format 
Available output formats:
   fas:     aligned fasta; all gaps '-'
   a2m:     aligned fasta; inserts: lower case, matches: upper case, deletes: '-', gaps aligned to inserts: '.'
   a3m:     like a2m, but gaps aligned to inserts are omitted
   sto:     Stockholm format; sequences in just one block, one line per sequence
   psi:     psiblast-readable format (like sto with -M first -r: remove all columns without residue in first sequence)
   clu:     clustal format
If no input or output format is given the file extension is interpreted as format specification ('aln' as 'clu')

Options:
  -v int    verbose mode (0:off, 1:on)
  -num      add number prefix to sequence names: '0|name', '1|name' etc
  -noss     remove secondary structure sequences (beginning with >aa_ or >ss_)
  -M int    make all columns with less than X% gaps match columns (for output format a2m or a3m)
  -M first  make all columns with residue in first seuqence match columns (for output format a2m or a3m)
  -r        remove all lower case residues (insert states) (AFTER -M option has been processed)
  -r int    remove all lower case columns with more than X% gaps
  -g ''     suppress all gaps
  -g '-'    write all gaps as '-'
  -uc       write all residues in upper case (AFTER all other options have been processed)
  -lc       write all residues in lower case (AFTER all other options have been processed)
  -l        number of residues per line (for CLUSTAL, FASTA, A2M, A3M formats) (default=$numres)

Examples: reformat.pl 1hjra.a3m 1hjra.a2m   (same as reformat.pl a3m a2m 1hjra.a3m 1hjra.a2m)
          reformat.pl test.a3m test.fas -num -r 90
          reformat.pl fas sto '*.fasta' .stockholm
\n");
#  clu:  clustal format (hmmer output)
#  sel:  Selex format
#  phy:  Phylip format
    exit(1);
}

my $informat="";
my $outformat="";
my $infile="";
my $outfile="";
my $num=0;                    # don't use sequence number as prefix: '>n|name'
my $noss=0;                   # don't remove secondary structure
my $line;
my $options="";
my @names;   # names of sequences read in
my @seqs;    # residues of sequences read in
my $n;       # number of sequences read in
my $k;       # counts sequences for output
my $remove_inserts=0;
my $remove_gapped=0;
my $matchmode="no";   # do not change capitalization
my $match_gaprule=0;  # columns with less than this percentage of gaps will be match columns
my $v=2;
my $update=0;

my @informats=  ("fas","a2m","a3m","sto","psi","clu");
my @outformats= ("fas","a2m","a3m","sto","psi","clu","ufas");
my $found;
my $element;
my $gap="default";
my $case="default";

# Process options
for (my $i=0; $i<$ARGC; $i++) {$options.=" $ARGV[$i]";}
if ($options=~s/ -i\s+(\S+)//g) {$infile=$1;}
if ($options=~s/ -o\s+(\S+)//g) {$outfile=$1;}
if ($options=~s/ -num//g)    {$num=1;}
if ($options=~s/ -noss//g)   {$noss=1;}
if ($options=~s/ -g\s+\'?(\S*)\'?//g) {$gap=$1;}
if ($options=~s/ -r\s+(\d+)//g) {$remove_gapped=$1;}
if ($options=~s/ -r//g)     {$remove_inserts=1;}
if ($options=~s/ -lc//g)    {$case="lc";}
if ($options=~s/ -uc//g)    {$case="uc";}
if ($options=~s/ -v\s*(\d+)//g) {$v=$1;}
if ($options=~s/ -v//g) {$v=2;}
if ($options=~s/ -M\s+(\d+)//g) {$matchmode="gaprule"; $match_gaprule=$1;}
if ($options=~s/ -M\s+first//g) {$matchmode="first";   $match_gaprule=$1;}
if ($options=~s/ -u//g) {$update=1;}
if ($options=~s/ -l\s+(\S+)//g) {$numres=$1;}
 
# Assign informat, outformat, infile, and outfile
if ($outfile eq "") {
    if ($options=~s/(\S+)\s*$//) {
	$outfile=$1;
    } else {
	die("Error: no output file given: '$options'\n");
    }
}
if ($infile eq "") {
    if ($options=~s/(\S+)\s*$//) {
	$infile=$1;
    } else {
	die("Error: no input file given: '$options'\n");
    }
}
if ($options=~s/(\S+)\s*$//) {
    $outformat=$1;
} else {
    if ($outfile=~/\S*\.(\S+?)$/) {
	$outformat=lc($1);  
	if    ($outformat eq "aln")   {$outformat="clu";}
	elsif ($outformat eq "fa")    {$outformat="fas";}
	elsif ($outformat eq "fasta") {$outformat="fas";}
    } else {
	print ("Using FASTA output format: '$options'\n"); $outformat="fas";
    }
}
if ($options=~s/(\S+)\s*$//) {
    $informat=$1;
} else {
    if ($infile=~/\S*\.(\S+?)$/) {
	$informat=lc($1);
	if    ($informat eq "aln")   {$informat="clu";}
	elsif ($informat eq "fa")    {$informat="fas";}
	elsif ($informat eq "fasta") {$informat="fas";}
    } else {
	print ("Using FASTA input format: '$options'\n"); $informat="fas";
    }
}


# Warn if unknown options found
if ($options!~/^\s*$/) {
    $options=~s/^\s*(.*?)\s*$/$1/g; 
    print("\nWARNING: unknown options '$options'\n");
}

# Check if input and output formats are valid
$found=0;
foreach $element (@informats) {if ($informat eq $element) {$found=1; last;}} 
if(!$found) {die("\nError: $informat is not a valid input format option\n");}
$found=0;
foreach $element (@outformats) {if ($outformat eq $element) {$found=1; last;}} 
if(!$found) {die("\nError: $outformat is not a valid output format option\n");}

#if($outformat eq "psi") {
#   $remove_inserts=1;
#}
if($outformat eq "ufas") {$gap=""; $outformat="fas";}


if ($infile=~/\*/ || $outfile=~/^\./) # if infile contains '*' or outfile is just an extension
{
    $outfile=~/.*\.(\S*)$/;
    my $outext=$1;
    my @infiles=glob($infile);
    printf("%i files to reformat\n",scalar(@infiles));
    foreach $infile (@infiles) 
    {
	if ($infile!~/(\S+)\.\S+/) {$infile=~/(\S+)/}
	$outfile="$1.$outext";
	if ($update && -e $outfile) {next;}
	if ($v>=3) {print("Reformatting $infile from $informat to $outformat ...\n");} 
	&reformat($infile,$outfile);
    }
    exit;
}
else
{
    if ($v>=3) {print("Reformatting $infile from $informat to $outformat ...\n");} 
    exit &reformat($infile,$outfile); 
}


################################################################################################
# Reformat a single file
################################################################################################
sub reformat()
{
    $infile=$_[0];

################################################################################################
# Input part
################################################################################################

    open (INFILE,"<$infile") or die ("ERROR: cannot open $infile: $!\n");
    my $skip=0;

    # Read a2m, a3m, fas format
    if ($informat eq "fas" || $informat eq "a2m" || $informat eq "a3m")
    { 
	$n=-1;
	while ($line = <INFILE>) 
	{
	    if ($line=~/^>/) # nameline detected
	    {
		if ($noss && $line=~/^>(aa|ss)_/) { # do not read in >ss_ and >aa_ sequences
		    $skip=1;
		} else {
		    $skip=0; 
		    $line=~/^>(.*)/;
		    $n++;
		    $names[$n]=$1;
		    $seqs[$n]="";
		}
	    }
	    elsif ($line=~/^\#/) {next;} # commentary detected
	    elsif (!$skip)
	    {
		chomp($line);
		$seqs[$n].=$line;
	    }
	}
	$n++;
    }
    
    # Read Stockholm format
    elsif ($informat eq "sto")
    {
	my %seqhash;
	my $name;
	my $first_block=1;
	
	$n=0;
	while ($line = <INFILE>)
	{
	    $line=~s/\s+/ /g;
	    if ($line=~/^\#/) {next;}                    # skip commentary and blank line
	    if ($line=~/^\/\//) {last;}                   # reached end of file
	    if ($line=~/^\s*$/){$first_block=0; next;}             # new sequence block starts
	    if ($line!~/^\s*(\S+)\s+(\S+)/) {
		die ("\nERROR found in stockholm format: $!");
	    }
	    if (!(exists $seqhash{$1})) 
	    {
		$line=~/^\s*(\S+)\s+(\S+)/;
		$names[$n]=$1;
		$seqs[$n]=$2;
		$seqhash{$1}=$n++;
		$first_block=1;
	    }
	    else
	    {
		if ($first_block) {die ("\nERROR: sequence $1 appears more than once per block\n");}
		$seqs[$seqhash{$1}].=$2;
	    }
#	    printf("%3i %s\n",$n-1,$seqs[$n-1]);
	}
    }

    elsif ($informat eq "clu")
    {
	my $residues_per_line=60; # default number of characters for namefield 
                             # (only needed if no gap between name and residues in first sequence -> SMART)
	my $block=1;         # number of block currently read 
	my $name;
	my $residues;
	$n=0;                # sequence number of first block
	$k=0;                # sequence index to zero for first block

	while ($line = <INFILE>)
	{
#	    print($namefieldlen.$line);
	    if ($line=~/CLUSTAL/i) {next;}
	    if ($line=~/^\#/) {next;}                     # skip commentary and blank line
	    if ($line=~/^\/\//) {last;}                   # reached end of file
	    if ($line=~/^\s*$/){                          # new sequence block starts
		if ($k) {
		    if ($n && $n!=$k) {die("\nError: different number of sequences in blocks 1 and $block of $infile\n");} 
		    $block++;
		    $n=$k;
		    $k=0;  # set sequence index to zero for next block to read
		}
		next;
	    }
	    if ($line!~/^(\S+)\s+([ a-zA-Z.-]+?)\s*\d*$/) {
		if ($line=~/^[*.: ]*$/) {next;}
		chomp($line);
		if ($line!~/^(\S+)([a-zA-Z.-]{$residues_per_line})\s*\d*$/) {
		    die ("\nError found in CLUSTAL format in $infile, line $.: '$line'\n");
		} 
		$name=$1;
		$residues=$2;
		print("WARNING: Found no space between name and residues in $infile, line $.: '$line'\n");
	    } else {
		$name=$1;
		$residues=$2;
		$residues=~tr/ //d;
		$residues_per_line=length($residues);
	    }
	    if ($block==1) {
		$names[$k]=$name;
		$seqs[$k]=$residues;
	    } else {
		$seqs[$k].=$residues;
		if ($names[$k] ne $name) {
		    print("WARNING: name of sequence $k in block 1 ($names[$k]) is not the same as in block $block ($name) in $infile\n");
		}
	    }
#	    printf("%3i %3i %-16.16s %s\n",$block,$k,$names[$k],$seqs[$k]);
	    $k++;
	}
	if ($k && $n && $n!=$k) {die("\nError: different number of sequences in blocks 1 and $block of $infile\n");} 
	if (!$n) {$n=$k;}
    }

   # Read psi format
    elsif ($informat eq "psi")
    {
	my $name;
	my $first_block=1;
	
	$n=0;
	while ($line = <INFILE>) #scan through PsiBlast-output line by line
	{
	    if ($line=~/^\s*$/) {next;}
	    if ($line!~/\s*(\S+)\s+(\S+)/) {
		die ("\nERROR found in psi format: $!");
	    }
	    if ($1 eq "aa_dssp" || $1 eq "ss_dssp" || $1 eq "aa_pred" || $1 eq "ss_pred" || $1 eq "ss_conf") {next;}
	    $names[$n]=$1;
	    $seqs[$n]=$2;
	    $n++;
	}
    }
    
    
    close INFILE;
    
    
    # Empty input file?
    if ($n==0) {die("\nERROR: input file $infile contains no sequences\n");}
 

################################################################################################
# Removing non-alphanumberic symbols
################################################################################################
    for ($k=0; $k<$n; $k++) {
	$seqs[$k]=~tr/A-Za-z0-9.-//cd;
    }
    
    
################################################################################################
# Filling up with gaps '.' or deleting gaps
################################################################################################
    
    if ($informat eq "a3m" && !$remove_inserts)
    {
	my @len_ins;   # $len_ins[$j] will count the maximum number of inserted residues after match state $i.
	my $j;       # counts match states
	my @inserts; # $inserts[$j] contains the insert (in small case) of sequence $i after the $j'th match state
	my $insert;
	
	# Determine length of longest insert 
	for ($k=0; $k<$n; $k++)
	{
	    # split into list of single match states and variable-length inserts
	    # ([A-Z]|-) is the split pattern. The parenthesis indicate that split patterns are to be included as list elements
	    # The '#' symbol is prepended to get rid of a perl bug in split
	    @inserts = split(/([A-Z]|-|~|[0-9])/,"#".$seqs[$k]."#"); 
	    $j=0;
#	printf("Sequence $k: $seqs[$k]\n");
#	printf("Sequence $k: @inserts\n");
	    
	    # Does sequence $k contain insert after match state j that is longer than $ngap[$j]?
	    foreach $insert (@inserts) 
	    {
		if( !defined $len_ins[$j] || length($insert)>$len_ins[$j]) {$len_ins[$j]=length($insert);}
		$j++;
#	    printf("$insert|");
	    }
#	printf("\n");
	}
	my $ngap;
	
	# Fill up inserts with gaps '.' up to length $len_ins[$j]
	for ($k=0; $k<$n; $k++)
	{
	    # split into list of single match states and variable-length inserts
	    @inserts = split(/([A-Z]|-|~|[0-9])/,"#".$seqs[$k]."#");
	    $j=0;
	    
	    # append the missing number of gaps after each match state
	    foreach $insert (@inserts) 
	    {
		for (my $l=length($insert); $l<$len_ins[$j]; $l++) {$insert.=".";}
		$j++;
	    }
	    $seqs[$k] = join("",@inserts);
	    $seqs[$k] =~ tr/\#//d; # remove the '#' symbols inserted at the beginning and end
	}
    }    
    
    
################################################################################################
# Match state assignment
################################################################################################

    # Use gap rule for match state assignment?
    if ($matchmode eq "gaprule") {

	my @gaps=();
	my $residues;
	my @residues;
	
	# Determine number of gaps per column
	for ($k=0; $k<$n; $k++) {
	    @residues=unpack("C*",$seqs[$k]);
	    for (my $l=0; $l<@residues; $l++) {
		if ($residues[$l]==46 || $residues[$l]==45) {
		    if (defined $gaps[$l]) {$gaps[$l]++;} else {$gaps[$l]=1;}
		}
	    }
	} 
	
	# Set columns with more gaps than $match_gaprule to lower case,
	for ($k=0; $k<$n; $k++) {
	    @residues=unpack("C*",$seqs[$k]);
	    $residues="";
	    for (my $l=0; $l<@residues; $l++) {
		if (!defined $gaps[$l] || $gaps[$l]<0.01*$match_gaprule*$n) {
		    if ($residues[$l]==46) {
			$residues .= "-";
		    } else {
			$residues .= uc(chr($residues[$l]));
		    }
		} else {
		    if ($residues[$l]==45) {
			$residues .= ".";
		    } else {
			$residues .= lc(chr($residues[$l]));
		    }
		}
		$seqs[$k]=$residues;
	    }
	}
    }
	    
    # Use first sequence for match state assignment?
    if ($matchmode eq "first") {
	
	my @match=();
	my $residues;
	my @residues;
	
	# Determine which columns have a gap in first sequence
	@residues=unpack("C*",$seqs[0]);
	for (my $l=0; $l<@residues; $l++) {
	    if ($residues[$l]==46 || $residues[$l]==45) {$match[$l]=0;} else {$match[$l]=1;}
	}

	# Set columns without residue in first sequence to upper case,
	for ($k=0; $k<$n; $k++) {
	    @residues=unpack("C*",$seqs[$k]);
	    $residues="";
	    for (my $l=0; $l<@residues; $l++) {
		if ($match[$l]) {
		    if ($residues[$l]==46) {
			$residues .= "-";
		    } else {
			$residues .= uc(chr($residues[$l]));
		    }
		} else {
		    if ($residues[$l]==45) {
			$residues .= ".";
		    } else {
			$residues .= lc(chr($residues[$l]));
		    }
		}
		$seqs[$k]=$residues;
	    }
	}
    }
    

################################################################################################
# Remove gaps etc.
################################################################################################

    # Remove columns with too many gaps?
    if ($remove_gapped) {

	my @gaps=();
	my $residues;
	my @residues;
	
	# Determine number of gaps '.' or '-' per column
	for ($k=0; $k<$n; $k++) {
	    @residues=unpack("C*",$seqs[$k]);
	    for (my $l=0; $l<@residues; $l++) {
		if ($residues[$l]==45 || $residues[$l]==46) {
		    if (defined $gaps[$l]) {$gaps[$l]++;} else {$gaps[$l]=1;}
		}
	    }
	} 
	
	# Remove columns with too many gaps
	for ($k=0; $k<$n; $k++) {
	    @residues=unpack("C*",$seqs[$k]);
	    $residues="";
	    for (my $l=0; $l<@residues; $l++) {
		if (!defined $gaps[$l] || $gaps[$l]<0.01*$remove_gapped*$n) {
		    $residues .= chr($residues[$l])
		}
		$seqs[$k]=$residues;
	    }
	}
    }
	    
    # Remove/transliterate gaps? 
    for ($k=0; $k<$n; $k++) {$seqs[$k]=~tr/ //d;}
    if ($remove_inserts) {
	for ($k=0; $k<$n; $k++) {
	    $seqs[$k]=~tr/a-z.//d;
#	    printf("%s\n",$seqs[$k]);
	}
    }	
    

    # Remove sequences which contain only gaps
    my $nin=$n;
    for ($k=0; $k<$n; $k++) {
	if (($seqs[$k]=~tr/a-zA-Z/a-zA-Z/==0)) {
	    if ($v>=2) {print("Sequence contains only gaps and is removed: $names[$k]\n");}
	    splice(@seqs,$k,1);
	    splice(@names,$k,1);
	    $k--; $n--;
	}
    }
    

    # Cut down sequence names to $DESCLEN characters maximum
    for ($k=0; $k<$n; $k++) {$names[$k]=substr($names[$k],0,$DESCLEN);}

    if ($outformat eq "a3m") {
	# Remove all '.' gaps
	for ($k=0; $k<$n; $k++) {$seqs[$k]=~tr/.//d;}
    } elsif ($outformat eq "fas" || $outformat eq "clu" || $outformat eq "sto" || $outformat eq "psi" ) {
	# Substitute all '.' by '-'
	for ($k=0; $k<$n; $k++) {$seqs[$k]=~tr/./-/;}
    }
    if ($gap ne "default") {
	for ($k=0; $k<$n; $k++) {$seqs[$k]=~s/\./$gap/g; $seqs[$k]=~s/-/$gap/g;}
    }
    if ($case eq "uc") {
	for ($k=0; $k<$n; $k++) {$seqs[$k]=~tr/a-z/A-Z/;}
    } elsif ($case eq "lc") {
	for ($k=0; $k<$n; $k++) {$seqs[$k]=~tr/A-Z/a-z/;}
    }
    
    

    ####################################################
    # Check that sequences have same length
    ####################################################
    if ($outformat ne "a3m") {
	my $len=length($seqs[0]);
	for($k=1; $k<$n; $k++) {
  	    if (length($seqs[$k])!=$len) {
		printf("\nError: Sequences in $infile do not all have same length, e.g. %s (len=%i) and %s (len=%i)\n",
		       $names[0],$len,$names[$k],length($seqs[$k]));
		exit 1;
	    }
	}
    }



################################################################################################
# Output part
################################################################################################
    
    if ($num) {
	for ($k=0; $k<$n; $k++) {$names[$k]=~s/^(\d+\|)?(.*)/$k|$2/;}
    }

    #Write sequences into output file
    open (OUTFILE, ">$outfile") or die ("cannot open $outfile:$!\n");
    if ($outformat eq "sto" || $outformat eq "psi") {
	for ($k=0; $k<$n; $k++) {
	    $names[$k] =~ s/\s*(\S+).*/$1/;
	    printf(OUTFILE "%-32.32s %s\n",$names[$k],$seqs[$k]);
	}
	
    } elsif ($outformat eq "clu") {
	printf(OUTFILE "CLUSTAL\n\n\n");
	while($seqs[0] ne "") {                         # While there are still residues left
	    for ($k=0; $k<scalar(@names); $k++) {  # go through all sequences
		$names[$k] =~ s/\s*(\S+).*/$1/;
		$seqs[$k]=~s/(\S{1,$numres})//;           # remove the leftmost (up to) 60 residues from sequence $nseq
		printf(OUTFILE "%-15.15s  %s\n",$names[$k],$1); # and print them after the sequence name
	    }
	    print(OUTFILE "\n");                            # leave a blank line after each block of 60 residues
	}
    } else {
	for ($k=0; $k<$n; $k++) {
	    $seqs[$k]=~s/(\S{$numres})/$1\n/g;
	    printf(OUTFILE ">%s\n%s\n",$names[$k],$seqs[$k]);
	}
    }
    
    close OUTFILE;
    if ($v>=2) {
	if ($nin==1) {print("Reformatted $infile with 1 sequence from $informat to $outformat and written to file $outfile\n");}
	else {
	    if (!$nin==$n) {printf("Removed %i sequences which contained no residues\n",$nin-$n); }
	    print("Reformatted $infile with $n sequences from $informat to $outformat and written to file $outfile\n");
	}
    }
    
# return number of hits in one byte (0-255)
if    ($n<110)  {return $n;}
elsif ($n<1100) {return 100+int($n/10);}
elsif ($n<5500) {return 200+int($n/100);}
else                {return 255;}              

}
