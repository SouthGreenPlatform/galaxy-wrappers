#!/usr/bin/perl -w
use strict;

###PCOILS PROGRAM###

#run PCoils not only over the sequence, but over all sequences in the corresponding alignment (as to compare it with Coils)
#only authorized regions according to the ones in the input file are considered
#only windows without gaps are considered

#THIS PROGRAM NEEDS THE EXECUTABLE NCOILS_PLAY (based on ncoils_play.c and read_matrix_play.c) or some comparable program
#ENVIRONMENT VARIABLE COILSDIR NEEDS TO BE SET

#output: coils roc curve or sensitivity/specificity curve
#input: preferably a file in fasta format

#OPTIONS:
my $window_size=21; 
print (STDERR "PCoils window size is $window_size\n");

unless(defined($ARGV[0])) {die "usage: perl PCoils.pl <input_file> <database>\nThe database is ideally the NR90 filtered by Pfilt with the -c option (leaving the coiled-coil regions unfiltered)\n";}
unless(defined($ARGV[1])) {die "usage: perl PCoils.pl <input_file> <database>\nThe database is ideally the NR90 filtered by Pfilt with the -c option (leaving the coiled-coil regions unfiltered)\n";}

my $extension;
my $name_counter;
my $plus_counter=0;
my $minus_counter=0;
my $seq_long_enough;
my $name;
my $line;
my $beginning;
my $end;
my $chain;
my $complete_fasta_cc;
my $length_complete_fasta_cc;
my $helical_part_of_socket_sequence;
my $ignore_for_helix;
my $ignore_for_reg;
my $starting_point;
my @register;
my @complete_sequence;
my $limit;
my $svm_counter;
my $svm_line_counter;
my $i;
my $print_cc_out;
my $low;
my $high;
my $play_line;
my $svm_line;
my $threshold;
my $false_positives;
my $true_positives;
my $positives;
my $negatives;
my $j;
my $k;
my $pre_printer;
my $post_printer;
my @authorised;
my $play_counter;
my $last;
my $intermediate_line;
my $log_extend_line;
my $position_i;
my $position_counter=0;
my $chomped_line;
my $length_of_sequence;
my $cc_line;
my @sequence;
my @authorization;
my $field_counter;
my $valid_window;
my $window;
my %window_hash;
my %sequence_hash;
my @problematic_lines;
my $again_svm_line;
my $again_svm_line_counter;
my $length_of_line;
my $number_of_problems=0;
my $intermediate_counter=0;
my $helpbuffer_line;
my $helpmulticoilline;
my $hbmline;
my $okay_flag=0;
my $win_line;
my $centralresiduenumber;
my $coils_index;
my $negative_counter=0;
my $coilsline;
my $coils_counter;
my @coils_output;
my $counter=0;
my $hashline;
my $gi_number;
my $cpo_line;
my $true_negatives;
my $false_negatives;
my $sensitivity;
my $specificity;
my $classification;
my $code;
my $skip_this_residue;
my $weighting;

$extension=($window_size-1)/2;
$centralresiduenumber=$extension+1;
$coils_index=$window_size+1;

sub profile_for_roc()
{
    #INPUT: alignment in PSI/FAS-BLAST format
    #OUTPUT: PCoils score (same format as run_play_full_without_window_maximization)
    
    my $m;
    my $n;
    my $coils_line;
    my $residue_counter=0;
    my $sequence_counter=0;
    my @coils_output_al;
    my $max_res=0;
    my @average_score;
    my $sequence_line;
    my @seq;
    my @gap;
    my $first_sequence=1;
    my $valid_sequence_counter;

    system("perl manipulate_pdbaa.pl temp_files/hhfilter_output.fas > temp_files/hhfilter_output.fas_cleaned"); #NEW
    system("mv temp_files/hhfilter_output.fas_cleaned temp_files/alignment.fas"); #NEW

    $sequence_counter=0;
    #read in alignment
    open(A,"temp_files/alignment.fas") || die "error opening\n";
    while($line=<A>)
    {
	#this is always the start of a new sequence
	if($line=~/>/)
	{
	    #put gi number and sequence in a buffer file
	    open(BUF,">temp_files/buffer_file") || die "error opening\n";
	    #print gi number in file
	    print BUF $line;
	    #get sequence
	    $line=<A>;
	    #this is needed for chomping without altering the original variable
	    $sequence_line=$line;
	    #get rid of newline
	    chomp $sequence_line;
	    #make sure arrays are empty
	    undef @seq;
	    #put sequence in array
	    @seq=split(//,$sequence_line);	    
	    
	    $line=~s/-/X/g; #gaps are transformed into neutral residues
	    #print sequence in file
	    print BUF $line;
	    close(BUF);

	    #run Coils over the sequence in the buffer file
	    if($sequence_counter==0) #only run PCoils once (for the first sequence in the alignment)
	    {
		system("./PCoils_core -nw -win $window_size -prof temp_files/myhmmmake_PCoils_output < temp_files/buffer_file > PCoils.output");
	    }
	}
	$sequence_counter++;
    }
    close(A);
}

#MAIN part of program

system("./blast/blastpgp -d $ARGV[1] -i $ARGV[0] -o temp_files/alignhits_input -e 1E-4 -j 2");
system("perl alignhits.pl -psi -b 1.0 -e 1E-4 -q $ARGV[0] temp_files/alignhits_input temp_files/alignhits_output.psi");
system("perl reformat.pl -uc temp_files/alignhits_output.psi temp_files/hhfilter_input.a3m");
system("./hhfilter -i temp_files/hhfilter_input.a3m -qid 40 -cov 20 -o temp_files/hhfilter_output.a3m");
system("perl reformat.pl -uc temp_files/hhfilter_output.a3m temp_files/hhfilter_output.fas");
system("perl myhmmmake_PCoils.pl"); #creates ~/write/myhmmmake_PCoils_output that is used by PCoils_core

profile_for_roc();
