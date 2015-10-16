#!/usr/local/bin/perl -w

=head1 NAME

    logo.pm - organizes data in FASTA and CLUSTAL formats into height data.

=head1 SYNOPSIS

    Perl module

=head1 DESCRIPTION

    logo.pm: Takes in strings of aligned sequences and sorts them vertically
             based on height as assigned by the following equations found in
             Schneider and Stephens paper "Sequence Logos: A New Way to Display
             Consensus Sequences":
    
                 height = f(b,l) * R(l)                        (1)

             where f(b,l) is the frequency of base or amino acid "b" at position
             "l". R(l) is amount of information present at position "l" and can
             be quantified as follows:

                 R(l) for amino acids = 2 - (H(l) + e(n))      (2a)
                 R(l) for bases = log(20) - (H(l) + e(n))      (2b)

             where log is taken base 2, H(l) is the uncertainty at position "l", 
             and e(n) is the error correction factor for small "n". H(l) is
             computed as follows:

                 H(l) = - (Sum f(b,l) * log[ f(b,l) ])         (3)

             where again, log is taken base 2. f(b,l) is the frequency of base
             "b" at position "l". The sum is taken over all amino acids or
             bases, depending on which the data is.

             Currently, logo.pm uses an approximation for e(n), given by:

                 e(n) = (s-1) / (2 * ln2 * n)                  (4)

             Where s is 4 for nucleotides, 20 for amino acids ; n is the number
             of sequences in the alignment. e(n) also  gives the height of error
             bars.

=cut

package logo;

use strict;

################################################################################
######                             SOME VARIABLES                         ######
################################################################################

my $AA = 0;
my $NA = 1;

my %BASES = ("a" => "adenine",
	     "t" => "thymine",
	     "g" => "guanine",
	     "c" => "cytosine",
	     "u" => "uracil");

# does not include B or Z
my %AMINOACIDS = ("a" => "", "c" => "", "d" => "", "e" => "", "f" => "",
		  "g" => "", "h" => "", "i" => "", "k" => "", "l" => "",
		  "m" => "", "n" => "", "p" => "", "q" => "", "r" => "",
		  "s" => "", "t" => "", "v" => "", "w" => "", "y" => "");

my @data;
my $kind;
my ($seqs_r, $desc_r);


################################################################################
######                             SOME FUNCTIONS                         ######
################################################################################

=head1 APPENDIX

=cut

=head2 getHeightData()

 Usage   : my ($height_data_r, $description_data_r, $kind) =
              logo::getHeightData($input_data_r, $params);
 Returns : * REFERENCE TO array of height data
           * REFERENCE TO array of input description strings
           * $AA if the data is amino acid, $NA otherise
 Args    : $input_data_r : input data in CLUSTAL or FASTA formats
         : $params       : hash of parameters

 getHeightData is the entry point into the logo.pm module. $input_data_r is a
 reference to  an array of strings containing FASTA or CLUSTAL data, where all
 lines whose first character is "#" is considered a comment line.

 $params is a hash of parameters with the following keys:
   * smallsampletoggle : 0 to turn off small sample correction, otherwise
                         small sample correction is on
   * input_kind : 0 for amino acids, 1 for nucleic acids; if undefined,
                  logo.pm will attempt to automatically detect whether the
                  input consists of amino or nucleic acid data. If
                  $input_kind is defined, only those residues  defined by
                  $input_kind will be in the output -- all other residues are
                  considered as spaces. For example, if $input_kind is $NA,
                  the residue "I" or "i" are considered spaces, since "I" and
                  "i" are not nucleic acid residues.
   * stretch : stretch all characters so they are flush at the maximum number
               of bits allowed

 Sample use:

  # get FASTA data
  open (FASTA, "$fastafile");
  my @inputdata = <FASTA>;
  close (FASTA);

   my %heightparams = (
		       smallsamplecorrection => 0,
		       input_kind => 0,
		       stretch => 0
		       );

  # get height data
  my ($heightdata_r, $desc_r, $kind) = logo::getHeightData(\@inputdata, \%heightparams);

=cut

# entry point into module
sub getHeightData {

    # $smallsampletoggle is toggle to turn off small sample correction (1 to turn off)
    # $input_kind can be $AA or $NA or undef
    my ($input_data_r, $params) = @_;

    # gets sequences, sets $kind temporarily
    my ($goodlength, $maxreslength, $badline);
    ($seqs_r, $desc_r, $maxreslength, $goodlength, $badline) = getSeqs(@$input_data_r, $params->{input_kind});

    # check for bad length
    if (!$goodlength) {
	return (undef, undef, undef, $goodlength, $badline);
    }

    # reset $kind if in $input_kind
    if (defined $params->{input_kind} && isLegalKind($params->{input_kind}) ) {
	$kind = $params->{input_kind};
    }

    # build data
    buildData(@$seqs_r, $params->{smallsampletoggle}, $params->{stretch}, $maxreslength);

    return (\@data, $desc_r, $kind, $goodlength, $badline);
}

sub isLegalKind {
    return ($_[0] =~ /^[01]$/);
}

sub getSeqs {
    # skip all comment chars and lines of all spaces
    while ( ($_[0] =~ /^\s*\#/) || ($_[0] =~ /^\s*$/) ) {
	shift @_;
    }

    if ($_[0] =~ />/) {
	return getSeqs_FASTA(@_);
    } else {
	return getSeqs_CLUSTAL(@_);
    }
}

################################################################################
######                          FORMATTING FUNCTIONS                      ######
################################################################################

sub getSeqs_CLUSTAL {
    my @returnVal;
    my @desc;
    my $seqCount=0;
    my $reslength = 0;
    my ($name, $seq);

    my $input_kind = pop @_;
    my $confidence_limit = 0.90;
    my $NA_count = 0;
    my $total_residues = 0;
    my ($prevlinelength, $linelength) = (0,0);

    foreach (@_) {
	chomp;

	# skip if it is a comment character -- first character is "#"
	next if (/^\s*\#/);

	# if spaces or just "*" and "."
	if (/^[\*\.\s]*$/) {
	    $seqCount=0;
	    $prevlinelength=0;
	    next;
	}

	($name,$seq) = (/^\s*(\S+)\s+(\S+)\s*$/);

	# add new entry
	if (!defined $desc[$seqCount]) {
	    $desc[$seqCount] = $name;
	    $returnVal[$seqCount] = "";
	}

	my @chars = split(//,$seq);
	my $char;
	foreach (@chars) {
	    if ($seqCount == 0) {
		$reslength++;     # all sequences have same residue length, so only count first one
	    }

	    $total_residues++;
	    $linelength++;

	    # set $char
	    if (defined $input_kind) {
		if ($input_kind == $AA) {
		    $char = (isAA($_)) ? $_ : "-";
		} else { # == $NA
		    $char = (isNA($_)) ? $_ : "-";
		}
	    } else {
		$char = $_;
		if (isNA($char)) {
		    $NA_count++;
		}
	    }
	    
	    $returnVal[$seqCount] .= $char;
	}

	if ($seqCount == 0) {
	    $prevlinelength = $linelength;
	} elsif ($prevlinelength != $linelength) {  # different number of residues, so complain
	    return (undef, undef, undef, 0, $name);
	}
	$linelength=0;

	$seqCount++;
    }

    # determine whether to use $NA or $AA
    if (!defined $input_kind) {
	if ($NA_count / $total_residues >= $confidence_limit) {
	    $kind = $NA;
	} else { 
	    $kind = $AA;
	}
    }

    return (\@returnVal, \@desc, $reslength, 1, undef);
}

# if $input_kind is defined, residues that are not defined are set to space
sub getSeqs_FASTA {
    my @returnVal;
    my @desc;
    my $count=-1;
    my $newElem=0;

    my $input_kind = pop @_;

    my $confidence_limit = 0.90;
    my $NA_count = 0;
    my $total_residues = 0;
    my $reslength = 0;
    my $maxreslength = 0;
    
    my ($goodlength, $currline, $prevline);


#    # skip all lines that are all spaces
#    while ($_[0] =~ /^\s*$/) {
#	shift @_;
#    }

    foreach (@_) {

	# skip if it is a comment character -- first character is "#"
	next if (/^\s*\#/);

	# skip all lines that are all spaces
	next if (/^\s*$/);

	$_ =~ s/\s+$//;  # cut trailing white space
	$_ =~ s/^\s+//;  # cut leading white space
	if (/>/) {
	    $currline = $_;
	    ($desc[scalar @desc]) = ($_ =~ />\s*(.+)$/);

	    if (not $newElem) {		
		$count++;
		$newElem = 1;
	    }
	} else {
	    if ($newElem) {
		$maxreslength = $reslength if $maxreslength == 0;
		if (($maxreslength != 0) && ($maxreslength != $reslength)) {
		    return (undef, undef, undef, 0, $prevline);
		}

		$maxreslength = $reslength;
		$reslength = 0;
	    }

	    my @chars = split(//,$_);
	    my $char;
	    foreach (@chars) {
		$reslength++;
		$total_residues++;

		# set $char
		if (defined $input_kind) {
		    if ($input_kind == $AA) {
			$char = (isAA($_)) ? $_ : "-";
		    } else { # == $NA
			$char = (isNA($_)) ? $_ : "-";
		    }
		} else {
		    $char = ($_ =~ /[a-zA-Z]/) ? $_ : "-";  # if not alpha char, use space
		    if (isNA($char) && !isSpace($char)) {
			$NA_count++;
		    }
		}

		if ($newElem) {
		    $returnVal[$count] = $char;
		} else {
		    $returnVal[$count] .= $char;
		}
		$newElem = 0;
	    }
	    $prevline = $currline if $currline =~ />/;
	}
    }

    # check if last is biggest
    if (($maxreslength != 0) && ($maxreslength != $reslength)) {
	return (undef, undef, undef, 0, $prevline);
    }
#    $maxreslength = ($reslength > $maxreslength) ? $reslength : $maxreslength;

    # determine whether to use $NA or $AA
    my $ratio = $NA_count / $total_residues;
    if (!defined $input_kind) {
	if ($NA_count / $total_residues >= $confidence_limit) {
	    $kind = $NA;
	} else { 
	    $kind = $AA;
	}
    }

    return (\@returnVal, \@desc, $maxreslength || $reslength, 1, undef);  # 1 for good lengths
}

sub isSpace {
    return $_[0] =~ /[ \-]/;
}

sub isAA {
    return (defined $AMINOACIDS{lc $_[0]});
}

sub isNA {
    return (defined $BASES{lc $_[0]});
}

################################################################################
######                       DATA BUILDING FUNCTIONS                      ######
################################################################################


# arguments: takes reference to array and lines aligned sequences of bases or
#            amino acids
# returns: updated reference array to reflect contents of sequences sorted
#          vertically by height as described by (1)
sub buildData {
    
    my $currentx = 0;
    my $h;
    my $count=0;
    my $maxreslength = pop (@_);
    my $stretch = pop(@_);
    my $smallsampletoggle = pop (@_);
    my $totalsize = $#_+1;

    while ($currentx < $maxreslength) {       #(length $_[0])) {
	my $allspaces = 1;
	my $spaceCount=0;

	# get vertical sequence
	my @vert=();
	foreach (@_) {  # foreach sequence
	    my $currentchar;

	    # set currentchar, set to " " if $_ is not long enough
	    if ($currentx >= (length $_)) {
		$currentchar = " ";
	    } else {
		$currentchar = substr($_,$currentx,1);
	    }

	    # in all sequences, "-" is considered as a space
	    # don't count " " and "-"
	    if (($currentchar ne "-") && ($currentchar ne " ")) {
		$vert[scalar @vert] = uc substr($_,$currentx,1);
		$allspaces = 0;
	    } else {
		$spaceCount++;
	    }
	}

	if ($allspaces) {
	    # build @vert
	    @vert = (" 0", ">0");

	    # place in @data
	    $data[scalar @data] = \@vert;

	    $currentx++;
	    next;
	}

	my $error;
	if ($smallsampletoggle) {
	    $error=getError($kind,$totalsize - $spaceCount);
	} else {
	    $error = 0;
	}

	# sort vertical sequence by amino acid or base
	@vert = sort(@vert);
	my $total = $#vert + 1;

	# find H(l) -- must be done before collapsing
	$h = getH(@vert);

	# collect like terms
	@vert = collapse(@vert);

	# get R(l)
	my $r;
	if (!defined $stretch || !$stretch) {
	    $r= getR($kind, $h, $error);
	} else {
	    $r = ($kind == $NA) ? 2 : (log(20) / log(2));
	}

	# place heights
	my $count=0;
	my $height;
	my $elem;
	foreach $elem (@vert) {
	    $height = getHeight(substr($elem, 2) / $total,$r);
	    $vert[$count] = substr($elem,0,1) . $height;
	    $count++;
	}

	# sort by height
	@vert = sort height_sort @vert;

	# put in error bar size
	$vert[$count] = ">$error";

	# place in @data
	$data[scalar @data] = \@vert;

	$currentx++;
    }
}

# uses error approximation given by:
#             e :=  (s-1) / (2 * ln2 * ntrue);
sub getError {
    return ((($_[0] == $NA) ? 4 : 20) - 1) / (2 * log(2) * $_[1]);
}

sub height_sort {
    my ($lettera, $heighta) = ($a =~ /^(.{1})(\S+)$/); #substr($a,1);
    my ($letterb, $heightb) = ($b =~ /^(.{1})(\S+)$/); #substr($b,1);
    
    # compare by height first, then letter
    if ($heighta <=> $heightb) {
	return $heighta <=> $heightb;
    } else {
	return $letterb cmp $lettera;  #reversed for some reason...
    }
}

sub collapse {
    my @returnVal;
    my $current = $_[0];
    my $count=0;
    my $freq;

    foreach (@_) {
        if ($current eq $_) {
            $count++;
        } else {
	    $returnVal[scalar @returnVal] = "$current $count";

            $current = $_;
            $count=1;
        }
    }

    # add last element
    $returnVal[scalar @returnVal] = "$current $count";

    return @returnVal;
}

# arguments : $_[0] : list of bases or amino acids
sub getH {
    my $h = 0;
    my (@vert) = @_;  # vertical sequence (comparing multiple sequences)

    my $current = $vert[0];
    my $count=0;
    my $freq;

    foreach (@vert) {
	if ($current eq $_) {
	    $count++;
	} else {
	    $freq = $count / ($#vert + 1);
	    $h += $freq * (log ($freq) / log(2));

	    $current = $_;
	    $count=1;
	}
    }

    # add last element
    $freq = $count / ($#vert + 1);
    $h += $freq * (log ($freq) / log(2));

    return -$h;
}

# arguments : $_[0] : $AA or $NA
#             $_[1] : uncertainty = H(l)
#             $_[2] : error correction factor for small number of sequences
#                     = e(n) ; assummed to be 0 if not given
sub getR {
    my $max = ($_[0] == $NA) ? 2 : (log(20) / log(2));
    my $e = (defined $_[2]) ? $_[2] : 0;
    return ($max - ($_[1] + $e));
}

# arguments: $_[0] : frequency = f(b,l)
#            $_[1] : R(l)
sub getHeight {
    return $_[0] * $_[1]
}

1;
