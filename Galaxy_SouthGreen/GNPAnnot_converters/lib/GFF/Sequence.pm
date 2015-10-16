
#FUNCTION  GFF subsequence extraction module

=head1 NAME

    GFF::Sequence.pm - GFF subsequence extraction module

=head2 AUTHORSHIP

    Copyright (c) 2000  Richard Bruskiewich (rbsk@sanger.ac.uk)
    Sanger Centre, Wellcome Trust Genome Campus, Cambs, UK
    All rights reserved.

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation

=head2 SYNOPSIS

    use GFF::Sequence ;

=head2 DESCRIPTION

 This module provides access to the real, corresponding 
 subsequences of GFF Perl module objects (see GFF.pm modules at
 http://www.sanger.ac.uk/Software/GFF/GFF.shtml)

=cut

package GFF::Sequence;

use Carp;
use strict;
use Cwd ;

# I assume that these modules are available locally?
use GFF ;
use Bio::SeqIO ;
use Bio::PrimarySeq ;

require Exporter ;
use vars qw($AUTOLOAD @ISA @EXPORT @EXPORT_OK);

use vars qw() ;

#
# @ISA has our inheritance.
#
@ISA = qw(Exporter GFF) ;
#
# Place functions/variables you want to *export*, ie be visible from the caller package
#
@EXPORT = qw() ;
@EXPORT_OK = qw() ;

=head2 METHODS

 The general protocol is to:

 1. Open the sequence file associated with the GFF object
 2. Retrieve sequence components defined by the GFF object
 3. Close the sequence file

=over 4

=item open( -seqname=>$seqname, -file=>$file, [-format=>$format, origin=>$origin] )

This method opens a '$seqname', in file '$file' of a given '$format'.

The relative offset between sequence and GFF is specified by the
$origin argument, which is the GFF position corresponding to 
base position 1 in the sequence (default: 1).

Either one or both of $seqname and $file must/can be given. 

If only the $seqname is given then the method attempts 
to open a file '$seqname.seq' in the current working directory. 

If only the $file is given, the $seqname name is set to whatever 
sequence name is defined in the file.

The method throws an exception if the file cannot be opened.

The '$format' defaults to 'Fasta', but may be optionally set to
any of those available to the version of Bioperl being used.

=back

=cut

sub next ;

sub open() {
    my $ref     = shift;
    my $class   = ref($ref) || $ref ;

    my %args=(
	      '-format' => 'Fasta',
	      '-origin' => 1,
	      @_,
	      ) ;

    my $seqname = $args{'-seqname'}; 
    my $file    = $args{'-file'};

    confess "GFF::Sequence::open(): needs either a -sequence or a -file argument!\n" 
	unless(defined($seqname) or defined($file)) ;

    $file = (Cwd::cwd)."/$seqname.seq" if !$file ;
    my $format = $args{'-format'};
    my $origin = $args{'-origin'};

    my $self = {} ;
    eval {
       $self->{'_SEQIOH'} = Bio::SeqIO->new( -file => "$file", '-format' => $format );
    } ;
    if($@) { 
	confess "GFF::Sequence::open() ERROR: I could not read in sequence from file '$file': $@\n" ;
    }
    
    bless $self, $class;

    $self->next($seqname,$origin) ;

    return $self;
}

=pod

=over 4

=item next([-seqname => $seqname, -origin => $origin])

If the sequence file previously opened contains multiple sequences, 
then invocation of this method advances to the next sequence in the file.

If the optional '$seqname' value is given, then the given sequence
is given that name, overriding any display id given for the current
sequence read out of the file.

The relative offset between sequence and GFF is specified by the
optional $origin argument, which, if defined, is the GFF position
corresponding to base position 1 in the sequence (default: 1).

The method returns the Bio::Seq object reference of the sequence read in;
otherwise, a null value when no further sequences are available.

=back

=cut

sub next {
    my $self = shift ;
    my %args=(
	      '-origin' => 1,
	      @_,
	      ) ;

    my $seqname = $args{'-seqname'}; 
    my $origin  = $args{'-origin'};

    confess "GFF::Sequence::next() ERROR: no sequence file is open to read?"
	unless(defined($self->{'_SEQIOH'})) ;
 
    $self->{'_SEQOBJ'} = $self->{'_SEQIOH'}->next_seq() ;

    if($self->{'_SEQOBJ'}) {
	
	if($seqname) {
	    $self->{'_SEQNAME'} = $seqname ;
	} else {
	    $self->{'_SEQNAME'} = $self->{'_SEQOBJ'}->display_id() ;
	}
	$self->{'_OFFSET'} = $origin-1 ;

    } else {
	undef $self->{'_SEQOBJ'} ;
	undef $self->{'_SEQNAME'} ;
	undef $self->{'_OFFSET'} ;
    }

    return ($self->{'_SEQOBJ'}) ;
}

=pod

=over 4

=item subseq($gf)

This method returns a string corresponding to the 
nucleotide subsequence of the currently loaded sequence (see open/next above) 
corresponding to the GFF feature object given as its argument.

=back

=cut

sub subseq {
    my $self = shift ;
    my $gf   = shift ;

    confess "GFF::Sequence::subseq() ERROR: sequence object undefined?\n"
	unless(defined($self->{'_SEQOBJ'})) ;

    GFF::GeneFeature->verify($gf) ;

    my $offset = $self->{'_OFFSET'} ;

    my $subseq = $self->{'_SEQOBJ'}->subseq($gf->start-$offset,$gf->end-$offset) ;
    my $gfseq  = Bio::PrimarySeq->new( -seq => $subseq,
				       -moltype => 'dna', ) ;
    if( $gf->strand eq '-') {
	$gfseq = $gfseq->revcom ;
    }
    return($gfseq->seq) ;
}

=pod

=over 4

=item translation($gf)

This method returns the amino acid translation of the 
subsequence of the currently loaded sequence (see open/next above)
corresponding to the GFF feature object given as its argument.

Note that this method expects that the strand and frame (phase)
of the GFF record are properly set.  Note that one simple exon
boundary GFF, the frame may indicate that the first one or two
nucleotides are ignored in the translation, and the last codon
may also be so undefined.

=back

=cut

sub translation {
    my $self = shift ;
    my $gf   = shift ;

    confess 'GFF::Sequence::translation() ERROR: sequence object undefined?\n'
	unless(defined($self->{'_SEQOBJ'})) ;

    GFF::GeneFeature->verify($gf) ;

    my ($start,$end,$strand,$frame) ;

    if($gf->strand !~ /^[\+-]$/) {
	carp "GFF::Sequence::translation() WARNING: undefined or invalid feature strand? Assuming '+'\n" ;
	$strand = '+' ;
    } else {
	$strand = $gf->strand ;
    }

    if($gf->frame eq '.') {
	carp "GFF::Sequence::translation() WARNING: undefined feature frame/phase? Assuming '0'\n" ;
	$frame = 0 ;
    } else {
	$frame = $gf->frame+0 ;
    }

    if( $strand eq '+') {
	$start = $gf->start+$frame ;
	$end   = $gf->end ;
    } else { # '-'
	$start = $gf->start ;
	$end   = $gf->end-$frame ;
    }

    my $offset = $self->{'_OFFSET'} ;

    my $subseq = $self->{'_SEQOBJ'}->subseq($start-$offset,$end-$offset) ;
    my $dna  = Bio::PrimarySeq->new( -seq => $subseq,
				     -moltype => 'dna', ) ;
    if( $gf->strand eq '-') {
	$dna = $dna->revcom ;
    }
    my $peptide = $dna->translate() ;
    return($peptide->seq) ;
}

=pod

=over 4

=item codons($gf,$position,$alleles)

A specialized version of 'subseq', this method returns an
anonymous array of the codons overlapping the specified 
$position in the currently loaded sequence corresponding to the 
GFF feature object '$gf' given as its argument.

If $position is undefined, then a list of all the codons
in the given GFF subsequence is returned.

If the optional '$alleles' argument is given, then a list of
codons is returned, corresponding to the given variants at the 
$position.  The $alleles are specified as any combination of 
'A', 'T', 'C' or 'G' in a single string (case insensitive).

Note: if $alleles is given, then $position should be given!

Note that this method expects that the $gfc strand and frame (phase)
of the GFF record are properly set.  Note that one simple exon
boundary GFF, the frame may indicate that the first one or two
nucleotides are ignored in the translation, and the last codon
may also be so undefined.

=back

=cut

sub codons {
    my $self = shift ;
    my $gf   = shift ;

    confess "GFF::Sequence::codons() ERROR: sequence object undefined?\n"
	unless(defined($self->{'_SEQOBJ'})) ;

    GFF::GeneFeature->verify($gf) ;

    my $position = shift ;
    if($position and $position<$gf->start or $position>$gf->end) {
	carp "GFF::Sequence::codons() ERROR: position specified is outside gene feature coding bounds?\n" ;
	return [] ;
    }

    my $alleles  = shift ;

    if($alleles and !defined($position)) {
	confess "GFF::Sequence::codons() ERROR: undefined position for alleles?\n" ;
    }

    my $strand ;
    if($gf->strand !~ /^[\+-]$/) {
	carp "GFF::Sequence::codons() WARNING: undefined or invalid feature strand? Assuming '+'\n" ;
	$strand = '+' ;
    } else {
	$strand = $gf->strand ;
    }

    my $frame ;
    if($gf->frame eq '.') {
	carp "GFF::Sequence::codons() WARNING: undefined feature frame/phase? Assuming '0'\n" ;
	$frame = 0 ;
    } else {
	$frame = $gf->frame+0 ;
    }

    my ($start,$end) ;
    if( $strand eq '+') {
	$start = $gf->start+$frame ;
	$end   = $gf->end ;
    } else { # '-'
	$start = $gf->start ;
	$end   = $gf->end-$frame ;
    }

    my $gfc = $gf->copy() ; # for safety, use a copy
    my @codons ;
    if($position) {
	# identify the single codon hit overlapping this position
	my $cpos ;
	if( $strand eq '+') { # forward strand codon scan
	    for($cpos=$start; $cpos<=$end; $cpos+=3) {
		if($position == $cpos) {
		    $gfc->frame('0') ;
		    last ;
		} elsif($position == $cpos+1) {
		    $gfc->frame('1') ;
		    last ;
		} elsif($position == $cpos+2) {
		    $gfc->frame('2') ;
		    last ;
		}
	    }
	    $gfc->start($cpos) ;
	    $gfc->end($cpos+2) ;
	} else { # reverse strand codon scan
	    for($cpos=$end; $cpos>=$start; $cpos-=3) {
		if($position == $cpos) {
		    $gfc->frame('0') ;
		    last ;
		} elsif($position == $cpos-1) {
		    $gfc->frame('1') ;
		    last ;
		} elsif($position == $cpos-2) {
		    $gfc->frame('2') ;
		    last ;
		}
	    }
	    $gfc->start($cpos-2) ;
	    $gfc->end($cpos) ;
	}
	my $subseq = $self->subseq($gfc) ;
	if($alleles) {
	    my $spos = $gfc->frame() ;
	    my @c = split //, lc($subseq) ; # assumed to be only 3 bp long
	    foreach my $a (split //,lc($alleles)) {
		my $variant ;
		if($spos eq '0') {
		    $variant = $a ;
		} else {
		    $variant = $c[0] ;
		}
		if($spos eq '1') {
		    $variant .= $a ;
		} else {
		    $variant .= $c[1] ;
		}
		if($spos eq '2') {
		    $variant .= $a ;
		} else {
		    $variant .= $c[2] ;
		}
		push(@codons,$variant) ;
	    }
	} else {
	    push(@codons,$subseq) ;
	}

    } else {
	# return all codons in given gene feature?
	$gfc->start($start) ;
	$gfc->end($end) ;

	my $subseq = $self->subseq($gfc) ;

	# split result into codons
	for(my $i=0 ; $i+3<=length($subseq);$i+=3) {
	    push(@codons,substr($subseq,$i,3)) ;
	}
    }
    return(@codons) ;
}

=pod

=over 4

=item amino_acids($gfc,[$position,$alleles,$codontable])

A specialized version of the 'codon' method, this method returns an
anonymous array of the amino acids corresponding to the codons
which would be returned by the codons' method, with the given
method arguments (see also the 'codon' method above).

If provided, the optional $codontable is a Bio::Tools::CodonTable reference.
Otherwise, the Standard (id=1) Bio::Tools::CodonTable is used.

=back

=cut

use Bio::Tools::CodonTable ;

sub amino_acids {
    my $self = shift ;
    my $gf   = shift ;
    my $position = shift ;
    my $alleles  = shift ;

    my $codontable = shift ;

    $codontable = Bio::Tools::CodonTable->new() 
	unless(defined($codontable));

    my @codons = $self->codons($gf,$position,$alleles) ;

    grep { $_ = $codontable->translate($_); } @codons ;

    return @codons ;
}

=pod

=over 4

=item close()

This method closes a sequence handle previous initialized by 'open'.

=back

=cut

sub close {
    my $self = shift ;

    if(defined($self->{'_SEQIOH'})) {
	$self->{'_SEQIOH'}->close() ;
	undef $self->{'_SEQIOH'} ;
    } else {
	carp "GFF::Sequence::close() WARNING: no sequence file handle to close?" ;
    }
}


1;

=head2 REQUIRES

 In order to be run, this module requires:

   - Bioperl Release 0.6.* from http://www.bioperl.org

=head2 REVISIONS

 2/7/2000 (rbsk) - Creation

=cut
