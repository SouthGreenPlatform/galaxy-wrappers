
#FUNCTION - Perl Object Modules for GFF

=pod

=head1 NAME

    GFF.pm - Perl Object Modules for GFF Annotation Protocol

=head2 AUTHORSHIP

    Copyright (c) 1999, 2000

    Richard Bruskiewich <rbsk@sanger.ac.uk>
    Tim Hubbard         <th@sanger.ac.uk>

    Sanger Centre, Wellcome Trust Genome Campus, Cambs, UK
    All rights reserved.

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation

=head2 SYNOPSIS

    use GFF ;

=head1 DESCRIPTION

=head2 Overview

    GFF.pm is a Perl Object base class/module for the GFF annotation format.

    A GFF object is a base class for GeneFeatureSet, GFF::GeneFeature and
    GFF::HomolGeneFeature objects.  

=head2 How to Read Method Protocols

    Normal Perl data type notations are used for argument declarations in
    the method protocols. A backslash denotes argument passing by
    reference.  Class methods are invoked using the
    'class->method(args)' or 'method class args' Perl call formats.

=head2 SOURCE CODE

    The most current release of the Perl source code for this module is
    available at the Sanger FTP site. All bug reports may be submitted 
    to Richard Bruskiewich (rbsk@sanger.ac.uk).

=head1 DETAILS

=head2 What is GFF?

    A simple flat file exchange format for genomic DNA sequence feature 
    descriptions, into its second specification generation (Version 2).

    Can be dumped by ACEDB (and AcePerl, as GFF.pm objects)

    Full information about the protocol is
    available at http://www.sanger.ac.uk/Software/formats/GFF 

=head2 How can GFF be used?

    Convenient for simple bioinformatic computations which 
    compare overlap and clustering of features on sequences.

    One can compare similar outputs from different sources 
    (e.g. gene predictions against experimental data).

    One can combine the outputs of various genome signal 
    analysis packages, as input for higher level analysis 
    programs, and compute statistics on GFF using simple 
    grouping algorithms, also output results as GFF 
    (e.g. human chromosome 22 whole chromosome analyses).

=head2 Architecture.pm

    The whole set of GFF.pm related modules provide object oriented manipulation
    and graphical representation of GFF V.1/2 records and record sets (files).
    
    The module set consists of the following modules:
    
    The 'Core' set (imported all by 'use GFF;'):

	GFF                   - some class level methods (version, verify, trace)

	GFF::GeneFeatureSet   - collection class for (a file of) GFF records

	GFF::GeneFeature      - basic gene feature object

        GFF::Sequence         - module for sequence extraction using GFF
	                        (requires BioPerl to work)

        Aside from setting up the basic objects for manipulation, these
        modules provide rudimentary operations to manipulate GFF (see
        specific module docs for further information).

    Optional Special Purpose Modules (some selective for EXPORT symbols):
    
	GFF::Analysis         - some utilities for GFF analysis: e.g. high 
				level functions to manipulate GFF features 
				and associated sequences...

	GFF::GifGFF           - graphical representation of GFF specified 
				sequence features, including homology information

	GFF::Graph            - graphical plots of GFF scores over sequence 
				(using Raphael Leplae's Curve_plot.pm module).

	GFF::CWA              - chromosome wide analysis module for 
	                        SuperMap graphics

	GFF::CWN              - WWW navigation utilities for SuperMap's on the WWW

=head2 GFF Structure

    Flat file format, one gene feature per record, 
    with nine defined, tab delimited fields:
    
    <Seqname>:    host sequence containing the feature
    
    <Source>: 	  analytic source determining the feature
    
    <Feature>:	  primary label for the feature
    
    <Start>: 	  start coordinate (relative to Seqname)
    
    <End>: 	  end coordinate
    
    <Score>: 	  score (if any) from Source analysis of the feature
    
    <Strand>: 	  strandedness (if pertinent) of the feature
    
    <Frame>:	  reading frame (if pertinent)
    
    <Group>:      group/attribute field (optional; tag-values in GFF V.2)
    
    #comments and special ##meta-comments

=head2 (GFF V.2) Example

    ##version 2
    ##date 1999-08-04
    ##sequence-region 130N4  1  35023
    130N5  EMBL      exon       -293   240     .     +    0    Sequence "GUFFAW.1"
    130N5  GD_mRNA   exon       1100  1430     .     +    1    Sequence "130N5.1.1"
    130N5  Genscan   CDS        1100  1430   0.9     +    1    Sequence "130N5.GENSCAN.1"
    130N5  BLASTX    similarity   40    60  86.2     -    .    Target   "GUFFAW" 20 40; E_Value 0.05

=head2 EXAMPLE SCRIPT

    # Sample script using core functionality
    # See other special purpose modules for 
    # their extensions to this functionality

    #!/usr/local/perl
    use GFF ;  

    # Inputting a GFF object from a GFF input file

    my $gffi = new GFF::GeneFeatureSet;  
    open(GFFIN,"<myfilein.gff") 
		    || die "Can't open GFF input file";

    $input_filter = sub { # the filter is optional
	    my $self = shift ;   # expects a gene feature object reference...
	    ($self->source() =~/^(genscan|supported_gene)$/i) &&
	    ($self->feature() =~/^CDS$/i);  # returns non-null (include) or not
    }
    my $count = $gffi->read(\*GFFIN, $input_filter ) ; 
    print "Found $count features\n" ;
    close(GFFIN) ;
    
    # Doing some more GFF queries/computations

    my ($seqname,$start,$end) = $gffi->region ;
    print "Maximum score for $seqname is: ",$gffi->maxScore() ;

    foreach my $gf ($gffi->eachGeneFeature) {
	 if($gf->score() >= 0.7) {
	     $gf->feature('High_Scoring_Feature') ; # resets the feature field
	     print $gf->Sequence, ":", $gf->start,",",$gf->end ;
	 }
    }

    # Doing some more sophisticated GFF computations

    # filter for true genes
    $truegenes = sub{my $s = shift ;($s->source()=~/^supported_gene $/i) ;}

    # filter for predicted genes
    $predictions = sub{ my $s = shift;($self->source()=~/^genscan$/i) ;}

    my $gfft = filter($truegenes) ;
    my $gffp = filter($predictions) ;

    # test for overlap between true and predicted genes
    my $tolerance = 0 ; my $strand = 1 ;
    my $gffo = 
	    $gffp->intersect_overlap_matches($gfft,$tolerance,0,$strand) ;

    # Output the results as GFF
    open(GFFOUT,">myfileout.gff") 
		    || die "Can't open GFF output file";

    $gffo->dump_header(\*GFFOUT) ;
    $gffo->dump(\*GFFOUT) ;     # or $gffo->dump_matches(\*GFFOUT) ; 
    close(GFFOUT) ;

=cut

package  GFF;
use Carp;
use strict;
use vars qw($AUTOLOAD @ISA @EXPORT @EXPORT_OK $VERSION $TRACE);
require Exporter;
#
# @ISA has our inheritance.
#
@ISA = qw( Exporter );
#
# Place functions/variables you want to *export*, ie be visible from the caller package into @EXPORT_OK
#
@EXPORT    = qw();
@EXPORT_OK = qw();

$VERSION = '3.05' ;

my %fields = (
              version  => undef, 
	      seqname  => undef, 
);

my $GFF_VERSION = '2' ; # rbsk: 2/99 GFF Version 2 in use, new default value for global flag control

=pod

=head1 GFF.pm METHODS

=head2 GFF Class Methods and Variables

=over 4

=item version( $version )

This method may actually be called as a class (static) method or
as a virtual (instance) method. When called as a class method,
the method sets and/or queries the package default GFF protocol
version. GFF objects internally record the current default
version value when they are created. If the method is called by
such a derived object, it sets and/or queries the GFF object's
version (thus, one can in principle create and manipulated GFF
objects of various versions). If not otherwise called, the
module currently assumes that Version 2 GFF objects will be
manipulated. If the format of the anticipated input is not the
default protocol version, then this method must be called first,
prior to calling other format sensitive methods, to set the
correct version. The method may be called with or without an
argument. In both cases, the method returns the current (or
newly set) version value, either the package ("class") default
or the GFF object value.

=back

=cut

sub version  {
    my $self     = shift ;
    my $version  = shift ;

    if( ref($self) and GFF->verify($self) ) { 
	if(defined($version) and $version) {
	    # invoked by a GFF object, as an instance method
	    if( $version =~/^(\d+)$/  ) { 
		$self->{'version'} = $1 + 0 ; # reset the GFF object version
	    } else {
		TRACE(2,ref($self)."->version() $version invalid: just returning existing object version...\n") ;
	    }
	}
        return $self->{'version'} ;
    } else { 
        # invoked as a class method!
	if(defined($version) and $version) {
	    if( $version =~/^(\d+)$/  ) { # can reset the default version
		$GFF_VERSION = $1 + 0 ; 
		TRACE(2,"GFF package default version reset to $GFF_VERSION\n") ;
	    } else {
		TRACE(2,ref($self)."->version() $version invalid: just returning module default version...\n") ;
	    }
	}
        return $GFF_VERSION ;
    }
}

=pod

=over 4

=item traceFile( \*FILEHANDLE )

Redirects trace output to specified file (default: \*STDERR)

=item trace( $level )

Setting $level to a number greater than 1 triggers a verbose GFF.pm et al.
module diagnostic mode output to STDERR. If $level argument is
omitted, then the trace is turned off. Examples: verify()
(above) reports success; GFF::read() provides user feedback
during file input. Generally speaking, the verbosity of the output
depends upon the $level given: 1 = minimum trace (most useful for
monitoring certain iterative operations like reading in a file),
2 = verbose (debug) trace (not so useful except for GFF.pm coders)

=back

=cut

my $TRACE = 0 ; # user can reset this...

my $TRACEFILE = \*STDERR ;

sub traceFile { $TRACEFILE = shift ; }

sub trace { my $state = shift ; $TRACE = (defined $state and $state)? $state : 0 ; }

sub TRACE {
    my $level = shift ; 
    if($TRACE >= $level) {
	if( $level == 1 ) { # level 1 is a simple trace mode
	    print $TRACEFILE @_  ;
	} else {
	    carp @_ ;
	}
    }
}

=pod

=over 4

=item verify( $ref, $nodie )

Any derived class can invoke this method, to verify that $ref is
a properly defined object reference in that class. If invoked by
the GFF class, then verifies that the $ref belongs to any one of
the GFF derived classes. Method 'dies' upon failure, unless the 
$nodie variable is non-null (default: null), in which case, the
function merely carps a warning and returns null for failure.

=back

=cut

sub verify {
    my $self  = shift ;
    my $class = ref($self) || $self ;
    $class = '(GFF::(GeneFeature|HomolGeneFeature))' if $class eq 'GFF::GeneFeature' ;
    $class = '(GFF(::(GeneFeatureSet|GeneFeature|HomolGeneFeature)))' if $class eq 'GFF' ;
    my $ref   = shift ;

    my $nodie = shift ;
    $nodie = 0 if !defined($nodie) ;

    if(!defined $ref) {
        if($nodie) {
	    carp "GFF::verify(): undefined object!" ;
	    return 0 ;
	} else {
	    confess "GFF::verify(): undefined object!"  ;
        }
    } elsif(!ref($ref)) {
        if($nodie) {
	    carp "GFF::verify(): not a reference!" ;
	    return 0 ;
	} else {
            confess "GFF::verify(): not a reference!" ;
        }
    } elsif( ref($ref) !~ /$class/ ) {
        if($nodie) {
	    carp "GFF::verify(): ".ref($ref)." is not a $class class object!"  ;
	    return 0 ;
	} else {
            confess "GFF::verify(): ".ref($ref)." is not a $class class object!"  ;
        }
    }

    TRACE(2,"$class\:\:verify($ref) OK!\n") ;
    return 1 ;
}

=head2 GFF Construction Methods

=over 4

=item new

Not generally used; derived classes generally override this class.

$version is GFF version of the object. $seqname is the associated 
sequence name of the object.

=back

=cut

sub new {
    my $ref     = shift ;
    my $class   = ref($ref) || $ref ;
    my $version = shift ;
    my $seqname = shift ;

    my $self = { '_permitted' => \%fields } ;
    
    if( defined $seqname ) {
        $self->{'seqname'} = $seqname ;
    }
    bless $self, $class ;

    if(defined($version) && $version =~ /\d+/) {
	$self->{'version'} = $version ;
    } else {
	$self->{'version'} = $GFF_VERSION ;
    }

    return $self ;
}

=head2 GFF Object Input/Output Methods

=over 4

=item AUTOLOAD

Not generally used; derived classes generally override this class.

Autoloads up my field name set/get.


=back

=cut

sub AUTOLOAD {
    my $self = shift;
    my $class = ref($self) || carp "\'$self\' is not an object of mine!";
    my $name = $AUTOLOAD;

    # don't propagate DESTROY messages...

    $name =~ /::DESTROY/ && return;

    $name =~ s/.*://; #get only the bit we want
    unless (exists $self->{'_permitted'}->{$name} ) {
	confess "In type $class, can't access $name - probably passed a wrong variable into $class  ";
    }
    if (@_) {
	return $self->{$name} = shift;
    } else {
	return $self->{$name}
    }
}


use lib 'GFF' ;			
use GFF::GeneFeature ;
use GFF::HomolGeneFeature ;
use GFF::GeneFeatureSet ;

1; # required for GFF

__END__

=pod

=head1 REQUIRES

    Bioperl 0.6.1+ is required for GFF::Sequence to work.

    Raphael Leplae's Curve_plot.pm module is 
    required for the GFF::Graph module to work

=head1 REVISIONS

 3.05    22/06/2000 rbsk: see GeneFeature.pm and HomolGeneFeature.pm for bug fixes

 3.03     7/4/2000  rbsk: various updates to subsidiary modules, in particular,
                         new memory efficient data implementation for GFF::GeneFeature.pm
 
 3.02    14/02/2000 rbsk: version() method bug: should return $self->version
                         when invoked as an object method
 3.01    02/02/2000 rbsk: novel, more memory efficient data structure 
                         implementation for GeneFeature and GeneFeatureSet's

 2.12    10/11/99 - rbsk: see GeneFeatureSet.pm, GifGFF.pm and Graph.pm
 2.11     4/10/99 - rbsk: see GeneFeatureSet.pm and GeneFeature.pm
 2.10    27/08/99 - rbsk: GeneFeature parse_group bug fix
 2.09    26/08/99 - rbsk: GeneFeature.pm changes...
 2.08     28/1/99 - rbsk: module renamed from GFFObject.pm => GFF.pm
                     overall GFF object module hierarchy remolded
                      GFF.pm              => GFF::GeneFeatureSet.pm
                      GeneFeature.pm      => GFF::GeneFeature.pm
                      HomolGeneFeature.pm => GFF::HomolGeneFeature.pm
                   all in GFF subdirectory
    22/4/99 - rbsk: added trace(), verify() and trace() methods
    16/3/99 - rbsk: created this class to abstract "version" method
                   from GFF, GeneFeature and HomolGeneFeature.


=head1 SEE ALSO

    http://www.sanger.ac.uk/Software/GFF/GeneFeature.shtml,
    http://www.sanger.ac.uk/Software/GFF/GeneFeatureSet.shtml
    http://www.sanger.ac.uk/Software/GFF/HomolGeneFeature.shtml,
    http://www.sanger.ac.uk/Software/GFF/GifGFF.shtml
    http://www.sanger.ac.uk/Software/GFF/Analysis.shtml

=cut
