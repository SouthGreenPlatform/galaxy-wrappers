=pod

=head1 NAME

GFF::Analysis.pm - Perl utility library for 'Gene Finding Format' (GFF) analysis routines
                   using the GFF.pm Perl Object libraries

=head1 SYNOPSIS

# include what functions you need... GFF::Analysis.pm contains an implicit 'use GFF ;'
use GFF::Analysis qw(constructGene makeGenes mRNA featureLengthStats
                     segregateGeneFeatures normalize mergeGeneFeatures 
                     cleanUpSeqName normalize_mRNA gfProximity);  

=head1 DESCRIPTION

GFF::Analysis (derived from GFF) is a utility
library for for Gene Finding Format, built upon
the GFF.pm module library.

Exports: 
    constructGene() 
    makeGenes()
    mRNA() 
    segregateGeneFeatures()
    cleanUpSeqName()  #  name cleanup protocol used by normalize
    normalize()
    mergeGeneFeatures()
    normalize_mRNA() # calls segregateGeneFeatures(), normalize(), and mergeGeneFeatures() sequentially
    gfProximity()

=head1 AUTHORSHIP

Copyright (c) 1999 
Created by Richard Bruskiewich <rbsk@sanger.ac.uk>

Sanger Centre, Wellcome Trust Genome Campus, Cambs, UK
All rights reserved.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation

=head1 SOURCE CODE

The most current release of the Perl source code 
for this module is available through a link at 
http://www.sanger.ac.uk/Software/GFF/
All bug reports may be submitted to Richard 
Bruskiewich (rbsk@sanger.ac.uk). Future releases may
be issued as Bioperl archived modules.

=cut

package GFF::Analysis ;

use Carp;
use strict;
use vars qw($AUTOLOAD @ISA @EXPORT @EXPORT_OK $VERSION);
require Exporter;
use GFF ;

$VERSION = '2.16' ; # Based upon Version 2.0 GFF
#
# @ISA has our inheritance.
#
@ISA = qw( Exporter GFF );
#
# Place functions/variables you want to *export*, ie be visible from the caller package into @EXPORT_OK
#
@EXPORT    = qw();
@EXPORT_OK = qw( constructGene 
                 makeGenes 
                 mRNA 
                 featureLengthStats
                 segregateGeneFeatures
                 normalize
                 mergeGeneFeatures 
                 cleanUpSeqName 
                 normalize_mRNA
                 gfProximity
               );
#
# @ISA has our inheritance.
#
@ISA = ( 'Exporter' );

=pod

=head1 GFF::Analysis Construction Methods

=over 4

=item new()

Class method to construct a new empty
GFF::Analysis object.

=back

=cut

my %fields = (
	      );
#
# creates a new object
#
sub new {
    my $ref     = shift;
    my $class   = ref($ref) || $ref ;

    my $self = {} ;

    $self->{'_permitted'} = \%fields ;

    bless $self, $class;
    return $self;
}

=pod

=head1 Methods to do calculations on GFF::GeneFeatureSet objects

Note: These methods are not 'object' invoked, but take a
GFF::GeneFeatureSet reference as their first argument.

=over 4

=item constructGene($gffi)

Given a GeneFeatureSet object ('$gfs') GeneFeatures with 
<feature> fields specifically labelled with 'exon' 
and possibly 'promoter', 'transcription_start' and/or
'polyA_signal' tags, and belonging to a single gene 
as defined by a common [group] field label, this method
returns an augmented GeneFeatureSet object fully describing
a 'gene' containing introns, UTRs and flanking sequences
inferred from the original GeneFeature set. 

GFF::GeneFeatures in the calling GeneFeatureSet object 
are assumed to all be on the same <strand>, from the same 
<seqname>, <source> and group_value('Sequence') named 
gene (i.e. makeGenes() clustering) thus, the method 
looks at the first GeneFeature encountered in the
 GeneFeatureSet for all these values! The <frame> and 
<score> are assumed to be irrelevant for all features 
added in this method (e.g. introns), and is thus set 
to '.'. The incoming GeneFeatures are also assumed to 
be non-overlapping, since this assumption drives the
identification of 'inter' GeneFeature gaps ('introns' 
et al.) Also, if the first (and/or last) GeneFeature 
start (end) does (do) not coincide with the start (end)
of the <seqname> region range, then the 5' (and 3') 
flanking regions are inferred and so labelled in the 
field. This latter labelling is also influenced by 
the presence of 'promoter', 'transcription_start' 
and/or 'polyA_signal' GeneFeatures.

=back

=cut

sub constructGene{
    my $self  = shift ;  # isa GFF::GeneFeatureSet reference
    my $noseq = shift ;  # suppress creation of a sequence record

    return undef unless(GFF::GeneFeatureSet->verify($self,1)) ; # nodie check
    return $self if $self->count() <= 1 ;

    my ($gene,$start,$end) = $self->region() ;
    my $inStart = defined $start ? $start : $self->start_range() ;
    my $inEnd   = defined $end   ? $end   : $self->end_range() ;
#
#   Determine default values for <source>, et al.
#   assumed to be uniform in the GFF::GeneFeatures of the invoking GeneFeatureSet object
#   hence, simply extracted from the first GFF::GeneFeature of the GeneFeatureSet object.
#
    my $gfa = $self->{'feature'} ; # is an array ref?
    my $gf = @{$gfa}[0] ;          # look at the first gf in the list
    my $seqname  = $gf->seqname() ;
    my $source   = $gf->source() ;
    my $strand   = $gf->strand() ;
#
#   Now attempt to identify the introns and extragenic regions
#
    my $newGF = sub {
        my $gfStart = shift ;
	my $gfEnd   = shift ;
	my $feature = shift ;
	my $gfi = GFF::GeneFeature->new($self->version()) ;
	$gfi->seqname($seqname) ;
	$gfi->source($source) ;
	$gfi->feature($feature) ;
	$gfi->start($gfStart) ;
	$gfi->end($gfEnd) ;
	$gfi->strand($strand) ;
	$gfi->group_value('Sequence',0,$gene) ;
        $self->addGeneFeature($gfi) ;
    };
    #
    # GFF::GeneFeature order is strand sensitive... 
    # Reverse sort order for negative stranded genes?
    #
    ##################### positive stranded gene #########################
    if( $strand eq '+' ) {
        my @gflist = $self->order_gf() ;
        $gf = $gflist[0] ; # look at the first gf...
#
#       Check for 5' flanking sequence found in front of first GFF::GeneFeature
# 
	if($gf->feature() =~ /(promoter|transcription_start)/) {
            if($inStart < $gf->start()) { # leader sequence?
	        $newGF->($inStart,$gf->start()-1,'5\'flanking_region') ;
	    }
	    $inStart = $gf->end() + 1 ; # leap over the gf end
	    shift @gflist ;
            my $gf2 = $gflist[0] ; # look at the second gf...
	    # and add a 5'UTR
	    $newGF->($inStart,$gf2->start()-1,'5\'UTR') ;
	    $inStart = $gf2->end() + 1 ; # leap over the current gf end & gf itselfq
	} else { # GFF::GeneFeature is an exon?
            if($inStart < $gf->start()) { # leader sequence == UTR?
                $newGF->($inStart,$gf->start()-1,'5\'UTR') ;
	    }
	    $inStart = $gf->end() + 1 ; # leap over the current gf end & gf itselfq
        }
	shift @gflist ;
#
#       ...then go through the rest of the GFF::GeneFeature list...
#       
        foreach $gf (@gflist) {
	    if($inStart < $gf->start()) {
	        if($gf->feature() eq 'polyA_signal') {
	            $newGF->($inStart,$gf->start()-1,'3\'UTR') ;
	        } else {
	            $newGF->($inStart,$gf->start()-1,'intron') ;
		}
	    }
	    $inStart = $gf->end() + 1 ; # leap over the gf end
        }
#
#       ...check for trailing sequence after final GFF::GeneFeature
#
        if(defined $end and $inStart < $end) {
	    if($gf->feature() eq 'polyA_signal') {
	        $newGF->($inStart,$end,'3\'flanking_region') ;
	    } else {
	        # last GFF::GeneFeature was an exon? assume trailer is a 3'UTR
	        $newGF->($inStart,$end,'3\'UTR') ;
            }
        }
    #
    ##################### negative stranded gene #########################
    #
    } else { 
        my @gflist = $self->order_gf(1) ;
        $gf = $gflist[0] ; # look at the first gf...
#
#       Check for 5' flanking sequence found in front of first GFF::GeneFeature
# 
	if($gf->feature() =~ /(promoter|transcription_start)/) {
            if($inEnd > $gf->end()) { # leader sequence?
	        $newGF->($gf->end()+1,$inEnd,'5\'flanking_region') ;
	    }
	    $inEnd = $gf->start() - 1 ; # leap over the gf end
	    shift @gflist ;
            $gf = $gflist[0] ; # look at the second gf...
	    # and add a 5'UTR
	    $newGF->($gf->end()+1,$inEnd,'5\'UTR') ;
	} else { # GFF::GeneFeature is an exon?
            if($inEnd > $gf->end()) { # leader sequence == UTR?
                $newGF->($gf->end()+1,$inEnd,'5\'UTR') ;
	    }
        }
	$inEnd = $gf->start() - 1 ; # leap over the current gf start & gf itself
        shift @gflist ;
#
#       ...then go through the rest of the list...
# 
        foreach $gf (@gflist) {
	    if($inEnd > $gf->end()) {
	        if($gf->feature() eq 'polyA_signal') {
	            $newGF->($gf->end()+1,$inEnd,'3\'UTR') ;
	        } else {
	            $newGF->($gf->end()+1,$inEnd,'intron') ;
		}
	    }
	    $inEnd = $gf->start() - 1 ; # leap over the gf end
        }
#
#       ...check for trailing sequence after final GFF::GeneFeature
#
        if(defined $start and $inEnd > $start) {
	    if($gf->feature() eq 'polyA_signal') {
	        $newGF->($start,$inEnd,'3\'flanking_region') ;
	    } else {
	        # last GFF::GeneFeature was an exon? assume trailer is a 3'UTR
	        $newGF->($start,$inEnd,'3\'UTR') ;
            }
        }
    }
#
#   Finally, construct a 'sequence' GFF::GeneFeature first
#   based upon the region range of the gene
#
    unless(defined($noseq) and $noseq) {
	$gf = GFF::GeneFeature->new($self->version()) ; # score defaults to '.' 
	$gf->seqname($seqname) ; 
	$gf->source($source) ;
	$gf->feature('sequence') ;
	$gf->start($start) ;
	$gf->end($end) ;
	$gf->strand($strand) ;
	$gf->group_value('Sequence',0,$gene) ;
	$self->addGeneFeature($gf) ;
    }
#
#   the new gene GeneFeatureSet...
#
    return $self ;
}

=pod

=over 4

=item makeGenes($gffi)

After clustering a GeneFeatureSet set of
predicted exons, promoter, polyA's etc. by 'gene'
groups (i.e. by Version 1 [group] tags or by
Version 2 [group] field 'Sequence' tag-values),
this method uses the
GeneFeatureSet::constructGene method to infer
additional GeneFeatures (e.g. introns, [5'|3']
UTRs and Flanking regions). The method then
returns all the GeneFeatures (old and new) in a
new GeneFeatureSet object.

=back

=cut

sub makeGenes {
   my $self = shift ;
   GFF::GeneFeatureSet->verify($self) ;

   my $newSet = $self->new($self->version(),$self->region()) ;
#
#  Cluster all GFF::GeneFeatures by 'Sequence' identities
#  into a hash of 'gene' GeneFeatureSet objects...
#
   my $gene_filter=sub{
       my $self = shift ;
       return $self->group_value('Sequence') ;
   };
   
   my %genes = $self->group($gene_filter) ;
#
#  Now, add inferred features to predicted genes...
#
   foreach my $name (keys %genes) {
       my $gffg = $genes{$name} ;
       constructGene($gffg) ;
       $newSet->addGFF($gffg) ;
   }
   return $newSet ;
}

=pod

=over 4

=item mRNA($gffi, $seq, $pattern, $offset)

Method to return a single string of a subsequence 
representing a mRNA or similar gapped entity 
represented by the gene features in the invoking 
object, whose <feature> field matches the $pattern.
The method expects a string '$seq' corresponding 
to the sequence from which the features are to be 
extracted. Returns a concatenated string of all 
the subsequences defined by the matching gene features.

Note that '$seq' is a Perl '0' based string, however
GFF (and biology) convention sets the first base
sequence coordinate to '1'...

If the GFF defining the gene structure is offset
in coordinates relative to the $seq, then the
$offset should also be provided, that is,
the $offset is the value to be added to base
position '1' of the sequence, to give the GFF
coordinate lining up with that position.

=back

=cut

sub mRNA {
    my $gffi    = shift ; # reference to a GFF description of gene
    GFF::GeneFeatureSet->verify($gffi) ;

    my $seq     = shift ; # reference to a string of sequence
    my $pattern = shift ; # simple identifier (e.g. 'exon') 
                          # or regex describing features of interest

    my $offset  = shift ;  # GFF bulk offset
    $offset=0 if !defined($offset) ;

    my $ts_filter = sub{
        my $self = shift ;
        my $match ;
        eval {  ($match) =( $self->feature() =~ /^$pattern$/ ) } ;
        if($@) { 
	    die "*** Error: Improper mRNA feature pattern: $pattern !?!" ;
        }
        if($match) {
	    return 1 ;
        } else {
            return 0 ;
	}
    } ;

    my $gffts = $gffi->filter($ts_filter) ;
    my $subseq = "" ; # isa string...
    foreach my $gf ($gffts->eachGeneFeature()) {
	# remember that Perl strings begin a '0' whereas first sequence string base is 1
	# also may need to offset the GFF to the sequence origin
	$subseq .= substr( ${$seq}, $gf->start()-$offset-1, ($gf->end() - $gf->start() + 1) ) ;
    }
    return $subseq ; # isa string..., forward stranded only
}


=pod

=over 4

=item featureLengthStats($gffo,\%statab,$label)

This method applies the GFF::GeneFeatureSet::lengthStats()
method to a given GeneFeatureSet, returning the results in
a primary hash (passed by reference,\%statab ) indexed under the 
given $label, and returned as a secondary level hash reference.

    i.e.  \%statab->{"$label"}->{'<data_label>'}

The '<data_label> secondary hash keys for these statistics 
are as follows:

    'M'    == mean length
    'SD'   == standard deviation of lengths
    'N'    == total number of features
    'Min'  == minimum length
    'Max'  == maximum length
    'Cov'  == total sum of lengths ('coverage')
    'Cov2' == total sum of lengths squared
    'LenC' == reference to an array of feature incidence 
	      counts for each class of length
    
A side effect of the setting of the table is that these values
are returned (in an array context) as a list, in the order
indicated above.

=back

=cut

sub featureLengthStats {
    my $gffo   = shift ; # a reference to a GFF::GeneFeatureSet object
    GFF::GeneFeatureSet->verify($gffo) ;

    my $statab = shift ; # isa reference to a hash
    my $label  = shift ;

    $statab->{"$label"} = {} ;

    ( $statab->{"$label"}->{'M'},
      $statab->{"$label"}->{'SD'},    
      $statab->{"$label"}->{'N'}, 
      $statab->{"$label"}->{'Min'},    
      $statab->{"$label"}->{'Max'},    
      $statab->{"$label"}->{'Cov'},
      $statab->{"$label"}->{'Cov2'}, 
      $statab->{"$label"}->{'LenC'}
    ) = $gffo->lengthStats() ;    # length statistics for '$label' features
}

=pod

=over 4

=item segregateGeneFeatures($gffi,\%statab,$trace)

Given a GeneFeatureSet containing 'sequence', 'exon' and 'CDS' records,
this method returns a list of four $gff object pointers representing each
segregated subsets for each of 'genes',  'pseudogenes', 'exons'
and 'CDSs', respectively. 

If a suitable reference '\%statab' to a hash table is given, then the 
method also compiles statistics for each subset into that table using 
the featureLengthStats() method (see above). That is, the table is of 
the form:

        \%statab->{('Gene'|'Exon'|'CDS')}->{'<data item>'}

where '<data item>'s are statistics as returned by the 
featureLengthStats().

Another side effect of the method is that mRNA/coding_exon and 
CDS/exon redundancies are filtered out of the dataset.

The optional '$trace' boolean flag, when non-null, turns on GFF module tracing.

=back

=cut

sub segregateGeneFeatures {

    my $gffi   = shift ;
    GFF::GeneFeatureSet->verify($gffi) ;

    my $statab = shift ;
    my $trace  = shift ;

    
    croak "GFF::Analysis::segregateGeneFeatures(): invalid statistics hash reference argument: $statab\n" 
	if $statab && ref($statab) !~ /HASH/ ;

    my $echo ;
    if(defined($trace) and $trace) {
        GFF::trace(1) ;
	$echo = \*STDERR ;
    } else {
	$echo = 0 ;
    }

    ################################################################################################
    # I don't want to use (GD_mRNA or supported_mRNA) and coding_exon records
    # since they are already properly covered by *CDS coding_exons
    # I also don't want to use (GD_CDS or supported_CDS) and exon records
    # since they should be already properly covered by *mRNA exons
    ################################################################################################
    my $Redundancyfilter = sub {
        my $self = shift ;
        if( !($self->source()  =~ /(Pseudogene|GD_mRNA|supported_mRNA|GD_in_progress|GD_not_confirmed)/i and 
              $self->feature() =~ /^(coding_exon|CDS)$/   or
              $self->source()  =~ /(GD_CDS|supported_CDS)/i and 
              $self->feature() =~ /^exon$/
	     )
          ) 
	{
	    return 1 ;
        } else {
	    GFF::TRACE(1,"*** Redundant: sequence(",$self->group_value('Sequence'),
		       "), source(",$self->source(),"), feature(",$self->feature(),")\n") ;
            return 0 ;
        }
    } ;

    print STDERR "\nEntering segregateGeneFeatures() with a total of ",$gffi->count(),
                 " gene features\n\tFiltering out redundancy...\n" if defined($trace) and $trace ;
    my $gffnr = $gffi->filter($Redundancyfilter) ;

    ########################################################
    # Separate out transcript sequence, exon and CDS records
    ########################################################
    my $truegene_filter = sub {
        my $self = shift ;
	my $type = $self->group_value('Type') ;
        if( $self->feature() =~ /^sequence$/i and
	    $self->source()  !~ /pseudogene/i and
	    (!defined($type) or $type !~  /^pseudogene$/i)) {
	    return 1 ;
        } else {
            return 0 ;
        }
    } ;

    my $pseudogene_filter = sub {
        my $self = shift ;
	my $type = $self->group_value('Type') ;
        if( $self->feature() =~ /^sequence$/i and
	   (($self->source()  =~ /pseudogene/i) or
	    (defined($type) and $type =~  /^pseudogene$/i))) {
	    return 1 ;
        } else {
            return 0 ;
        }
    } ;

    my $exon_filter = sub {
        my $self = shift ;
        if( $self->feature() =~ /^exon$/i) {
	    return 1 ;
        } else {
            return 0 ;
        }
    } ;

    my $cds_filter = sub {
        my $self = shift ;
        if( $self->feature() =~ /^(CDS|coding_exon)$/i) {
	    return 1 ;
        } else {
            return 0 ;
        }
    } ;

    #################################################################
    # Segregate transcript sequences, exons && CDS's ######################
    #################################################################
    print STDERR "\tSegregating transcript components...\n" if defined($trace) and $trace ;

    my $gffg = $gffnr->filter($truegene_filter) ;
    print STDERR "\t\t",$gffg->count()," gene sequences\n" if defined($trace) and $trace ;

    my $gffp = $gffnr->filter($pseudogene_filter) ;
    print STDERR "\t\t",$gffp->count()," pseudogene sequences\n" if defined($trace) and $trace ;

    my $gffe = $gffnr->filter($exon_filter) ; # gene+pseudogene
    print STDERR "\t\t",$gffe->count()," exons\n" if defined($trace) and $trace ;

    my $gffc = $gffnr->filter($cds_filter) ;  # gene+pseudogene?
    print STDERR "\t\t",$gffc->count()," CDSs\n" if defined($trace) and $trace ;


    ###############################
    # Eliminate exact duplicates #######
    ###############################
    print STDERR "\tEliminate exact sequence feature duplications...\n" if defined($trace) and $trace ;
    my $gffndg = $gffg->self_overlap($echo,1,1,'Sequence') ; 
    print STDERR "\t\t",$gffndg->count()," gene sequences left\n" if defined($trace) and $trace ;

    my $gffndp = $gffp->self_overlap($echo,1,1,'Sequence') ; 
    print STDERR "\t\t",$gffndp->count()," pseudogene sequences left\n" if defined($trace) and $trace ;

    my $gffnde = $gffe->self_overlap($echo,1,1,'Sequence') ;
    print STDERR "\t\t",$gffnde->count()," exons left\n" if defined($trace) and $trace ;

    my $gffndc = $gffc->self_overlap($echo,1,1,'Sequence') ;
    print STDERR "\t\t",$gffndc->count()," CDSs left\n" if defined($trace) and $trace ;


    ##############################################################
    # Compute some (unnormalized feature) length statistics ##############
    ##############################################################
    if($statab) {
	print STDERR "\tComputing length statistics...\n" if defined($trace) and $trace ;
	featureLengthStats($gffndg,$statab,'Gene') ;
	featureLengthStats($gffndp,$statab,'Pseudogene') ;
	featureLengthStats($gffnde,$statab,'Exon') ;
        featureLengthStats($gffndc,$statab,'CDS') ;
    }

    print STDERR "Exiting segregateGeneFeatures()\n\n" if defined($trace) and $trace ;

    GFF::trace(0) if defined($trace) and $trace ;

    return ($gffndg, $gffndp, $gffnde, $gffndc) ;
}

=pod

=over 4

=item mergeGeneFeatures($gffg,$gffp,$gffe,$gffc,\%statab,$trace,$raw)

This method performs the reverse operation to that of segregateGeneFeatures(),
taking four types of features - 'gene (sequence)', 'pseudogene (sequence)', 'exon' and 'CDS' records -
and remerging them into cohesive 'gene' sets based upon overlapping GFF coordinates.
The method then returns  a list of references to the each of the 'true' and 'pseudo'
GeneFeatureSets.

Along the way, a further set of statistics may be (optionally) computed and stored into a 
dereferenced hash table \%statab passed to the function. These statistics pertain
rather to 'exons' and 'CDSs' per (pseudo)gene, and are stored in the hash at the primary 
level under the 'Gene' key, and at the secondary level in dereferenced hashes
with 'Exon' and 'CDS' keys, then the specific data items, e.g.

    \%statab->{Gene}->{(Exon|CDS|Transcript|Translation)}->{<data item>}

The '<data item>'s are as follows:

    'M'     == mean features per gene
    'SD'    == standard deviation of features per gene
    'X'     == total sum of given features
    'X2'    == total sum of given features squared
    'GeneC' == reference to an array of gene incidence 
               counts for each class of 'per gene' values

The optional '$trace' boolean flag, when defined and non-null, turns on GFF module tracing.

The optional '$raw' boolean flag, when defined and non-null, suppresses the 
addition of 'intron' and UTR features to the reconstructed gene descriptions.

=back

=cut

sub mergeGeneFeatures {

    my $gffg   = shift ;
    GFF::GeneFeatureSet->verify($gffg) ;

    my $gffp   = shift ;
    GFF::GeneFeatureSet->verify($gffp) ;

    my $gffe   = shift ;
    GFF::GeneFeatureSet->verify($gffe) ;

    my $gffc   = shift ;
    GFF::GeneFeatureSet->verify($gffc) ;


    my $statab = shift ;
    my $trace  = shift ;
    my $raw    = shift ;

    croak "GFF::Analysis::mergeGeneFeatures(): invalid statistics hash reference argument: $statab\n" 
	if $statab && ref($statab) !~ /HASH/ ;

    GFF::trace(1) if defined($trace) and $trace ;
    $raw=0 if !defined($raw) ;

    my $pseudogene_features = sub {
	my $gf = shift ;
	if($gf->source() =~ /pseudogene/i) {
	    1 ;
	} else {
	    0 ;
	}
    } ;

    ########################################################
    # Gather corresponding exons/CDS's into sequence sets
    ########################################################
    print STDERR "\nEntering GFF::Analysis::mergeGeneFeatures()...\n" if defined($trace) and $trace ;

    ############################
    # First, for 'true' genes...
    ############################
    my $gffgenes = new GFF::GeneFeatureSet(2,$gffg->region()) ; # limitation: can't copy over sequences?
    $gffgenes->addGFF($gffg) ;

    my ($Ng, $Xeg,$Xcg,$X2eg,$X2cg,$Xts,$Xcs,$X2ts,$X2cs,%egGeneC,%cgGeneC,%ts_cov,%cs_cov) ;
    $Ng = $Xeg = $Xcg = $Xts = $Xcs = $X2eg = $X2cg = $X2ts = $X2cs = 0 ;

    my $gffte = $gffe->exclude($pseudogene_features) ;
    my $gfftc = $gffc->exclude($pseudogene_features) ;

    foreach my $gf ($gffg->eachGeneFeature()) {
	$Ng++ ;
	my $seqname = $gf->group_value('Sequence') ;
	if(!defined($seqname)) {
	    carp "### GFF::Analysis::mergeGeneFeatures(): unknown seqname for gene:".$gf->dump_string()."\n" ;
	    next ;
	}

	############
	# Exons ############
        ############
	my $gfset = $gffte->intersect_range($gf->start,$gf->end,1,$gf->strand) ; # should be exact and strand matching...
	$gfset->region($seqname) ;
	if(my $count = $gfset->count()) {
	    $gfset = $gfset->rewriteField('GROUP','Sequence', $seqname, undef,'Old_Sequence'  ) ;

	    if($statab) {
		# exon in gene data
		print STDERR "### GFF::Analysis::mergeGeneFeatures() warning: gene sequence '$seqname' has no exons?\n" if !$count ;
		if($count) {
		    $Xeg  += $count ;
		    $X2eg += $count*$count ;
		}
		if(!exists($egGeneC{"$count"})) {
		    $egGeneC{"$count"} = 1  ;
		} else {
		    $egGeneC{"$count"}++  ;
		}
		my $transcript_size = 0 ;
		if($count) {
		    $transcript_size = $gfset->mask_length() ;
		    $Xts  += $transcript_size ;
		    $X2ts += $transcript_size*$transcript_size ;
		}
		if(!exists($ts_cov{"$transcript_size"})) {
		    $ts_cov{"$transcript_size"} = 1  ;
		} else {
		    $ts_cov{"$transcript_size"}++  ;
		}
	    }
	    $gfset = constructGene($gfset,1) unless($raw) ; # add introns et al. w/o a sequence record
	    $gffgenes->addGFF($gfset) ;
	}

        ############
	# CDS's ############
        ############
	$gfset = $gfftc->intersect_range($gf->start,$gf->end,1,$gf->strand) ; # should be exact and strand sensitive...
	$gfset->region($seqname) ;
	if(my $count = $gfset->count()) {
	    $gfset = $gfset->rewriteField('GROUP','Sequence', $seqname, undef,'Old_Sequence'  ) ;
	    if($statab) {
		# CDS in gene data
		print STDERR "### GFF::Analysis::mergeGeneFeatures() warning: gene sequence '$seqname' has no CDSs?\n" if !$count ;
		if($count) {
		    $Xcg  += $count ;
		    $X2cg += $count*$count ;
		}
		if(!exists($cgGeneC{"$count"})) {
		    $cgGeneC{"$count"} = 1  ;
		} else {
		    $cgGeneC{"$count"}++  ;
		}
		my $translation_size = 0 ;
		if($count) {
		    $translation_size = $gfset->mask_length() ;
		    $Xcs  += $translation_size ;
		    $X2cs += $translation_size*$translation_size ;
		}
		if(!exists($cs_cov{"$translation_size"})) {
		    $cs_cov{"$translation_size"} = 1  ;
		} else {
		    $cs_cov{"$translation_size"}++  ;
		}
	    }
	    $gfset = constructGene($gfset,1) unless($raw) ; # add introns et al. w/o a sequence record
	    $gffgenes->addGFF($gfset) ;
	}
    } # foreach 'true' gene...

    ##########################
    # then for pseudogenes...
    ##########################
    my $gffpseudogenes = new GFF::GeneFeatureSet(2,$gffp->region()) ; # limitation: can't copy over sequences?
    $gffpseudogenes->addGFF($gffp) ;

    my ($Np, $Xep,$Xcp,$X2ep,$X2cp,%epGeneC,%cpGeneC) ;
    $Np = $Xep = $Xcp = $X2ep = $X2cp = 0 ;

    my $gffpe = $gffe->filter($pseudogene_features) ;
    my $gffpc = $gffc->filter($pseudogene_features) ;

    foreach my $gf ($gffp->eachGeneFeature()) {
	$Np++ ;
	my $seqname = $gf->group_value('Sequence') ;
	if(!defined($seqname)) {
	    carp "### GFF::Analysis::mergeGeneFeatures(): unknown seqname for pseudogene:".$gf->dump_string()."\n" ;
	    next ;
	}

	############
	# Exons ############
        ############
	my $pfset = $gffpe->intersect_range($gf->start,$gf->end,1,$gf->strand) ; # should be exact and strand matching...
	$pfset->region($seqname) ;
	if(my $count = $pfset->count()) {
	    $pfset = $pfset->rewriteField('GROUP','Sequence', $seqname, undef,'Old_Sequence' ) ;

	    if($statab) {
		# exon in gene data
		print STDERR "### GFF::Analysis::mergeGeneFeatures() warning: pseudogene sequence '$seqname' has no exons?\n" if !$count ;
		if($count) {
		    $Xep  += $count ;
		    $X2ep += $count*$count ;
		}
		if(!exists($epGeneC{"$count"})) {
		    $epGeneC{"$count"} = 1  ;
		} else {
		    $epGeneC{"$count"}++  ;
		}
	    }

	    $pfset = constructGene($pfset,1) unless($raw) ; # add introns et al. w/o a sequence record
	    $gffpseudogenes->addGFF($pfset) ;
	}
        ############
	# CDS's ############
        ############
	$pfset = $gffpc->intersect_range($gf->start,$gf->end,1,$gf->strand) ; # should be exact and strand sensitive...
	$pfset->region($seqname) ;
	if(my $count = $pfset->count()) {
	    $pfset = $pfset->rewriteField('GROUP','Sequence', $seqname, undef,'Old_Sequence' ) ;

	    if($statab) {
		# CDS in gene data
		print STDERR "### GFF::Analysis::mergeGeneFeatures() warning: pseudogene sequence '$seqname' has no CDSs?\n" if !$count ;
		if($count) {
		    $Xcp  += $count ;
		    $X2cp += $count*$count ;
		}
		if(!exists($cpGeneC{"$count"})) {
		    $cpGeneC{"$count"} = 1  ;
		} else {
		    $cpGeneC{"$count"}++  ;
		}
	    }
	}
	$pfset = constructGene($pfset,1) unless($raw) ; # add introns et al. w/o a sequence record
        $gffpseudogenes->addGFF($pfset) ;

    } # foreach pseudogene...

    if($statab) {
	# compute the exon and CDS count stats (Ng, Np already set above?)

        # First, for genes...

	$statab->{'Gene'} = {} if(!defined($statab->{'Gene'})) ;
	croak "GFF::Analysis::mergeGeneFeatures(): invalid statistics hash argument: ".($statab->{'Gene'})."\n" 
	    if( ref($statab->{'Gene'}) !~ /HASH/ ) ;

	$statab->{'Gene'}->{'Exon'} = {} if(!defined($statab->{'Gene'}->{'Exon'})) ;
	croak "GFF::Analysis::mergeGeneFeatures(): invalid statistics hash argument: ".($statab->{'Gene'}->{'Exon'})."\n" 
	    if( ref($statab->{'Gene'}->{'Exon'}) !~ /HASH/ ) ;

	if($Ng) {
	    my $Meg = $Xeg/$Ng ;

	    $statab->{'Gene'}->{'Exon'}->{'M'}     = sprintf("%.2f",$Meg) ;
	    $statab->{'Gene'}->{'Exon'}->{'SD'}    = sprintf("%.2f",( $X2eg/$Ng - $Meg**2 )**0.5) ;
	    $statab->{'Gene'}->{'Exon'}->{'X'}     = $Xeg ;
	    $statab->{'Gene'}->{'Exon'}->{'X2'}    = $X2eg ;
	    $statab->{'Gene'}->{'Exon'}->{'GeneC'} = \%egGeneC ;
  
	    $statab->{'Gene'}->{'CDS'} = {} if(!defined($statab->{'Gene'}->{'CDS'})) ;
	    croak "GFF::Analysis::mergeGeneFeatures(): invalid statistics hash argument: ".($statab->{'Gene'}->{'CDS'})."\n" 
		if( ref($statab->{'Gene'}->{'CDS'}) !~ /HASH/ ) ;

	    my $Mcg = $Xcg/$Ng ;
	    $statab->{'Gene'}->{'CDS'}->{'M'}     = sprintf("%.2f",$Mcg) ;
	    $statab->{'Gene'}->{'CDS'}->{'SD'}    = sprintf("%.2f",( $X2cg/$Ng - $Mcg**2 )**0.5) ;
	    $statab->{'Gene'}->{'CDS'}->{'X'}     = $Xcg ;
	    $statab->{'Gene'}->{'CDS'}->{'X2'}    = $X2cg ;
	    $statab->{'Gene'}->{'CDS'}->{'GeneC'} = \%cgGeneC ;

	    my $Mts = $Xts/$Ng ;

	    $statab->{'Gene'}->{'Transcript'}->{'M'}     = sprintf("%.2f",$Mts) ;
	    $statab->{'Gene'}->{'Transcript'}->{'SD'}    = sprintf("%.2f",( $X2ts/$Ng - $Mts**2 )**0.5) ;
	    $statab->{'Gene'}->{'Transcript'}->{'X'}     = $Xts ;
	    $statab->{'Gene'}->{'Transcript'}->{'X2'}    = $X2ts ;
	    $statab->{'Gene'}->{'Transcript'}->{'GeneC'} = \%ts_cov ;
  
	    my $Mcs = $Xcs/$Ng ;
	    $statab->{'Gene'}->{'Translation'}->{'M'}     = sprintf("%.2f",$Mcs) ;
	    $statab->{'Gene'}->{'Translation'}->{'SD'}    = sprintf("%.2f",( $X2cs/$Ng - $Mcs**2 )**0.5) ;
	    $statab->{'Gene'}->{'Translation'}->{'X'}     = $Xcs ;
	    $statab->{'Gene'}->{'Translation'}->{'X2'}    = $X2cs ;
	    $statab->{'Gene'}->{'Translation'}->{'GeneC'} = \%cs_cov ;
	} else {
	    print STDERR "### GFF::Analysis::mergeGeneFeatures() ",
	                 "Warning: Number of true genes encountered is zero!?\n" if defined($trace) and $trace ;
	}

	# then for pseudogenes...

	$statab->{'Pseudogene'} = {} if(!defined($statab->{'Pseudogene'})) ;
	croak "GFF::Analysis::mergeGeneFeatures(): invalid statistics hash argument: ".($statab->{'Pseudogene'})."\n" 
	    if( ref($statab->{'Pseudogene'}) !~ /HASH/ ) ;

	$statab->{'Pseudogene'}->{'Exon'} = {} if(!defined($statab->{'Pseudogene'}->{'Exon'})) ;
	croak "GFF::Analysis::mergeGeneFeatures(): invalid statistics hash argument: ".($statab->{'Pseudogene'}->{'Exon'})."\n" 
	    if( ref($statab->{'Pseudogene'}->{'Exon'}) !~ /HASH/ ) ;

	if($Np) {
	    my $Mep = $Xep/$Np ;

	    $statab->{'Pseudogene'}->{'Exon'}->{'M'}     = sprintf("%.2f",$Mep) ;
	    $statab->{'Pseudogene'}->{'Exon'}->{'SD'}    = sprintf("%.2f",( $X2ep/$Np - $Mep**2 )**0.5) ;
	    $statab->{'Pseudogene'}->{'Exon'}->{'X'}     = $Xep ;
	    $statab->{'Pseudogene'}->{'Exon'}->{'X2'}    = $X2ep ;
	    $statab->{'Pseudogene'}->{'Exon'}->{'GeneC'} = \%epGeneC ;
  
	    $statab->{'Pseudogene'}->{'CDS'} = {} if(!defined($statab->{'Pseudogene'}->{'CDS'})) ;
	    croak "GFF::Analysis::mergeGeneFeatures(): invalid statistics hash argument: ".($statab->{'Pseudogene'}->{'CDS'})."\n" 
	    if( ref($statab->{'Pseudogene'}->{'CDS'}) !~ /HASH/ ) ;

	    my $Mcp = $Xcp/$Np ;
	    $statab->{'Pseudogene'}->{'CDS'}->{'M'}     = sprintf("%.2f",$Mcp) ;
	    $statab->{'Pseudogene'}->{'CDS'}->{'SD'}    = sprintf("%.2f",( $X2cp/$Np - $Mcp**2 )**0.5) ;
	    $statab->{'Pseudogene'}->{'CDS'}->{'X'}     = $Xcp ;
	    $statab->{'Pseudogene'}->{'CDS'}->{'X2'}    = $X2cp ;
	    $statab->{'Pseudogene'}->{'CDS'}->{'GeneC'} = \%cpGeneC ;
	} else {
	    print STDERR "### GFF::Analysis::mergeGeneFeatures() ",
	                 "Warning: Number of pseudogenes encountered is zero!?\n" if defined($trace) and $trace ;
	}
    }

     print STDERR "\tReturning ",$gffgenes->count(), " gene feature sets (sequences, exons, introns & CDSs) and\n",
                   $gffpseudogenes->count()," pseudogene feature sets (sequences, exons, introns & CDSs), respectively\n\n",
                   "Exiting GFF::Analysis::mergeGeneFeatures()\n" if defined($trace) and $trace ;

    GFF::trace(0) if defined($trace) and $trace ;

    return ($gffgenes, $gffpseudogenes) ;
}


=pod

=over 4

=item cleanUpSeqName(name,keepIsoforms)

Default namefilter used by normalize(), which removes:

- CDS/mRNA suffixes
- 5'/3' suffixes
- alphabet or digit isoform suffixes (unless keepIsoforms)
- invalid characters

If the 'keepIsoforms' flag is defined and non-null,
then isoforms are kept.

The cleaned up name is returned.

=back

=cut

sub cleanUpSeqName {
    my $name = shift ;
    my $keepIsoforms = shift ;

    return undef if(!(defined($name) and $name)) ;
    $name =~ s/\;\!"\$\%\^\&\*\(\)\{\}\[\]\@\~\#'\?\<\>\|//g ; # filter out suspicious characters
    # $name =~ s/^em://i ;              # no EMBL 'em:' prefixes - but will removing these cause problems later?
    $name =~ s/(-|_|\.)(CDS|mRNA)$//i ; # no CDS/mRNA suffixes
    $name =~ s/-(5|3)\'$//i ;           # no 5'/3' suffixes
    unless(defined($keepIsoforms) and $keepIsoforms) {
       $name =~ s/\.?([a-z])$//i ;         # no alpha isoforms
       $name =~ s/(\.\d+)\.\d+$/$1/i ;     # no numeric isoforms
    }
    return $name ;
}

=pod

=over 4

=item normalize($gffndg, $gffndp, $gffnde, $gffndc, \&nameFilter, \%statab, $trace)

Method to return a normalized set of gene descriptions in which
all isoforms and exons have been merged into distinct, non-overlapping
non-duplicated sets of data. 

&nameFilter should take a $name string as $input, perform 
clean up of name decorations, then return the $name string.
If nameFilter is undef or NULL, then name clean-up is suppressed.
If nameFilter is defined, non-NULL but not a reference point to 
a subroutine, then a standard cleanup of names is performed.
Otherwise, the user supplied subroutine is used.

=back

=cut

sub normalize {
    my $gffndg = shift ;
    GFF::GeneFeatureSet->verify($gffndg) ;

    my $gffndp = shift ;
    GFF::GeneFeatureSet->verify($gffndp) ;

    my $gffnde = shift ; 
    GFF::GeneFeatureSet->verify($gffnde) ;

    my $gffndc = shift ;
    GFF::GeneFeatureSet->verify($gffndc) ;

    my $namefn = shift ;
    my $statab = shift ;
    my $trace  = shift ;

    GFF::trace(1) if(defined($trace) and $trace) ;

    $namefn = \&cleanUpSeqName if defined $namefn && !ref($namefn) ; # if $namefn is only a flag

    #################################################################
    # Exact, strand sensitive, normalized merge
    # of (separated) transcript sequences ('isoforms') => 'Genes'
    #################################################################
    print STDERR "\tPerform self overlap merge of components...\n" if defined($trace) and $trace ;
    my $gffg = $gffndg->self_overlap_merge(0,1,'Overlaps') ;
    my $gffp = $gffndp->self_overlap_merge(0,1,'Overlaps') ;

    #################################################################
    # Conceptually, the merged transcript sequences are 'genes'
    # so take some more length stats of these
    #################################################################
    if($statab) {
	featureLengthStats($gffg,$statab,'Gene') ;
	featureLengthStats($gffp,$statab,'Pseudogene') ;
    }

    #################################################################
    # Exact, strand sensitive, normalized merge
    # of (separated) exons && CDS's
    #################################################################
    my $gffne  = $gffnde->self_overlap_merge(0,1) ; # don't bother recording exon...
    my $gffnc  = $gffndc->self_overlap_merge(0,1) ; # and CDS overlap merges...

    #####################################################################
    # Use an explicit loop and substitutions for the Sequence names here
    # because of the complexity of the rewrite rule...
    #####################################################################
    if(defined($namefn) && ref($namefn)) {
	print STDERR "\tCleaning up sequence names...\n" if defined($trace) and $trace ;
	foreach my $gf ($gffg->eachGeneFeature()) {
	    my $oldname = $gf->group_value('Sequence') ;
	    next if !(defined($oldname) and $oldname) ;
	    my $newname = &{$namefn}($oldname) ;
	    $gf->group_value('Sequence',0,$newname) ;
	}
	foreach my $gf ($gffp->eachGeneFeature()) {
	    my $oldname = $gf->group_value('Sequence') ;
	    next if !(defined($oldname) and $oldname) ;
	    my $newname = &{$namefn}($oldname) ;
	    $gf->group_value('Sequence',0,$newname) ;
	}
    }

    GFF::trace(0) if defined($trace) and $trace ;

    return ($gffg,$gffp,$gffne,$gffnc) ;

}

=pod

=over 4

=item normalize_mRNA($gffi, \&nameFilter, $source, \%statab, $trace )

Invokes segregate, normalize and merge routines (see above) to 
normalize transcript GFF.

The optional $source argument is used to rewrite the <source> field
of the file to a uniform value (default: 'Gene').

Providing a defined reference to an empty hash, '\%statab' triggers
the compilation of statistics about the file, as generated in the
segregateGeneFeatures() and mergeGeneFeatures() methods (see above). 

Separate statistics are generated for each of transcripts, genes,
exons  and CDS's. (Note that 'genes' are defined as the
normalized transcript sequence spans output by the method, whereas
'transcripts' are the sequence records before merge overlaps are done).
Note that the gene sequence count 'N' are after normalization, but all other
feature counts 'N' are unnormalized numbers.

A defined and non-null $trace flag turns on runtime tracing of
the normalization method.
 
=back

=cut

sub normalize_mRNA {

    my $gffi   = shift ;
    GFF::GeneFeatureSet->verify($gffi) ;

    my $namefn = shift ;
    my $source = shift ; # labelling for features

    # if defined, should be reference to an empty %hash 
    # within which the method will place some computed statistics
    my $statab = shift ; 

    my $trace  = shift ;

    $source = 'Gene'           if !defined($source) ;
    $statab = 0                if !defined($statab) ;

    croak "GFF::Analysis::normalize_mRNA(): invalid statistics hash reference argument: $statab\n" 
	if $statab && ref($statab) !~ /HASH/ ;

    my ($gffndg, $gffndp, $gffnde, $gffndc) = segregateGeneFeatures($gffi,$statab,$trace) ;

    my ($gffg,$gffp,$gffne,$gffnc) = normalize($gffndg, $gffndp, $gffnde, $gffndc, $namefn, $statab,$trace) ;

    my ($gffgenes,$gffpseudogenes) = mergeGeneFeatures($gffg,$gffp,$gffne,$gffnc,$statab,$trace) ;

    ######################################################
    #  Rewrite the true transcript set <source> field, 
    #  to take into account the various ways 
    #  true transcripts are labelled...
    ######################################################
    print STDERR "\tRelabel transcript set...\n" if defined($trace) and $trace ;
    my $gffgset = $gffgenes->rewriteField('SOURCE','*',$source, 0,'Source') ;
    my $gffpset = $gffpseudogenes->rewriteField('SOURCE','*','Pseudogene',0,'Source') ;

    print STDERR "Exiting normalize_mRNA()\n\n" if defined($trace) and $trace ;

    GFF::trace(0) if defined($trace) and $trace ;

    return ($gffgset,$gffpset);

}

=pod

=over 4

=item gfProximity($gff1, $gff2, $source, $feature, $t55, $t53, $t35, $t33, $incr, \%statab, $label )

This routine annotates one GFF feature set ($gffi) with proximity/overlap 
hits on a second GFF feature set (e.g. genes with CpG Islands) 
and concurrent computes related statistics.

$source and $feature are used to filter out those matches specifically 
annotated by the current overlap comparison (i.e. where several different 
matches may be annotated in $gff1). The 'source' field is also used to 
index the \%statab entries, unless the optional '$label' argument is given
for this purpose.

$t55, $t53, $t35, $t33 are the basepair overlap threshold boundary conditions 
as passed to the GFF::GeneFeatureSet::intersect_overlap_matches() method
as the anonymous array '[$t55, $t53, $t35, $t33]'.  These all default to zero.
If they are all zero, then gene feature match() is passed '-1' tolerance, that is,
accepts all overlaps.

The $binsize argument specifies the bin size (in base pairs) and defaults to 1 kb.

A defined reference to a hash, '\%statab' must be provided for
the compilation of statistics about the overlaps (the method is
rather pointless otherwise!). The overlaps are computed in $binsize kilobase
increments between and inclusive of all the boundary thresholds indicated,
for the 5' and 3' distances. Negative distance index values indicate that
the $gff2 position is upstream (5') of $gff1; positive values, downstream.

The statistics are returned as:

         $statab->{$source}->{(5'|3')}->{$distance}

where $distance ranges over the specified $t* tolerance values.
 
=back

=cut

sub max{
    my($i,$j)=@_;
    if($i>$j){
	$i;
    }else{
	$j;
    }
}

sub gfProximity { 
    my $gff1    = shift ;
    GFF::GeneFeatureSet->verify($gff1) ;

    my $gff2    = shift ;
    GFF::GeneFeatureSet->verify($gff2) ;

    my $source  = shift ;
    my $feature = shift ;
    my $t55     = shift ;
    my $t53     = shift ;
    my $t35     = shift ;
    my $t33     = shift ;
    my $binsize = shift ;
    my $statab  = shift ;
    my $label   = shift ;
    $label = $source if !defined($label) ;

    croak "GFF::Analysis::gfProximity(): at least one overlap tolerance must be defined!" 
	if !(defined($t55) or defined($t53) or defined($t35) or defined($t33)) ;

    $t55 = 0 if !defined($t55) ;
    $t53 = 0 if !defined($t53) ;
    $t35 = 0 if !defined($t35) ;
    $t33 = 0 if !defined($t33) ;
    my ($distbounds,$b5,$b3) ;
    if($t55 or $t53 or $t35 or $t33) {
	$distbounds = [$t55,$t53,$t35,$t33] ;
	# maximum bounds of 5' and 3' tolerance values?
	$b5=max($t55,$t35) ;
	$b3=max($t53,$t33) ;
    } else {
	$distbounds = -1 ;
	$b5=$b3=1 ;
    }

    $binsize=1000 if !defined($binsize) ;

    croak "GFF::Analysis::gfProximity() needs a statistics table hash reference!" 
	if !(defined($statab) and ref($statab) =~ /HASH/) ;

    my $gffm = $gff1->intersect_overlap_matches( $gff2, $distbounds ) ;

    $statab->{$label} = {} ;
    for(my $dist = -$b5 ; $dist <= $b3; ++$dist) {
	$statab->{$label}->{"5'"}->{"$dist"} = 0 ;
	$statab->{$label}->{"3'"}->{"$dist"} = 0 ;
    }
    $statab->{$label}->{'Total'} = $gff2->count() ; # total hits

    foreach my $gf ($gffm->eachGeneFeature()) {
	# getMatches() returns a reference to a hash keyed on $gff2 features
	foreach my $match (values %{$gf->getMatches($source,$feature)}) {
	    my ($hit,$e5,$e3) = @{$match} ;
	    if(my $dist5 = int($e5/$binsize)) {
		if( $gf->strand() eq '+' ) {
		    $statab->{$label}->{"5'"}->{$dist5}++ ;
		} else { # reverse stranded gene...
		    $statab->{$label}->{"5'"}->{"-$dist5"}++ ;
		}
	    } else {
		$statab->{$label}->{"5'"}->{"0"}++ ;
	    }
	    if(my $dist3 = int($e3/$binsize)) {
		if( $gf->strand() eq '+' ) {
		    $statab->{$label}->{"3'"}->{$dist3}++ ;
		} else { # reverse stranded gene...
		    $statab->{$label}->{"3'"}->{"-$dist3"}++ ;
		}
	    } else {
		$statab->{$label}->{"3'"}->{"0"}++ ;
	    }
	}
    }
}

1;  # says use was ok

__END__

=head1 Revision History

2.16 (25/05/2000) rbsk: $offset argument added to mRNA() method

2.15 (31/03/2000) rbsk: featureLengthStats() now also returns min and max length stats.

2.14 (10/01/2000) rbsk: copy over '##sequences' in makeGenes

2.13 (03/01/2000) rbsk: treat 'pseudogene' as mRNA in redundancy filter 
                        in segregateGeneFeatures() to avoid CDS gene feature duplication

2.12 (29/12/99) - rbsk: mergeGeneFeatures() logical error: constructGene in wrong place
                        added '$raw' argument to suppress constructGene invocation here

2.11 (15/12/99) - rbsk: small bug in constructGene() method

2.10 (19/11/99) - rbsk: added the gfProximity() method

2.09 (19/10/99) - rbsk: segregateGeneFeatures() and mergeGeneFeatures() handle
                        pseudogenes separately and explicitly;

2.08 (16/10/99) - rbsk: for greater flexibility, I extracted out normalize_mRNA()  
                        functionality, into externally visible methods:

                        - featureLengthStats()
                        - segregateGeneFeatures()
                        - mergeGeneFeatures()

2.07 (13/10/99) - rbsk: added $statab argument to normalize_mRNA()
2.06 (9/10/99)  - rbsk: now exporting cleanUpSeqName(); keep 'em:' prefixes for now...
2.05 (27/9/99)  - rbsk: need to make all normalization procedures strand sensitive!
2.04 (21/7/99)  - rbsk: normalize_mRNA() should not consider source(*CDS) records
                        with feature(sequence) to be redundant, in case the specific
                        'sequence' only has a CDS specified (but no mRNA);
                        This method now also relabels <source> fields to 'TranscriptSet'
2.03 (16/7/99)  - rbsk: removed draw_graph() from here (into GFF::Graph()) because
                        of Curve_plot.pm usage, which won't be universal outside Sanger
                        (at least until Raphael decides to release it for general use?)
2.02 (14/7/99)  - rbsk: transferred draw_graph() from GeneFeatureSet to here
                        and generalized to multiple GFF plot (API changed)
2.01 (12/7/99)  - rbsk: creation from miscellaneous GFF analysis code;
                        transferred methods makeGenes(), constructGene()
                        and mRNA from GFF::GeneFeatureSet to this module
