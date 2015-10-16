####################################################################
#
# HomolGeneFeature.pm - Perl Object GFF Homology Gene Feature objects
#
# Copyright (c) 1999 
# Created by   Tim Hubbard <th@sanger.ac.uk>
# Augmented by Richard Bruskiewich <rbsk@sanger.ac.uk>
#
# Sanger Centre, Wellcome Trust Genome Campus, Cambs, UK
# All rights reserved.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation
#
# History:
#   18/8/99 - rbsk: coded explicit primary field access methods, rather than 
#                   relying upon AUTOLOAD (i.e. to gain efficiency)
#   28/4/99 - rbsk: HomolGeneFeature.pm => GFF::HomolGeneFeature.pm
#   22/2/99 - rbsk: added Version 2 GFF code including &version class function
#
####################################################################
package GFF::HomolGeneFeature;
use Carp;
use strict;
use vars qw($AUTOLOAD @ISA @EXPORT @EXPORT_OK $VERSION);
require Exporter;

$VERSION = '3.004' ;
#
# @ISA has our inheritance.
#
@ISA = qw( Exporter GFF::GeneFeature );
#
# Place functions/variables you want to *export*, ie be visible from the caller package into @EXPORT_OK
#
@EXPORT    = qw();
@EXPORT_OK = qw();

=pod

##### old implementation #########

# need these for the permitted thing..

my %fields = (
              version  => undef, # gene feature objects need to know what GFF Version they are?
	      seqname  => undef, 
	      source   => undef,
	      feature  => undef,
	      start    => undef,
	      end      => undef,
	      score    => undef,
	      strand   => undef,
	      frame    => undef,
	      group    => undef,
              comment  => undef, 
	      matches  => undef,
	      member   => undef,
	      start2    => undef,
	      percentid => undef,
	      end2      => undef,
              exp_value => undef,  # Version 2 GFF
);

=cut

################### 02/02/2000 attempt at increasing memory efficiency #####################

# Refer to GFF::GeneFeature.pm too for concept

# Fixed length array indices
use constant MATCHES => 11 ; # this value should be the same as GFF::GeneFeature::MATCHES
use constant MEMBER  => 12 ; # this value should be the same as GFF::GeneFeature::MEMBER
use constant HSTART  => 14 ; # this value should be GFF::GeneFeature::COMMENT+1
use constant HEND    => 15 ;
use constant PCID    => 16 ;
use constant EVALUE  => 17 ;

################### 02/02/2000 attempt at increasing memory efficiency #####################


sub new {
    my $ref = shift;
    my $class = ref($ref) || $ref;
    my $version = shift ;
    
    if(ref($ref) && GFF->verify($ref)) {
        $version = $ref->version() if !defined $version ;
    }

    my $self = GFF::GeneFeature->new($version) ;

    bless $self, $class ;
    return $self;
}

###################################################################
# Access methods - more efficient to be explicit than AUTOLOAD'ed #
# See also inheritance from GeneFeature.pm ########################
###################################################################

#
# Sets and/or returns homology target
# GFF Version 1: just a group name? Tim coded things this way in Version 1
#
sub target {
    my $self   = shift ;
    my $class  = ref($self) || carp "\'$self\' is not an object of mine!";
    my $target = shift ;
    if(defined($target)) {
	if( $self->version() == 1 ) {
	    return $self->group($target) ; 
	} else {
	    return $self->group_value('Target',0,$target) ;
	}
    } else {
	if( $self->version() == 1 ) {
	    return $self->group() ; 
	} else {
	    return $self->group_value('Target') ;
	}
    }
}

sub start2 {
    my $self = shift;
    my $class = ref($self) || carp "\'$self\' is not an object of mine!";
    my $start2 = shift ;
    if (defined($start2)) {
	if($start2 !~ /^[\+-]?\d+$/) {
	    carp "GFF::HomolGeneFeature::start2 must be an integer\n" ;
	    return undef ;
	} else {
	    if( $self->version() == 1 ) {
		return $self->[HSTART] = $start2;
	    } else {
		return $self->group_value('Target',1,$start2);
	    }
	}
    } else { 
	if( $self->version() == 1 ) {
	    return $self->[HSTART];
	} else {
	    return $self->group_value('Target',1);
	}
    }
}

sub end2 {
    my $self = shift;
    my $class = ref($self) || carp "\'$self\' is not an object of mine!";
    my $end2 = shift ;
    if (defined($end2)) {
	if($end2 !~ /^[\+-]?\d+$/) {
	    carp "GFF::HomolGeneFeature::end2 must be an integer\n" ;
	    return undef ;
	} else {
	    if( $self->version() == 1 ) {
		return $self->[HEND] = $end2;
	    } else {
		return $self->group_value('Target',2,$end2);
	    }
	}
    } else {
	if( $self->version() == 1 ) {
	    return $self->[HEND];
	} else {
	    return $self->group_value('Target',2);
	}
    }
}

sub percentid {
    my $self = shift;
    my $class = ref($self) || carp "\'$self\' is not an object of mine!";
    my $pcid=shift ;
    if (defined($pcid)) {
	if($pcid !~ /^[\+-]?(\d+(\.\d+)?|\.)$/) {
	    carp "GFF::HomolGeneFeature::percentid must be a integer/float number or '.'\n" ;
	    return undef ;
	} else {
	    if( $self->version() == 1 ) {
		return $self->[PCID] = $pcid;
	    } else {
		return $self->group_value('Percent_Id',0,$pcid);
	    }
	}
    } else {
	if( $self->version() == 1 ) {
	    return $self->[PCID];
	} else {
	    return $self->group_value('Percent_Id');
	}
    }
}

sub exp_value {
    my $self = shift;
    my $class = ref($self) || carp "\'$self\' is not an object of mine!";
    my $evalue=shift ;
    if (defined($evalue)) {
	if($evalue !~ /^[\+-]?(\d+(\.\d+)?|\.)$/) {
	    carp "GFF::HomolGeneFeature::percentid must be a integer/float number or '.'\n" ;
	    return undef ;
	} else {
	    if( $self->version() == 1 ) {
		return $self->[EVALUE] = $evalue;
	    } else {
		return $self->group_value('E_Value',0,$evalue);
	    }
	}
    } else {
	if( $self->version() == 1 ) {
	    return $self->[EVALUE];
	} else {
	    return $self->group_value('E_Value');
	}
    }
}

sub fromGeneFeature {
    my $self    = shift ;

    # Version 2 GFF can get these parameters explicitly from Target tag-values
    my $target  = shift ; # Version 1/2; may be null in Version 2
    my $start2  = shift ; # Version 1 parameter
    my $end2    = shift ; # Version 1 parameter
    my $pcid    = shift ;
    my $evalue  = shift ;

    bless $self,'GFF::HomolGeneFeature' ;

    if( $self->version() == 1 ) {
        $self->group($target) ;  # also just a group name? Tim coded things this way in Version 1
        $self->[HSTART] = $start2;
        $self->[HEND]   = $end2;
        $self->[PCID]   = $pcid if defined($pcid) ;
        $self->[EVALUE] = $evalue if defined($evalue);
    } else {
        # assume a Version 2 tag-value specification of the homology feature
        # generally: "Target <$target> <$start2> <$end2>; Percent_Id <$pcid>; E_Value <$exp_value>"
        if( defined $target and $target ) { 
	     $self->parse_group($target) ; 
        }
    }
    return $self ;
}

#
# Class method to parse a tab delimited line into a GeneFeature object
#

sub new_from_line {
    my $ref     = shift;
    my $class   = ref($ref) || $ref;
    my $string  = shift ;
    my $version = shift ;

    my $self = HomolGeneFeature->SUPER::new_from_line( $string, $version ) ;

    $self->fromGeneFeature() ; # convert $self GeneFeature reference into Homol...

    return $self ;
}

#
# Method to parse a line from an MSP file from MSPcrunch -d output
#
sub new_from_msp {
    my $ref      = shift ;
    my $class    = ref($ref) || $ref ;   
    my $string   = shift ;
    my $source   = shift ;
    my $function = shift ;
    my $version  = shift ;

    my ($self,$group);
    my @array;

    $self = GFF::HomolGeneFeature->new( $version ) ; 

    # predefined
    $self->source($source);
    $self->feature('similarity');
    # undefined
    $self->strand('.'); # rather than 0
    $self->frame('.');  # rather than 0

    chomp($string);

    # remove leading spaces
    $string=~s/^\s+//;

    @array = split(/\s+/,$string);

    # load in order
    # raw score
    $self->score(shift(@array));
    # ignore %id
    $self->percentid(shift(@array));

    # FIXME
    # save start-end start2-end2 sensibly
    # it is possible that strand should be set if the object
    # matching does in reverse and if of known direction (protein)
    # but this is not done here.

    my $start=shift(@array);
    my $end=shift(@array);
    $self->seqname(shift(@array));
    my $start2=shift(@array);
    my $end2=shift(@array);
    
    # invert
    if($start>$end){
	my $t;
	$t=$start;
	$start=$end;
	$end=$t;
	$t=$start2;
	$start2=$end2;
	$end2=$t;
    }

    $self->start($start);
    $self->end($end);
    # allow special rule to be used to parse name
    $group=shift(@array) ;
    if($function){
	$group = &$function( $group, \@array ) ;
    }

    if( $self->version() == 1 ) {
        $self->group($group);
    } else {
        $self->group_value('Target',0,$group);
        $self->group_value('Description',0,join(' ',@array));
    }

    # I set these after 'Target'->[0] now
    # note that this means that start2>end2 is possible
    $self->start2($start2);
    $self->end2($end2);

    return $self;
}

#
# Method to dump a feature as a line of MSPcrunch -d output
#
sub dump_msp {
    my $self     = shift;
    my $file     = shift;

    if( !defined $file ){
	$file = \*STDOUT;
    }

    my @array;

    # raw score
    my $score=$self->score();
    croak "HomolGeneFeature::dump_msp: missing score?" if !$score;
    push(@array,$score) if $score;

    # ignore %id
    my $pid=$self->percentid();
    croak "HomolGeneFeature::dump_msp: missing percent id?" if !$pid;
    push(@array,$pid) if $pid;

    my $start=$self->start(); 
    croak "HomolGeneFeature::dump_msp: missing query start?" if !$start;
    push(@array,$start);

    my $end=$self->end();
    croak "HomolGeneFeature::dump_msp: missing query end?" if !$end;
    push(@array,$end) if $end;

    my $seqname=$self->seqname();
    croak "HomolGeneFeature::dump_msp: missing seqname?" if !$seqname;
    push(@array,$seqname) if $seqname;

    my $start2=$self->start2();
    croak "HomolGeneFeature::dump_msp: missing target start?" if !$start2;
    push(@array,$start2) if $start2;

    my $end2=$self->end2();
    croak "HomolGeneFeature::dump_msp: missing target end?" if !$end2;
    push(@array,$end2) if $end2;

    if( $self->version() == 1 ) {

	my $group=$self->group();
	$group='?' if !defined($group);
        push(@array,$group);

    } else {

	my $target=$self->group_value('Target');
	croak "HomolGeneFeature::dump_msp: missing target name?" if !$target;
	push(@array,$target);

	my $title=$self->group_value('Title');
	croak "HomolGeneFeature::dump_msp: missing description?" if !$title;
	push(@array,$title);
    }

    print $file join(' ',@array)."\n";
}

# for reading a line parsed by a function
sub new_from_parse {
    my $ref      = shift ;
    my $class    = ref($ref) || $ref ;
    my $string   = shift ;
    my $function = shift ;
    my $version  = shift ;

    my $self ;

    $self = new GFF::HomolGeneFeature( $version ) ; 

# if function 
    if($function){
	if(&$function($self,$string)){
	    return $self;
	}else{
	    return '';
	}
    }else{
	croak "function required by new_from_parse in HomolGeneFeature";
    }
}

#
# Method to make a new gf from an old one
#

sub copy {
    my $self    = shift ;
    my $version = shift ;

    my $other ;

    $other = new GFF::HomolGeneFeature( $version ) ; 

    $other->seqname($self->seqname);
    $other->source($self->source);
    $other->feature($self->feature);
    $other->start($self->start);
    $other->end($self->end);
    $other->score($self->score);
    $other->strand($self->strand);
    $other->frame($self->frame);

    $other->start2($self->start2);
    $other->percentid($self->percentid);
    $other->end2($self->end2);

    %{$other->[MATCHES]} = ( %{$self->[MATCHES]} );
    %{$other->[MEMBER]} = ( %{$self->[MEMBER]} ) ; 

    if( $self->version() == 1 ) {
        $other->group($self->group()) ;
    } else { # Version 2 or greater has tag-value matches to copy over?
        my $tag ;
        foreach $tag (keys %{$self->group()}) {
            my $values = [] ;

            # $self->group()->{$tag} is assumed to be an array reference
            push @{$values}, @{$self->group()->{$tag}} ;
  
	    $other->group()->{$tag} = \@{$values} ;
	}
    }    
    return $other;
}

#
# Method to return a formatted output of a GeneFeature object to a filehandle.
# Uses dump_string
#

sub dump {
    my $self = shift;
    my $file = shift;
    my $tab = shift;
    
    if( !defined $file ){
	$file = \*STDOUT;
    }

    print $file $self->dump_string($tab)."\n";
}

#
# Method to return a formatted output of a HomolGeneFeature object to a string
#

sub dump_string {
    my $self    = shift;
    my $tab     = shift ;
    my $newline = shift ;
    my $tag     = shift ;

    my $string = $self->SUPER::dump_string($tab,$newline,$tag) ;

    # GFF version 2 are dumped as tag-values
    # so don't need special treatment
    if( $self->version() == 1 ) {
	$string .= sprintf("%s\t%s\t%5.2f\t%f",
                           $self->[HSTART],
	                   $self->[HSTART],
	                   $self->[PCID],
                           $self->[EVALUE]);
    }
}

#
# function to determine the gap between 2 gf objects
# assumptions:
# - order is self(start-end)...other(start-end)
#

sub gap{
    my $self  = shift ;
    my $other = shift ;
    my $warn  = shift ;

    my($unit,$type,$mgap,$ngap,$gap,$gf,$frameshift,$ed,$st);
    # check both are same match type
    $type=$self->source;
    if($type ne $other->source){
	croak "Fatal error in HomolGeneFeature::gap - different homol types: ".
	    $type.",".$other->source;
    }
    # check valid match type
    if($type eq 'tblastn'){
	$unit=3;
    }else{
	croak "Fatal error in HomolGeneFeature::gap - type $type not recognised";
    }

    # gap in master seq
    $mgap=($other->start)-($self->end)-1 ;

    # gap in match
    # check orientation
    if(($self->end2)<($self->start2)){
	# reverse match -> expect reverse order
	if(($self->start2)<($other->start2)){
	    print "Warn: Matches are not ordered - ignored: ".
		$self->seqname.': '.$self->start2.','.$self->end2.
		    ','.$other->start2.','.$other->end2."\n";
	}else{
	    $ngap=($other->start2)-($self->end2)-1;
	}
    }else{
	if(($self->start2)>($other->start2)){
	    print "Warn: Matches are not ordered - ignored: ".
		$self->seqname.': '.$self->start2.','.$self->end2.
		    ','.$other->start2.','.$other->end2."\n";
	}else{
	    $ngap=($other->start2)-($self->end2)-1;
	}
    }
    if($unit==3){
	if($frameshift=$ngap%3){
	    if($warn){
		print "Warn: Match offset not divisible by 3 ($ngap): "
		    .$self->seqname.' '.$self->group_value('Target')."\n";
	    }
	    # round towards zero
	    if($ngap>0){
		$ngap=$ngap-$ngap%3;
	    }else{
		$ngap=$ngap-$ngap%3+3;
	    }
	}
	$ngap=$ngap/3;
    }
    # calculate gap and return if 0
    $gap=$ngap-$mgap;
    return unless $gap;
    # real indels
    if($self->percentid>90.0 && $gap>2 && $warn){
	print "$gap $ngap $mgap ".$self->seqname.' '.$self->group_value('Target').' ('.$self->percentid.")\n";
	print $self->start.'-'.$self->end.','.$self->start2.'-'.$self->end2.
	    ' '.$other->start.'-'.$other->end.','.$other->start2.'-'.$other->end2."\n";
    }
    # gap is represented as
    # length 1 feature (-ve score) for insertions
    # length N features (+ve score) for deletions
    # have to extrapolated midpoint
    $st=($self->end)+1;
    $ed=($other->start)-1;
    $st=$ed=int(($st+$ed)/2);
    if($gap>1){
	$ed+=($gap-1)/2;
	if(($gap-1)%2){
	    $ed++;
	}
	$st-=($gap-1)/2;
    }
    
    # now make a gf and return it
    # (not really a correct use of a hgf)
    $gf = new GFF::HomolGeneFeature( $self->version() ) ;

    $gf->seqname($self->seqname);
    $gf->source($self->source);
    $gf->feature('indel');
    $gf->start($st);
    $gf->end($ed);
    $gf->score($gap);
    $gf->strand("."); # rather than 0
    $gf->frame($frameshift);
    if( $self->version() == 1 ) {
        $gf->group($self->group());
    } else {
        my (@values) = @{$self->group_value_list('Target')} ; # makes a copy?
        $gf->group_value_list('Target', \@values ) ;
    }
    $gf->percentid($self->percentid);
    return $gf;
}

1;  # says use was ok
__END__

=head1 NAME

GFF::HomolGeneFeature - Perl extension for GFF Homology Gene Features

=head1 SYNOPSIS

    use GFF ;  # performs an implicit 'use GFF::HomolGeneFeature;'

=head1 AUTHORS

Richard Bruskiewich email - B<rbsk@sanger.ac.uk>

Tim Hubbard email - B<th@sanger.ac.uk>

=head1 DESCRIPTION

GFF::HomolGeneFeature (derived from GFF) is a class of Perl object,
derived (subclassed) from the GFF::GeneFeature object class, that
encapsulates a single "homology" gene feature record ("line") in the
Gene Finding Feature ("GFF") format. A GFF::GeneFeatureSet Perl object
is the container object for a set of GFF::HomolGeneFeature objects.

B<How to Read Method Protocols>

Normal Perl data type notations are used for argument declarations in
the method protocols. A backslash denotes argument passing by
reference.  Arguments shown in italics are optional parameters to the
method calls. Class methods are invoked using the
'class->method(args)' or 'method class args' Perl call formats.

=head1 SOURCE CODE

The most current release of the Perl source code for this module is
available here. NOTE: This is a March 1999 beta release of the modules
incorporating the proposed Version 2 revisions to the GFF
specification. All bug reports may be submitted to Richard Bruskiewich
(rbsk@sanger.ac.uk). Future releases will likely be issued as Bioperl
archived modules.

=head1 GFF::HomolGeneFeature Construction Methods

For all construction methods, an optional "$version" argument may be
given which sets the created object to a specified GFF specification
version. If this argument is not given, then the GFF class (package)
default version is used (see GFF->version()).

=over 4

=item new( $version )

Class method to construct a new, empty GFF::HomolGeneFeature object.

=item new_from_line( $line, $version )

Class method to parse a "homology" gene feature $line string into a
GFF::HomolGeneFeature object (creates and returns the object
reference).

=item new_from_msp( $line, $source, \$name_parse, $version )

Class method to parse a line into a GFF::HomolGeneFeature from
MSPcrunch format (creates and returns the object reference). The
$source string argument is used to fill the GFF source field with a
tag name. The optional &name_parser argument is user-defined function
which takes two arguments $group and \@array, where @array is assumed
to be the (delimiter split) remainder of a line of data beyond the
group field.

=item new_from_parse( $line, \&parser, $version )

Class method to parse the $line string into the GFF::HomolGeneFeature
object using a user defined &parser function (creates and returns the
object reference). The &parser should expect the (empty)
GFF::HomolGeneFeature object reference as its first argument and the
input line (string) as its second argument. So given, the function
should perform the appropriate parsing of the input $line to load the
GFF::HomolGeneFeature object with data.

=item copy( $version )

Method to duplicate the invoking GFF::HomolGeneFeature object.

=item fromGeneFeature( $version (See below for other arguments) )

A "type cast" method to convert an GeneFeature object into a
GFF::HomolGeneFeature object.

=over 4

=item Version 1 GFF protocol

C<fromGeneFeature( 1, $target, $start2, $end2 )> where he $target
string is the [group] field name of the homologous sequence. The start
coordinate of the match ("$start2") and end coordinate of the match
("$end2") should also be given.

=item Version 2 GFF protocol or better

C<fromGeneFeature( 2, $group )> where if the (optional) $group string
argument is given, it is parsed as a string containing a Version 2
style semicolon delimited tag-value pair [group] field description of
the homology match. The recommended Version 2 format for such group
fields is described in the GFF homology feature specification. If the
$group argument is omitted or null, then the group() hash from the
GeneFeature object is dereferenced for 'Target' and 'E_value' tags
with associated GFF::HomolGeneFeature specific values.

=back

=head1 GFF::HomolGeneFeature Output Methods

=over 4

=item dump( \*OUTPUT, $tab, $flen )

Method uses dump_string() to write a formatted output of a GeneFeature
object to a filehandle, OUTPUT. If \*OUTPUT is not given, \*STDOUT is
used. The "$tab" argument is a boolean flag, where a non-null value
directs the use tab as the field delimiter in the output line;
otherwise, use blank space (flag is assumed null if not
specified). The "$flen" argument is a boolean flag, where a non-null
value stipulates that the length of the current output line should be
printed as an extra field at the end of the output line (assumed null
if not specified. Note: the extra length of this field is *not* added
to the displayed line size, but the extra field is tab delimited, if
$tab is set).

=item dump_string( $tab )

Method to return a formatted output of a GFF::HomolGeneFeature object
to a string. The optional "$tab" argument is a boolean flag, where a
non-null value directs the use tab as the field delimiter in the
output line; otherwise, use blank space (flag is assumed null if not
specified).

=back

=head1 GFF::HomolGeneFeature Access Methods

The various GFF record fields may be set or queried by (inherited)
GeneFeature access methods. Additional GFF::HomolGeneFeature fields
are available by the following access methods.  All the methods can
take a single string argument to set the variable. With or without an
argument, the methods return the current (or newly set) value, as a
string, except as specifically noted below:

=over 4

=item target()

(Version 2) homology target (is simply [group] name in Version 1 GFF).

=item start2()

start coordinate of the homologous match sequence.

=item end2()

end coordinate of the homologous match sequence.

=item percentid()

percentage identity between query and match sequence.

=item exp_value()

(Version 2) expectation value.

=back

=head1 GFF::HomolGeneFeature Analysis Methods

=over 4

=item gap( $HGF2, $warn )

Method to determine the gap between the invoking and a second
("$HGF2") GFF::HomolGeneFeature object. Note: The method assumes that
the invoking GFF::HomolGeneFeature object lies upstream (to the left
of) the other object. If the optional "$warn" flag is defined and set
to non-null, then the method prints additional warning diagnostics
onto STDOUT, e.g.  possible frameshifts and indel's between the two
GFF::HomolGeneFeature objects. Returns a reference to a newly created
GFF::HomolGeneFeature object representing the gap.

=back

=head1 Revision History

 3.004 13/7/2000 - th reported patch: dump_msp should croak if data is missing

 3.003 22/06/2000 - rbsk

 'new_from_msp': set 'Target' name first, then start2, $end2 values

 3.002 15/3/2000 - rbsk

 I finally got around to fixing the new memory model for HomolGeneFeatures...

 3.001 2/2/2000 - rbsk

 Memory consumption is out of this world for large feature sets,
 so I devised an novel implementation of GFF::GeneFeatures
                         
 29/11/99 - th

 Fixed bugs in new_from_msp function.  \$gf->start now always before
 \$gf->end but \$gf->start2 and \$gf->end2 may be reversed.

 4/5/99 - rbsk

 Renamed HomolGeneFeature.pm => GFF::HomolGeneFeature.pm 

 3/3/99 - rbsk

 extensively revised and improved the documentation added Version 2 GFF
 code, especially group() field management methods

=head1 SEE ALSO

L<GFF>, L<GFF::GeneFeature> and L<GFF::GeneFeatureSet>

