=head1 NAME

GFF::GeneFeatureSet.pm - Perl extension for GFF (Homol)GeneFeature Set Container

=head1 SYNOPSIS

use GFF ;    # contains an implicit 'use GFF::GeneFeatureSet ;'

=head1 AUTHORS

Richard Bruskiewich email - B<rbsk@sanger.ac.uk>

Tim Hubbard email - B<th@sanger.ac.uk>

Copyright (c) 1999 
Created by   Tim Hubbard <th@sanger.ac.uk>
Augmented by Richard Bruskiewich <rbsk@sanger.ac.uk>

Sanger Centre, Wellcome Trust Genome Campus, Cambs, UK
All rights reserved.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation

=head1 DESCRIPTION

GFF::GeneFeatureSet (derived from GFF) is a Perl
Object for Gene Finding Feature format. A
GFF::GeneFeatureSet object is a container object
for a set of GeneFeature (or HomolGeneFeature)
objects.

=head2 How to Read Method Protocols

Normal Perl data type notations are used for
argument declarations in the method protocols. A
backslash denotes argument passing by reference.
Arguments shown in italics are optional
parameters to the method calls. Class methods are
invoked using the 'class->method(args)' or
'method class args' Perl call formats.

=head1 SOURCE CODE

The most current release of the Perl source code 
for this module is available here.
All bug reports may be submitted to Richard 
Bruskiewich (rbsk@sanger.ac.uk). Future releases will likely
be issued as Bioperl archived modules.

=cut

package GFF::GeneFeatureSet ;

use Carp;
use strict;
use vars qw($AUTOLOAD @ISA @EXPORT @EXPORT_OK $VERSION);
require Exporter;

$VERSION = '2.101' ;
#
# @ISA has our inheritance.
#
@ISA = qw( Exporter GFF );
#
# Place functions/variables you want to *export*, ie be visible from the caller package into @EXPORT_OK
#
@EXPORT    = qw();
@EXPORT_OK = qw();
#
# @ISA has our inheritance.
#

@ISA = ( 'Exporter', 'GFF' );

my %fields = (
              version   => undef, # GFF version of the object
              date      => undef, # date of creation of GFF::GeneFeatureSet
              region    => undef, # sequence-region meta-comment data
              source_version => undef, # source-version meta-comment data
              comments  => undef, # general comments, conditionally echoed to output by dump_header
	      feature   => undef, # going to be an array of GFF::GeneFeatures
	      member    => undef,
	      );

=head1 GFF::GeneFeatureSet Construction Methods

=over 4

=item new( $version, $seqname, $start, $end, \$sequences, $format )

Class method to construct a new empty
GFF::GeneFeatureSet object of version "$version".
If $version is not specified, it is taken to be
the current default GFF version.

Optional values for $seqname, $start, $end,
\$sequences and $format may also be provided.

=back

=cut

sub new {
    my $ref     = shift;
    my $class   = ref($ref) || $ref ;
    my $version = shift ;
    my $name    = shift ;
    my $start   = shift ;
    my $end     = shift ;
    my $sequences = shift ;
    my $format    = shift ;

    my $self = GFF->new($version) ;

    $self->{'_permitted'} = \%fields ;

    my ($year, $month, $day) = (localtime)[5,4,3] ;
    $self->{'date'}    = [ $year+1900, $month+1, $day ] ;

    $self->{'region'}  = []; # makes region a reference to an anonymous array.

    $self->{'feature'} = []; # makes feature a reference to an anonymous array.
    $self->{'member'}  = {}; # makes member a reference to an anonymous array.

    bless $self, $class;

    $self->region($name,$start,$end) ;

    return $self;
}

=pod

=over 4

=item addGeneFeature( $GeneFeature, \&filter, $copy, $exclude )

Method to add a reference to a GeneFeature object
(first argument) to a GFF::GeneFeatureSet object,
possibly subject to an optional, user-defined
&filter function. This predicate (boolean)
function "&filter" tests the GeneFeature object
(given as an argument) for inclusion in the
GFF::GeneFeatureSet object set, based upon user
criteria. If a "&filter" is not provided, then
the GeneFeature is unconditionally included. 

If the $copy_version argument is defined and non-null, 
then the method adds a copy (not the original object references) 
of the source GeneFeature objects to the invoking 
object.

Setting the optional '$exclude' flag complements 
the outcome of the $function discriminator.

=back

=cut

sub addGeneFeature {
    my $self     = shift ;
    my $gf       = shift ;
    my $function = shift ;
    my $copy     = shift ;
    my $exclude  = shift ; # complement the outcome of the $function discriminator

    $function = 0 if !defined $function ;
    croak "GFF::addGeneFeature() - hmm... $function == 1... did you". 
          "forget your the $function before the $copy argument?\n" if $function == 1 ;

    $copy = 0    if !defined $copy ;
    $exclude = 0 if !defined $exclude ;

    if(!$function || ($exclude xor &$function($gf)) ) {
        if($copy) {
	    push( @{$self->{'feature'}}, $gf->copy($self->version())) ;
	} else {
	    push( @{$self->{'feature'}}, $gf) ;
	}
    }
}

=pod

=over 4

=item nextGeneFeature()

Method to complete remove the next ('head') GeneFeature 
from the GeneFeatureSet, returning it to 
the caller. The order of elements returned is FIFO
relative to addGeneFeature() method calls.

=back

=cut

sub nextGeneFeature {
    my $self = shift ;
    pop( @{$self->{'feature'}} ) ;
}

=pod

=over 4

=item addGFF( $GeneFeatureSet2, \&filter, $copy )

Method to append GeneFeatures in another
GFF::GeneFeatureSet object to the
GFF::GeneFeatureSet object invoking the method
(where "$GeneFeatureSet2" above is the object
reference of the second object, given as the
method argument). Use the "union" method if you
wish to merge two GFF::GeneFeatureSet objects
without duplication of GeneFeature objects. If a
"&filter" (see addGeneFeature above) is not
provided, then the GeneFeature is unconditionally
included. If the $copy argument is set (and
'true' == non-null), then copies (not the
original object references) of the source
GeneFeature objects should be appended to the
invoking object.

=back

=cut

sub addGFF {
    my $self     = shift ;
    my $other    = shift ;

    GFF::GeneFeatureSet->verify($other) ;

    my $function = shift ;
    my $copy     = shift ;

    $function = 0 if !defined $function ;
    croak "GFF::addGFF() - hmm... $function == 1... did you". 
          "forget your the $function before the $copy argument?\n" if $function == 1 ;

    foreach my $gf (@{$other->{'feature'}}){
        $self->addGeneFeature($gf,$function,$copy) ;
    }
}

=pod

=over 4

=item copy( $version )

Method to duplicate the invoking
GFF::GeneFeatureSet object. If the optional
'$version' argument is specified (and greater
than 0) then the new copy is cast into the
specified version.  This allows for GFF version
casting of GeneFeatureSets.

=back

=cut

sub copy {
    my $self    = shift ;
    my $version = shift ;  
    my $no_sequences = shift ;  

    $version = $self->version() if !(defined($version) and $version > 0) ;

    my $gfs = GFF::GeneFeatureSet->new($version,$self->region()) ;

    foreach my $gf (@{$self->{'feature'}}) {
        $gfs->addGeneFeature($gf,0,1) ; # makes a copy
    }
    return $gfs ;
}


=head1 GFF::GeneFeatureSet Input/Output Methods

=over 4

=item read_header( \*INPUT )

Reads in the GFF file header meta-comments from 
the top (head) of an input file, (re)setting the
meta-data of the current object accordingly.

Reading stops at the first non-comment field 
encountered, returning non-null if meta-comments
encountered.

=back

=cut

sub read_header {
    my $self      = shift ;
    my $file      = shift ;

    if( !defined $file ) {
	carp(" GFF::GeneFeatureSet::read_header(): sorry - can't read a GFF with no file!");
	return 0;
    }

    my $n = 0 ;
    seek $file,0,0 ; # reset to start, just-in-case?
    while( defined( my $line = <$file> ) ) {
	chomp( $line );
	next unless $line;
        $self->version($1),$n++,next        if $line =~ /^##gff-version\s+(\d+)/ ;      
        $self->date($1,$2,$3),$n++,next     if $line =~ /^##date\s+(\d{4})\-(\d{1,2})\-(\d{1,2})/ ; 
        $self->region($1,$3,$5),$n++,next   if $line =~ /^##sequence-region\s+(\S+)(\s+([\+-]?\d+))?(\s+([\+-]?\d+))?/ ;
        $self->source_version($1),$n++,next if $line =~ /^##source-version\s+(.+)$/ ;
	$self->comment($1),next             if $line =~ /^#\s*(.+)$/ ;       

	last ; # ignore remainder of file after initial header...
    }
    return $n ;     
}

=pod

=over 4

=item read( \*INPUT, \&filter, \&converter)

Create GeneFeature object for each line of a
stream from a GFF formatted file and add
references to them to the GeneFeatureSet object.
An optional "&converter" function may be used to
modify or filter input lines on the fly. This
function should take a $string and return a
$string; lines converting to an empty string are
skipped by the read. An optional, user-defined
predicate (boolean) function "&filter" tests the
resulting GeneFeature object (given as its
argument) for conditional inclusion in the
GeneFeatureSet object, based upon user criteria.
If a "&filter" is not provided, then the
GeneFeature is unconditionally included. Comment
lines (lines beginning with a "#") are also
skipped. Use the GFF::trace(1) command to 
have read input tracing to \*STDERR.

=back

=cut

sub read {
    my $self      = shift ;
    my $file      = shift ;
    my $function  = shift ;
    my $function2 = shift ;

    my ($gf,$count,$line);

    GFF::TRACE(1,"Entering GFF::GeneFeatureSet::read()\n") ;

    $count = 0;

    if( !defined $file ) {
	carp(" GFF::GeneFeatureSet::read(): sorry - can't read a GFF with no file!");
	return 0;
    }

    # loop over file, make new gf objects and add them

    seek $file,0,0 ; # reset to start, just-in-case?
    while( defined( $line = <$file> ) ) {
	$count++ ;
        # diagnostic user feedback during reading of really large files...
        GFF::TRACE(1,'.')  if !($count % 10) ;
	GFF::TRACE(1,"\n") if !($count %600) ;
        
	chomp( $line );
	next unless $line;
#
#       check for GFF version metalines
#
#       parse for meta-comment fields; reset if found
#
        $self->version($1),next        if $line =~ /^##gff-version\s+(\d+)/ ;      
        $self->date($1,$2,$3),next     if $line =~ /^##date\s+(\d{4})\-(\d{1,2})\-(\d{1,2})/ ; 
        $self->region($1,$3,$5),next   if $line =~ /^##sequence-region\s+(\S+)(\s+([\+-]?\d+))?(\s+([\+-]?\d+))?/ ;       
        $self->source_version($1),next if $line =~ /^##source-version\s+(.+)$/ ;
	$self->comment($1),next        if $line =~ /^#\s*(.*)$/ ;       

# convert file format using function
# can also cause line to be skipped if returns ''

	if($function2){
	    $line = &$function2( $line );
	    next unless $line;
	}

	$gf = GFF::GeneFeature->new_from_line( $line, $self->version() );

	$self->addGeneFeature( $gf, $function );
	
	$count++;
    }
    GFF::TRACE(1,"\nExiting GFF::GeneFeatureSet::read()\n\n") ;

    return $count;
}

=pod

=over 4

=item pipe( \*INPUT, \*OUTPUT, \&filter )

Pipe a GFF file from \*INPUT to \*OUTPUT.
An optional "&filter" function may be used to
modify or filter input features on the fly. This
function should take a $gf and return a
$gf or 0; features returning 0 are skipped.
Comment lines (lines beginning with a "#") are 
piped verbatim to output.

The merit of this method is that it does not
read the whole GFF file into memory, so one
can use a filter function to make small, simple
sequential modifications to a GFF file without
incurring a large memory overhead.

Use the GFF::trace(1) command to 
have read input tracing to \*STDERR.

=back

=cut

sub pipe {
    my $self      = shift ;
    my $file1     = shift ;
    my $file2     = shift ;
    my $filter  = shift ;

    my ($gf,$count,$line);

    GFF::TRACE(1,"Entering GFF::GeneFeatureSet::pipe()\n") ;

    $count = 0;

    if( !(defined $file1 && defined($file2)) ) {
	carp(" GFF::GeneFeatureSet::pipe(): sorry - can't pipe GFF without both a source and target file!");
	return 0;
    }

    # loop over file, make new gf objects and add them

    seek $file1,0,0 ; # reset $file1 to start, just-in-case? but just append to $file2...
    while( defined( $line = <$file1> ) ) {

        # diagnostic user feedback during reading of really large files...
        GFF::TRACE(1,'.')  if !($count++%10) ;
	GFF::TRACE(1,"\n") if !($count%600) ;
        
	chomp( $line );
	next unless $line;
#
#       check for GFF version metalines
#
#       parse for meta-comment fields; reset if found
#
	if($line =~ /^#/) {
	   print $file2 "$line\n" ; 
	   next ;
        }

	$gf = GFF::GeneFeature->new_from_line( $line, $self->version() );

	next if($filter && !($gf = &$filter($gf)))  ;

	$gf->dump( $file2 );
	
	$count++;
    }
    return $count ;
}

=pod

=over 4

=item read_msp( \*INPUT, $source, \&filter, \&name_parser )

Create HomolGeneFeature object for each line of a
stream from an MSP file (output from MSPcrunch)
and add references to them to the invoking
GeneFeatureSet object. User specifies the origin
of the msp as the $source string argument to the
function. The &name_parser argument is user
defined function which takes two arguments $group
and \@array, where @array is assumed to be the
(delimiter split) remainder of a line of data
beyond the group field. An optional, user-defined
predicate (boolean) function "&filter" tests the
resulting GeneFeature object (given as its
argument) for inclusion in the GeneFeatureSet
object, based upon user criteria. If a "&filter"
is not provided, then the GeneFeature is
unconditionally included. Comment lines (lines
beginning with a "#") are also skipped.

=back

=cut

sub read_msp {
    my $self      = shift ;
    my $file      = shift ;
    my $source    = shift ;
    my $function  = shift ;
    my $function2 = shift ;

    if( !defined $file ) {
	carp("sorry - can't read a GFF with no file!");
	return 0;
    }

    # loop over file, make new gf objects and add them
    my $count = 0;
    while(defined(my $line = <$file>)) {
	chomp($line);

	# parse comments
	$self->comment($1) if $line =~ /^#\s*(.+)$/ ; 
      
	# skip lines not starting with a number (e.g. "Gapped" from est_genome)
	next unless $line =~ /^\s*\d+/;
	my $gf = GFF::HomolGeneFeature->new_from_msp( $line, $source, $function2, $self->version() );
	$self->addGeneFeature($gf,$function);
	$count++;
    }
    return $count;
}

=pod

=over 4

=item read_parse( \*INPUT, \&parser, \&filter, $type )

Method uses the
GFF::GeneFeature::new_from_parse() or
HomolGeneFeature::new_from_parse() method to
create GeneFeature objects from each line of the
specified stream input filehandle, using a
user-defined parsing function "&parser" (see
new_from_parse for &parser protocol details). The
resulting GeneFeature object references are added
to the invoking GeneFeatureSet object. An
optional, user-defined predicate (boolean)
function "&filter" test the resulting GeneFeature
objects (given as its first argument) for
inclusion in the GeneFeatureSet object, based
upon user criteria. If a "&filter" is not
provided, then the GeneFeature is unconditionally
included. Comment lines (lines beginning with a
"#") are also skipped. The optional $type
parameter may be set to "HomolGeneFeature" or
"GeneFeature" (defaults to "GeneFeature", if not
specified).

=back

=cut

sub read_parse {
    my $self      = shift ;
    my $file      = shift ;
    my $function  = shift ;
    my $function2 = shift ;
    my $class      = shift ;

    my ($gf,$count,$line);

    if( !defined $file ) {
	carp("sorry - can't read a GFF with no file!");
	return 0;
    }

    # loop over file, make new gf objects and add them
    $count = 0;
    while( defined($line = <$file>) ) {
        chomp( $line );
#
#       parse for meta-comment fields; reset if found
#
        $self->version($1),next        if $line =~ /^##gff-version\s+(\d+)/ ;      
        $self->date($1,$2,$3),next     if $line =~ /^##date\s+(\d{4})\-(\d{1,2})\-(\d{1,2})/ ; 
        $self->region($1,$3,$5),next   if $line =~ /^##sequence-region\s+(\S+)(\s+([\+-]?\d+))?(\s+([\+-]?\d+))?/ ;       
        $self->source_version($1),next if $line =~ /^##source-version\s+(.+)$/ ;
	$self->comment($1),next        if $line =~ /^#\s*(.+)$/ ;  # blindly parse any other comment/metacomment

	if(defined $class && $class eq 'GFF::HomolGeneFeature'){
	    $gf = GFF::HomolGeneFeature->new_from_parse( $line, $function, $self->version() );
	}else{
	    $gf = GFF::GeneFeature->new_from_parse( $line, $function, $self->version() );
	}
	if($gf){
	    $self->addGeneFeature($gf,$function2);
	    $count++;
	}
    }
    return $count;
}

=pod

=over 4

=item eachGeneFeature()

Method to return an array of refs to the
GeneFeatures in a GeneFeatureSet object.

=back

=cut

sub eachGeneFeature {
    my $self = shift;
    my @array;
    my $gf;

    foreach $gf ( @{$self->{'feature'}} ){
	push(@array,$gf);
    }
    return @array;
}

=pod

=over 4

=item dump_header( \*OUTPUT, $no_comments )

Method to dump meta-comment fields that this
GeneFeatureSet object knows about, currently, the
##version, ##date, ##sequence-region and ##source-version
(if defined). Returns non-null if any such fields are
defined (and thus printed).

By default, GFF::GeneFeatureSet recorded comments
(from input or recorded by the 'comments' method)
are all dumped after the above meta-comments
are dumped, consolidated as a block irrespective
of their location in the original input file.
If the optional $no_comments flag is defined and non-null,
then the output of such comments is suppressed 

=back

=cut

sub dump_header {
    my $self = shift ;
    my $file = shift ;
    my $no_comments = shift ;

    my $n = 0 ;

    if( !defined $file ){
	$file = \*STDOUT;
    }

    print $file "##gff-version ",$self->version(),"\n"
        if defined $self->version() and ++$n ;      
    print $file "##date ",join('-',$self->date()),"\n" 
        if defined $self->date() and ++$n ;
    print $file "##sequence-region ",join(' ',$self->region()),"\n" 
        if defined $self->region() and ++$n ;
    print $file "##source-version ",$self->source_version(),"\n" 
        if defined $self->source_version() and ++$n ;

    unless(defined($no_comments) and $no_comments) {
	foreach my $comment (@{$self->{'comments'}}) {
	    print "# $comment\n" ; ++$n ;
	}
    }

    return $n ; # non-null if any fields printed      
}

=pod

=over 4

=item dump( \*OUTPUT, $tab, $newline, $flen, $inorder, $tag, $no_sequences )

Dump out a GeneFeatureSet object (via
GFF::GeneFeature::dump_string()). If \*OUTPUT is
not given, \*STDOUT is used. The method returns
'true' ("1") if non-empty GeneFeatureSet object;
'false' ("0") if the GeneFeatureSet object is
empty. 

The "$tab" argument is a boolean flag,
where a "true" (non-null) value directs the use
tab as the field delimiter in the output line.
Otherwise, blank space is used as the delimiter
(default is "true" if not specified). 

The "$newline" argument is passed to GeneFeature
dump_string, which passes it on to
GFF::GeneFeature::dump_group() to affect group
printing. 

The "$flen" argument is a boolean flag,
where a non-null value stipulates that the length
of the current output line should be printed as
an extra field at the end of the output line
(assumed null if not specified. Note: the extra
length of this field is *not* added to the
displayed line size, but the extra field is tab
delimited, if $tab is set). 

The "$inorder" argument is a boolean flag (default:
true) which forces sorting of the GeneFeatures
during dump by "start" coordinate order. Users
may wish to suppress sorting (i.e. explicitly set
$inorder to 0 "false"), for performance reasons,
when the GeneFeatureSet file is large.

The optional $tag argument controls dumping of
[group] fields (see GFF:GeneFeature::dump()).

The method returns the number of features dumped.

=back

=cut

sub dump {
    my $self    = shift ;
    my $file    = shift ;
    my $tab     = shift ;
    my $newline = shift ;
    my $flen    = shift ;
    my $inorder = shift ;
    my $tag     = shift ;

    my $gf;
    my $n=0 ;

    $inorder = 1 if !defined $inorder ; # default to "true" if omitted
    if( $inorder ) {
        foreach $gf ( $self->order_gf() ) {
	    $n++ ;
    	    $gf->dump($file,$tab,$newline,$flen,$tag);
        }
    } else {
        foreach $gf ( @{$self->{'feature'}} ) {
    	    $n++ ;
    	    $gf->dump($file,$tab,$newline,$flen,$tag);
	}
    }
    return $n ;
}

=pod

=over 4

=item dump_matches( \*OUTPUT, $tab, $show_nomatches )

Dump out a GeneFeatureSet object (via method to 
dump out a GeneFeature object) along with
information about (overlap) matching GeneFeatures. 
If \*OUTPUT is not given, \*STDOUT is used. The
"$tab" argument is a boolean flag, where a "true" 
(non-null) value directs the use tab as the field
delimiter in the output line (assumed "true" 
(non-null) if not specified). Otherwise, blank space is
used as the delimiter.

Normally, only features with matches are dumped. The
optional boolean flag '$show_nomatches' when defined and non-null,
directs that 'no match' records are reported too.

=back

=cut

sub dump_matches {
    my $self           = shift;
    my $file           = shift;
    my $tab            = shift;
    my $show_nomatches = shift ;
    my $gf;

    foreach $gf ( @{$self->{'feature'}} ) {
	next unless ( defined($show_nomatches) or
		      $gf->getMatches_logical() ) ; 
	$gf->dump_matches($file,$tab);
    }
    
}

=pod

=head1 GFF::GeneFeatureSet Access Methods

The various GFF::GeneFeatureSet object parameters 
may be set or queried by the following access methods.
All the methods can take arguments as noted to 
set the variable. With or without an argument, the methods
return the current (or newly set) values, as a 
list, except as specifically noted below:

=cut

=pod

=over 4

=item comment( $comment )

Adds a comment line to the GFF::GeneFeatureSet;
These are conditionally echoed to output by dump_header.

=back

=cut

sub comment {
    my $self  = shift ;
    my $comment = shift ;
    return if !defined($comment) ;
    push( @{$self->{'comments'}},$comment) ;
}

=pod

=over 4

=item version( $version )

object protocol version (see
GFF::GeneFeatureSetObject::version()).

=back

=cut

=pod

=over 4

=item date( $year, $month, $date )

date of the GeneFeatureSet file (meta-comment
##date line).

=back

=cut

sub date {
    my $self  = shift ;
    my $year  = shift ;
    my $month = shift ;
    my $day   = shift ;

    # limited validation of dates
    if( defined $year and $year and $year =~ /\d{4}/ ) { 
        $self->{'date'}->[0] = $year ;
    }
    if( defined $month and $month and
           $month =~ /(\d{1,2})/ and $1 >0 and $1<= 12 ) {
        $self->{'date'}->[1] = $month ;
    }
    if( defined $day and $day and 
           $day =~ /(\d{1,2})/ and $1 >0 and $1<= 31 ) {
        $self->{'date'}->[2] = $1 ; 
    }
    return @{$self->{'date'}} ;
}

=pod

=over 4

=item region( $sequence, $start, $end )

sequence region of the GeneFeatureSet file
(##sequence-region meta-comment line).

=back

=cut

sub region {
    my $self    = shift ;
    my $seqname = shift ; 
    my $start   = shift ;
    my $end     = shift ;

    # limited validation of data
    if( defined $seqname and $seqname and $seqname =~ /\w+/ ) { 
        $self->{'region'}->[0] = $seqname ;
    }
    if( defined $start and $start and $start =~ /^[\+-]?\d+$/ ) { 
        $self->{'region'}->[1] = $start ;
    }
    if( defined $end and $end and $end =~ /^[\+-]?\d+$/ ) { 
        $self->{'region'}->[2] = $end ;
    }
    return @{$self->{'region'}} ;
}

=pod

=over 4

=item source_version( $description )

##source_version meta-comment line.

=back

=cut

sub source_version {
    my $self = shift ;
    my $desc = shift ;

    if( defined $desc ) { 
        $self->{'source_version'} = $desc ;
    }
    return $self->{'source_version'} ;
}

=pod

=head1 GFF::GeneFeatureSet Simple Set Operations

Note - all set comparisons are by reference only - 
comparisons cannot check if two GeneFeature objects
with distinct reference pointers actually contain the same data!

=over 4

=item member( $GeneFeature )

Method to test if a GeneFeature object (object
reference "$GeneFeature") is a member of an
existing GeneFeatureSet object

=back

=cut

=pod

=over 4

=item union( $GeneFeatureSet2 )

Method to generate a new GeneFeatureSet object
which is a union of 2 GeneFeatureSet objects (the
one invoking the method plus the one specified by
the GeneFeatureSet object reference argument,
$GeneFeatureSet2, to the method).

=back

=cut

sub union {
    my $self  = shift;
    my $other = shift;

    GFF::TRACE(1,"Entering GFF::GeneFeatureSet::union()\n") ;

    my $new = $self->new($self->version(),$self->region()) ;
    my $gf;

# first add everything in 1st object
    GFF::TRACE(1,"Adding everything in the first set:\n") ;
    my $n = 0 ;
    foreach $gf (@{$self->{'feature'}}){
        GFF::TRACE(1,'.')  if !($n++%10) ;
	GFF::TRACE(1,"\n") if !($n%600) ;
	$new->addGeneFeature($gf);
    }

# second add everything in 2nd object but not in 1st object
    GFF::TRACE(1,"\nAdding additional features from the second set:\n") ;
    $n = 0 ;
    foreach $gf (@{$other->{'feature'}}){
        GFF::TRACE(1,'.')  if !($n++%10) ;
	GFF::TRACE(1,"\n") if !($n%600) ;
	unless($self->member($gf)){
	    $new->addGeneFeature($gf);
	}
    }
    GFF::TRACE(1,"Leaving GFF::GeneFeatureSet::union()\n") ;

    return $new;
}

=pod

=over 4

=item intersection( $GeneFeatureSet2 )

Method to generate a new GeneFeatureSet object
which is every GeneFeature in first (invoking)
GeneFeatureSet GeneFeature set that is also a
member of the second GeneFeatureSet GeneFeature
set (argument $GeneFeatureSet2 above).

=back

=cut

sub intersection {
    my $self = shift;
    my $other = shift;

    GFF::TRACE(1,"Entering GFF::GeneFeatureSet::intersection()\n") ;

    my $new = $self->new($self->version(),$self->region()) ;
    my $gf;

# add everything in 1st object only it is in the 2nd object
    my $n = 0 ;
    foreach $gf (@{$self->{'feature'}}){
        GFF::TRACE(1,'.')  if !($n++%10) ;
	GFF::TRACE(1,"\n") if !($n%600) ;
	if($other->member($gf)){
	    $new->addGeneFeature($gf);
	}
    }
    GFF::TRACE(1,"Leaving GFF::GeneFeatureSet::intersection()\n") ;

    return $new;
}

=pod

=over 4

=item difference( $GeneFeatureSet2 )

Method to generate a new GeneFeatureSet object
which is everything in first (invoking)
GeneFeatureSet GeneFeature set that is not in the
second GeneFeatureSet GeneFeature set (argument
$GeneFeatureSet2 above).

=back

=cut

sub difference {
    my $self = shift;
    my $other = shift;

    GFF::TRACE(1,"Entering GFF::GeneFeatureSet::difference()\n") ;

    my $new = $self->new($self->version(),$self->region()) ;
    my $gf;

# add everything in 1st object unless it is in the 2nd object
    my $n = 0 ;
    foreach $gf (@{$self->{'feature'}}){
        GFF::TRACE(1,'.')  if !($n++%10) ;
	GFF::TRACE(1,"\n") if !($n%600) ;
	unless($other->member($gf)){
	    $new->addGeneFeature($gf);
	}
    }
    GFF::TRACE(1,"Leaving GFF::GeneFeatureSet::difference()\n") ;

    return $new;
}

=pod

=head1 GFF::GeneFeatureSet Feature Partition Methods

This series of methods partition a GeneFeatureSet 
into subsets based upon specified attribute or feature
criteria.

=over 4

=item filter( \&function )

Method to generate a new GeneFeatureSet object
based on a filtered version of the object.
Filtering is carried out by passing a reference
to a subroutine, "&function", which is applied to
each GeneFeature object in the invoking
GeneFeatureSet object. This user-defined
&function should be designed to accept a
reference to a single GeneFeature object and to
return the predicate (boolean) outcome of some
test upon that object reflecting the user filter
criterion: 1 implying inclusion, 0 implying
exclusion of the GeneFeature from the new
GeneFeatureSet object set.

=back

=cut

sub filter {
    my $self     = shift;
    my $function = shift;
    my $gf;

    my $newSet = $self->new($self->version(),$self->region());
    
    foreach $gf (@{$self->{'feature'}}){
	$newSet->addGeneFeature($gf,$function);
    }
    return $newSet;
}

=pod

=over 4

=item exclude( \&function )

Method to generate a new GeneFeatureSet object
based on a filtered version of the object. This
method is simply the negation of filter(), in that
the discriminant function values 1 implies *exclusion*
and 0 implied inclusion. This may be handy
in that the same discriminant functions can
therefore be used to partition a set by
using filter() and exclude() sequentially.

=back

=cut

sub exclude {
    my $self     = shift;
    my $function = shift;
    my $gf;

    my $newSet = $self->new($self->version(),$self->region());
    
    foreach $gf (@{$self->{'feature'}}){
	$newSet->addGeneFeature($gf,$function,0,1);
    }
    return $newSet;
}

=pod

=over 4

=item annotate( \&function )

Method to apply a user defined annotation function 
to all every feature records in a GFF::GeneFeatureSet.
To do this, the reference pointer to each
GFF::GeneFeature is passed to the function for 
manipulation. 

=back

=cut

sub annotate {
    my $self     = shift;
    my $function = shift;

    croak "GFF::GeneFeatureSet::annotate(): needs a function pointer argument!\n"
	if !(defined($function) and ref($function) =~ /CODE/) ;

    foreach my $gf (@{$self->{'feature'}}){
	&$function($gf) ;
    }
}

=pod

=over 4

=item rewriteField( $field, $target, $rewrite, $returnAll, $record )

This method creates a new GeneFeatureSet object
containing GeneFeatures in which a designated
$field (specified by a string value 'SEQNAME,'
'SOURCE', 'FEATURE' or 'GROUP' - case
insensitive) in all GeneFeatures that match the
specified $target value (which can either be a
simple identifier or a Perl regular expression.

GFF Version :, 'GROUP' fields $target matches
the field itself.

GFF Version 2: 'GROUP' fields, $target should be a 
simple identifier matching a tag of some tag-value). 

If $field is 'SEQNAME', the return object 
'sequence-region' is renamed to the first 
field matched.

The specified field is overwritten with the
specified $rewrite value (which may be a simple
identifier or a full Perl search & replace
expression, namely, "s/<search>/<replacement>/"
). Note: the 's' should be the very first
character in the string and should be immediately
followed by a non-alphanumeric delimiter
character, for this to work properly).
Backreferenced ($1 et al.) replacement values and
the 'g','i', & 'x'  search modifiers are
permitted.

Note: if the 'target' is simply an asterix '*'
('wildcard'), then all fields of the designated
type are rewritten with $rewrite specification.

The new GeneFeatureSet object returned only
contains GeneFeatures which triggered a rewrite,
unless the optional 'returnAll' boolean flag
argument is defined and non-null, in which case a
COPY of all GeneFeatures is returned, whether or
not it was modified.

If the optional $record argument is set (to a
simple tag identifier), then the old $field value 
is recorded as the $record [group] field tag value.

=back

=cut

sub rewriteField {
    my $self        = shift ;
    my $field       = shift ;   # 'SEQNAME' 'SOURCE' or 'FEATURE' or 'GROUP'
    my $target      = shift ;   # can be a tag or Perl regex
    my $rewrite     = shift ;   # new field value to be used; may be Perl 'search and replace' expression 
    my $returnAll   = shift ;   # if defined and non-null, then all GeneFeatures are returned (modified or not)
    my $record      = shift ;   # save old field value under this [attribute] tag

    $returnAll = 0 if !defined $returnAll ;

    croak "GeneFeatureSet::rewriteField() needs a field specification\n"   if !(defined $field  and $field) ;
    croak "GeneFeatureSet::rewriteField() needs a target\n"                if !(defined $target and $target) ;
    croak "GeneFeatureSet::rewriteField() needs a new tag specification\n" if !(defined $rewrite and $rewrite) ;

    my $wildcard = ($target eq '*')? 1 : 0 ;

    my $ftype ; 
    if( $field =~ /^SEQNAME$/i )  {
	$ftype= 1 ;
    } elsif ( $field =~ /^SOURCE$/i ) {
	$ftype = 2 ;
    } elsif ( $field  =~ /^FEATURE$/i ) {
	$ftype = 3 ;
    } elsif ( $field  =~ /^GROUP$/i ) {
	$ftype = 4 ;
    } else {
       croak "GeneFeatureSet::rewriteField() field specification not recognized\n"  ;
    }
    my ($search, $replace, $flags) ;
    if( $rewrite =~ /^s(\W)([^\1]+)\1([^\1]*)\1([gix]+)/ ) { # search & substitute rewrite?
       $search  = $2 ;
       $replace = $3 ;
       $flags   = $4 ;
    }
    my @region = $self->region() ;
    my $newSet = $self->new($self->version(),@region) ;
    my $newtag = $region[0] if $ftype == 1 ; # in case $self is empty?

    foreach my $ogf (@{$self->{'feature'}}){
        my $gf = $ogf->copy() ;
        my ($oldtag, $match) ;
        my $save = 0 ;
	if( $ftype == 1 )  {
	    $oldtag = $gf->seqname() ;
            $match = $wildcard || eval { $oldtag =~ /^$target$/i } ;
        } elsif ( $ftype == 2 ) {
	    $oldtag = $gf->source() ;
            $match = $wildcard || eval { $oldtag =~ /^$target$/i } ;
        } elsif ( $ftype == 3 ) {
	    $oldtag = $gf->feature();
            $match = $wildcard || eval { $oldtag =~ /^$target$/i } ;
        } elsif ( $ftype  == 4 ) {
	    if($gf->version() == 1) {
                $oldtag = $gf->group() ;
	    } else {
	        $oldtag = $gf->tagValue($target) ;
	    }
            $match = 1 if defined $oldtag ;
	}
        if($@) { 
            croak "GeneFeatureSet::rewriteField() error: invalid target match pattern!\n" ;
        } elsif ($match) {
            if(defined($record)) {
		$gf->tagValue($record,0,$oldtag) ;
	    }
            #
            # Triggered: overwrite the oldtag with the new
            #
            $save = 1 ;
            if( defined($search) ) { # search & substitute rewrite?
                $newtag = $oldtag ; # to be modified...
                eval "\$newtag =~ s/$search/$replace/$flags" ; 
                if($@) { 
	            croak "GFF::GeneFeatureSet::rewriteField():  improper \$rewrite expression: '$rewrite'!" ;
                }                
	    } else {
		$newtag = $rewrite ; # simple identifier substitution
	    }
	    if( $ftype == 1 )  {
		$gf->seqname($newtag) ;
            } elsif (  $ftype == 2 ) {
		$gf->source($newtag) ;
            } elsif (  $ftype == 4 ) {
	        if($gf->version() == 1) {
                    $gf->group($newtag) ;
		} else {
		    $gf->tagValue($target,0,$newtag) ;
		}
            }  else { # a FEATURE...anything else was caught above
		$gf->feature($newtag);
	    }
	} # else do nothing...

	$newSet->addGeneFeature($gf) if $returnAll or $save ;
    }
    if($field =~ /^SEQNAME/i) {
        # I blissfully assume that the last SEQNAME seen is the new name...
        $newSet->region($newtag,$region[1],$region[2]) ; 
    } else {
        $newSet->region(@region) ; 
    }
    return $newSet ;
}

=pod

=over 4

=item cluster( \&comparison, $single )

Method to build an array of GeneFeatureSet
objects each containing a group of GeneFeature
objects sharing some shared attribute [pairwise
comparison]. The &comparison function reference
is the operational user definition of this shared
attribute, which when given the references to
each of two GeneFeature objects, returns a 1
or a cluster specific string (cluster name)
implying inclusion in a cluster, or 0 implying
exclusion from a cluster, based upon shared
attributes. If the $single flag is set to 1
(assumed 0 if omitted) then all singular
GeneFeature objects not assigned to a group are
added to the array as independent,
single-membered GeneFeatureSet object groups.

=back

=cut

sub cluster {
    my $self = shift;
    my $function = shift;
    my $single = shift;

    my (@gf);
    my ($i,$j,$gf,$gf1,$gf2,$ngf,
	$nc1,$nc2,$nc,%cluster,$ncluster,
	$cluster,@cluster,%cname,$cname);

    GFF::TRACE(1,"Entering GFF::GeneFeatureSet::cluster()\n") ;
#
#  Validate parameters
#
    if( !$function ) {
	return undef ; # undef better for error checking?
    }

# setup array of GFF::GeneFeature objects
    @gf=@{$self->{'feature'}};
    $ngf=scalar(@gf);
    $ncluster=0;

# loop over each
    for($i=0;$i<$ngf;$i++){
	$gf1=$gf[$i];

	GFF::TRACE(1,"\n") if !($i%600) ;
        GFF::TRACE(1,'.')  if !($i%10) ;

# loop over remainder
	for($j=$i+1;$j<$ngf;$j++){
	    $gf2=$gf[$j];

# test relationship between $gf[$i] and $gf[j]
	    if($cname=&$function($gf1,$gf2)){

# if both in existing cluster, merge clusters
		$nc1=$cluster{$gf1};
		$nc2=$cluster{$gf2};
		if($nc1 && $nc2){
# merge clusters if identical
		    if($nc1!=$nc2){
			push(@$nc1,@$nc2);
			foreach $cluster (keys %cluster){
			    if($cluster{$cluster}==$nc2){
			        $cluster{$cluster}=$nc1;
			    }
			}
			undef @$nc2;
			$ncluster--;
		    }
# nc1 exists -> add nc2
		}elsif($nc1){
		    push(@$nc1,$gf2);
		    $cluster{$gf2}=$nc1;
		}elsif($nc2){
		    push(@$nc2,$gf1);
		    $cluster{$gf1}=$nc2;
		}else{
		    $nc=[];
		    push(@cluster,$nc);
		    push(@$nc,$gf1);
		    push(@$nc,$gf2);
		    $cluster{$gf1}=$nc;
		    $cluster{$gf2}=$nc;
		    $ncluster++;
		    $cname{$nc} = $cname ; # unique cluster name?
		}
	    } # endif &$function
	}
    }

# make clusters out of each singleton
    if($single){
	foreach $gf (@gf){
	    unless($cluster{$gf}){
		$nc=[];
		push(@cluster,$nc);
		push(@$nc,$gf);
	    }
	}
    }

# build an array of new objects from set of clusters
    undef @gf;
    my $n = 0 ;
    foreach $cluster (@cluster){
	my $name ; ++$n ;
	if(exists($cname{$cluster}) and 
	   $cname{$cluster}!~ /^1$/) { # not just a boolean 1
	    $name = $cname{$cluster} ; 
	} else {
	    $name = "Cluster # $n" ; 
	}
	my $newSet = $self->new($self->version(),$name) ;
	foreach $gf (@$cluster){
	    $newSet->addGeneFeature($gf);
	}
	push(@gf,$newSet);
    }
    GFF::TRACE(1,"Exiting GFF::GeneFeatureSet::cluster()\n") ;
    @gf;
}


=pod

=over 4

=item features( \&discriminator )

Method to return a hash, key indexed by distinct
feature types and their occurrence. The keys of
the hash are based upon a user-defined
"&discriminator" function which takes a
GeneFeature object reference as its input
parameter and returns a (string) key value label
characteristic of a feature type of interest. If
the &discriminator function returns undef or a
null string for a given GeneFeature object
argument, then that GeneFeature is ignored by the
method. For GeneFeatures with a non-null
&discriminator return value, the method uses this
return value as a hash key to maintain a
cumulative count of the occurrence of this return
value ("feature type") in the current
GeneFeatureSet object.

=back

=cut

sub features {
    my $self = shift ;
    my $function = shift ;

    croak "GFF::GeneFeatureSet::features() requires a \$function argument!\n" if !(defined $function and $function) ;

    my (%hash,$gf,$tmp) ;
    foreach $gf (@{$self->{'feature'}}){
	my $tmp = &$function($gf) ;
        next if !defined $tmp or !$tmp ;

	if( !defined $hash{$tmp} ) { 
           $hash{$tmp} = 1 ;
       } else {
           $hash{$tmp}++ ;
       }
    }
    %hash ;
}

=pod

=over 4

=item theFeature( \&discriminator )

Method to test if the &discriminator function
returns a single value, based upon a user-defined
"&discriminator" function which takes a
GeneFeature object as its input parameter and
returns a (string) key value characterizing the
feature types of interest. If such a singular
feature is found, then it is returned with its
frequency of occurrence in the GeneFeatureSet
file.

=back

=cut

sub theFeature {
    my $self = shift;
    my $function = shift;
    my ($feature,%hash,$n,$featureName);
    $n = 0 ;
    %hash = $self->features($function);
    foreach $feature (keys %hash){
	$featureName = $feature;
	$n++;
    }
    if($n!=1){
	croak "Error in GFF::GeneFeatureSet::theFeature - discriminator function returns multiple ($n) keys";
    }
    $featureName,$hash{$featureName};
}

=pod

=over 4

=item the( \&discriminator )

Non-fatal method to test if the discriminant function 
returns a singular value when acting upon the given
GFF::GeneFeatureSet. Returns the single value if unique;
returns undef otherwise. This is good for obtaining data 
such as the common strand of a GeneFeatureSet 'gene' object.

=back

=cut


sub the {
    my $self     = shift;
    my $function = shift ; # a reference to a function?

    my ($gf,$value,%hash,$theLast);
#
#   Loop over each GFF::GeneFeature
#
    my $n = 0 ;
    foreach $gf (@{$self->{'feature'}}){
       my $value = &{$function}($gf) ;
       if( !$hash{$value} ) {
	   $hash{$value} = 1 ;
	   $theLast = $value ; # saves feature last seen
	   $n++ ;
       } else {
	   $hash{$value}++ ;
       }
    }
    if($n == 1){
	return $theLast, $hash{$theLast};
    } else {
	return ; # undefined...
    }
}

=pod

=over 4

=item group( \&discriminator, $copy )

Method to build an hash of GeneFeatureSet objects
each containing a group of GeneFeature objects
sharing some attribute [fixed, named],
operationally defined by a user &discriminator
function taking a GeneFeature object reference as
the input argument and returning a unique name
string labelling the attribute. The group method
returns a hash of GeneFeatureSet object
references, key indexed by the attribute name
strings. The (meta-comment) 'sequence-region'
start and end coordinates are set to the minimum
and maximum start and end respectively of the
GeneFeatures included into each group
GeneFeatureSet object.

If the (optional) $copy switch is 'true'
(non-null) then the new group GeneFeatureSet
objects (dereferenced by the hash) are composed
of copies of the original GeneFeatures. In other
words, modifications of these new GeneFeatures
will not modify the GeneFeatures objects in the
original object.

=back

=cut

sub group {
    my $self     = shift ;
    my $function = shift ;
    my $copy     = shift ;

    confess "GFF::GeneFeatureSet::group(): needs a defined $function argument!\n" if !defined $function ;

    $copy = 0 if !defined($copy) ;

    my ($gf,$name,%group);
#
#   Construct the clusters of GFF::GeneFeatures...
#
    foreach $gf (@{$self->{'feature'}}){
	if($name = &$function($gf)){
	    push(@{$group{$name}},$gf);
	}
    }
#
# build an array of new objects from set of clusters
#
    foreach $name (keys %group){
	my $groupGFF = $self->new($self->version()) ;
	foreach $gf (@{$group{$name}}){
            if($copy) {
	        $groupGFF->addGeneFeature($gf->copy());
	    } else {
		$groupGFF->addGeneFeature($gf);
            }
	}
#
#       Set region info of each GeneFeatureSet cluster
#
	my ($start,$end) = $groupGFF->min_max_range() ;
	$groupGFF->region($name,$start,$end) ;
#
#       Substitute the new GeneFeatureSet object reference for the array
#
	$group{$name} = $groupGFF ; 
    }
#
#   Return the set of GeneFeatureSet clusters...
#
    return %group ;
}

=pod

=over 4

=item group_value_string( $source, $feature, $tag )

Method to return a string constructed from the list 
of values associated with a given $tag of [group] 
tag-value pairs, from a specified GeneFeature 
record, of the invoking gene feature set, 
which matches the given $source and $feature
(which may be Perl regular expressions). $source and
$feature may also be undef or '*', designating that any
source or feature can match. The $tag should
be a simple tag identifier (*not* a Perl regex).
Only values from the first such $tag encountered in 
the gene feature set are returned. Returns 'undef' 
if no such tag-value list is found.

=back

=cut

sub group_value_string {
    my $self    = shift ;
    my $source  = shift ;
    my $feature = shift ;
    my $tag     = shift ;

    croak "*** Error: GFF::GeneFeatureSet::group_value_string needs a defined \$tag!\n" if !(defined $tag and $tag) ;
    $source  = ".+" if !defined $source or !$source or $source eq '*' ;  # but I let <source> be a wild card
    $feature = ".+" if !defined $feature or !$feature or $feature eq '*' ;

    foreach my $gf (@{$self->{'feature'}}) {
	if( $gf->source()  =~ /^$source$/  and
            $gf->feature() =~ /^$feature$/)
        {
	    my $values = $gf->tagValueList($tag) ;
            return join ' ',@{$values} if $values ;
        } 
    }
    return undef ;
}

=pod

=over 4

=item deleteTag($tag)

Invokes the GFF::GeneFeature::deleteTag() method
on every gene feature object in the current
GeneFeature set. Note: this operation directly 
modifies the feature objects concerned.

=back

=cut

=pod

=over 4

=item label( $membergroup )

Method to add a label ("$membergroup") to each
GeneFeature in this GeneFeatureSet, indexing a
reference pointing back to the invoking
GeneFeatureSet object.

=back

=cut

=pod

=over 4

=item label_pair( $membergroup )

Method to make a list of the GeneFeatureSet
objects with which the GeneFeatures in this
GeneFeatureSet are paired by some label
("$membergroup").

=back

=cut

=pod

=over 4

=item addMember( $ParentGeneFeatureSet,
$membergroup )

Method to add a member record to indicate the
parent GeneFeatureSet object for this particular
grouping ("$membergroup") of GeneFeatures.

=back

=cut

#
# Method to add a member record to indicate the parent GeneFeatureSet object
# for this particular grouping of GFF::GeneFeatures
#

sub addMember {
    my $self = shift;
    my $parent = shift;
    my $membergroup = shift;
    ${$self->{'member'}}{$membergroup} = $parent;
}

=pod

=over 4

=item getMember( $membergroup )

Method to get reference to parent GeneFeatureSet
object of GeneFeature object under this
particular grouping ("$membergroup").

Method to test whether or not GeneFeatureSet
contains members.

=back

=cut

sub getMember {
    my $self = shift;
    my $membergroup = shift;
    my ($gff);
    $gff = $self->{'member'}->{$membergroup};
}

=pod

=over 4

=item containsMembers()

Method to get test whether $membergroups has
members.

=back

=cut

sub containsMembers {
    my $self = shift ;
    return scalar($self->{'member'});
}

=pod

=over 4

=item getAllMembers()

Method to get hash reference of all
$membergroups.

=back

=cut

sub getAllMembers {
    my $self = shift ;
    return $self->{'member'} ; # returns a reference to a hash
}

=pod

=head1 GFF::GeneFeatureSet Geometric Partition
Methods

This series of methods partition a GeneFeatureSet
into subsets based upon coordinate (geometric)
criteria.

=over 4

=item complement($source,$feature,$strand,$tag,$append)

Method returns a set of 'gene features' constructed
from the geometric complement of the calling GeneFeatureSet,
that is, all the sequence intervals *NOT* spanned by the input
feature records. 

The $source, $feature and $strand arguments are used to label the 
<source> and <feature> fields respectively (default: 'GFF_Complement'
for <source> and/or <feature>; '.' for <strand>).

The '$tag' value is a [group] field tag, common to features
in the input GeneFeatureSet, which is used to annotate
the complementary features in the form [Between_$tag $feature1 $feature2]
where $feature1 and $feature2 are the features flanking the 
newly generate complement feature. If complementary features
are generated at the start and/or end of the host coordinate
range, then the special names 'Start' and 'End' are used
as feature names for the $tag labelling.

The '$append' argument, if defined and non-null, directs that
the new GeneFeatureSet of complemented features are appended 
to the invoking GeneFeatureSet, which itself is returned.

=back

=cut

sub complement{
    my $self    = shift ;
    my $source  = shift ;
    my $feature = shift ;
    my $strand  = shift ;
    my $tag     = shift ;
    my $append  = shift ;

    $source  =  'GFF_Complement' if !(defined($source) and $source) ;
    $feature  =  'GFF_Complement' if !(defined($feature) and $feature) ;
    $tag = 0 if !(defined($tag) and $tag) ;

    my ($seqname,$start,$end) = $self->region() ;

    my $gffc = GFF::GeneFeatureSet->new(2,$seqname,$start,$end) ;

    my $prevTV = 'Start' ;
    foreach my $gf ($self->order_gf()) {
	if($gf->start() <= $start) {
	    my $nexts = $gf->end() +1 ;
	    $start= $nexts if $nexts > $start+1 ;
	    $prevTV = $gf->tagValue($tag) if($tag) ;
	    next ;
	}
	my $gfc = new GFF::GeneFeature(2) ;
	$gfc->seqname($seqname) ;
	$gfc->source($source) ;
	$gfc->feature($feature) ;
	$gfc->start($start) ;
	$gfc->end($gf->start()-1) ;
	$gfc->strand($strand) if($strand) ;
	if($tag) {
	    my $newTV = $gf->tagValue($tag) ;
	    $gfc->group("Between_$tag",$prevTV,$newTV) ;
	    $prevTV = $newTV ;
	}
	$start = $gf->end()+1 ;
        $gffc->addGeneFeature($gfc) ;
    }
    # create residual fragment if necessary
    if($start < $end) {
	my $gfc = new GFF::GeneFeature(2) ;
	$gfc->seqname($seqname) ;
	$gfc->source($source) ;
	$gfc->feature($feature) ;
	$gfc->start($start) ;
	$gfc->end($end) ;
	$gfc->strand($strand) if($strand) ;
	$gfc->group($tag,$prevTV,'End') if($tag) ;
        $gffc->addGeneFeature($gfc) ;
    }
    if(defined($append) and $append) {
	$self->addGFF($gffc) ;
        return $self ;
    } else {
	return $gffc ;
    }
}

=pod

=over 4

=item self_overlap(\*OUTPUT, $strict, $exact, $tag)

Method to generate a new GeneFeatureSet object
based on a filtered version of the object based
on any pairwise overlap detected by the
GFF::GeneFeature::overlap_logical() method.
Option to report overlaps to filehandle reference
"\*OUTPUT" (output sent to STDOUT if $file is
omitted or undef; output is suppress if a null
(zero) value is passed to the method for
\*OUTPUT. Note: the method does not test whether
the features being merged match in kind. Thus the
user is responsible for making sure that the
GeneFeatureSet object only contains features
mergeable in the semantic sense (e.g. all the
features are exons predicted by a single
prediction algorithm).

If the optional $strict flag is set 'true'
(non-null) then only overlapping Gene Features
which match identically with respect to
<seqname>, <source>, <feature> and <strand>
(and [group] $tag -- see below) are deleted 
(i.e. only truly 'duplicate' records deleted). 
$strict defaults to 'false' if omitted.

The optional $exact flag specifies that an exact
match is required for overlaps (defaults to false
if the argument is omitted or undef).

If '$tag' is defined, then the specified [group]
tag-value must also match for $strict matches.

=back

=cut

sub self_overlap {
    my $self   = shift ;
    my $file   = shift ;
    my $strict = shift ;
    my $exact  = shift ;
    my $tag    = shift ;

    my (@gf,@del);
    my ($i,$j,$gf);

    $strict = 0 if !defined $strict ;
    $exact  = 0 if !defined $exact ;
    $tag    = 0 if !defined $tag ;

    GFF::TRACE(1,"Entering GFF::GeneFeatureSet::self_overlap()\n") ;
    @gf=$self->order_gf();
    $i = 0 ;
    foreach $gf (@gf) {
        $del[$i++] = 0 ;
    }
#
#   loop over every GFF::GeneFeature, skipping previously deleted elements
#
    my $gffo = $self->new($self->version(),$self->region()) ;

    Outer: for($i=0;$i<scalar(@gf);$i++){
	GFF::TRACE(1,"\n") if !($i%600) ;
        GFF::TRACE(1,'.')  if !($i%10) ;
        next if $del[$i] ;               # skip if already deleted
	$gffo->addGeneFeature($gf[$i]) ; # add element to the new set
	my $gf1=$gf[$i] ;
#
#       Loop over GFF::GeneFeatures after currently examined one,
#       marking 'overlaps' (by one's constraints) for deletion
#
	for($j=$i+1;$j<scalar(@gf);$j++){
            next if $del[$j] ; # skip if already deleted
	    my $gf2=$gf[$j] ;

	    next Outer if $gf1->end() < $gf2->start() ; # short circuit non-overlap test!

	    my $string ;
	    if( $exact ) {
		$string = $gf1->match_logical  ($gf2,$file) ;
	    } else {
		$string = $gf1->overlap_logical($gf2,$file) ;
            }
            my ($value1, $value2) ;
            if( $string and 
                (!$strict or 
                    ( $gf1->seqname() eq $gf2->seqname() and
                      $gf1->source()  eq $gf2->source()  and
                      $gf1->feature() eq $gf2->feature() and
                      $gf1->strand()  eq $gf2->strand()  and
                      ( !$tag or 
			 ( defined($value1 = $gf1->tagValue($tag)) and
                           defined($value2 = $gf2->tagValue($tag)) and
                           $value1 eq $value2 
			  )
                       )
                   ) 
                )
             ) {
		   print $file "  OVERLAP: ",$gf1->source(),":",$gf1->feature(),
                               " $string ",$gf2->source(),":",$gf2->feature(),
                               " (DELETE)\n" if $file ;

		   $del[$j] = 1 ;  # mark as deleted
	    }
	}
    }    
    GFF::TRACE(1,"\nExiting GFF::GeneFeatureSet::self_overlap()\n\n") ;
    return $gffo;
}

=pod

=over 4

=item annotate_overlaps($tag,$only_one)

Method to annotate the (self) overlap segments of
features in the the invoking GeneFeatureSet object,
based upon any pairwise overlap detected.

The overlaps are annotated as tag-values:

'Overlap_[Left|Right] <Sequence name> <position>'

where <position> is the first base pair position
shared (either on the left or the right)
with the adjacent overlapping feature.

The optional '$tag' argument can be used to
specify an alternate tag-value prefix instead
of 'Overlap'.

Note that the <position> is GFF relative, hence,
feature <strand> insensitive, i.e. Overlap_Left
for a forward stranded feature clips the 5' end 
of a feature, whereas for a reverse stranded
feature, it clips the 3' end... This could be
important to consider if the 'features' are
sequences in a link contig...

This method is typically intended to characterize
a clean, non-redundant minimal tiling path of 
features (e.g. link contig). Therefore, the method 
generally expects to only find one overlap on each
side of a given feature; however this is not strictly
enforced unless the '$single' flag is set non-null, 
in which case the method *only* annotates the first 
overlap it encounters.

The method also normally complains about features which are 
exactly aligned at either the 5' or 3' end, or situations 
in which one feature appears to be embedded in another. 
Such features are *not* annotated as overlaps.
However, if the above situation are not an issue,
then the optional boolean '$ignore' flag may be set,
in which case no complaints are made and every
overlap is annotated irrespective of its character.

=back

=cut

sub annotate_overlaps {
    my $self   = shift ;
    my $tag    = shift ;
    $tag = 'Overlap' if !defined $tag ;
    my $single = shift ;
    my $ignore = shift ;

    GFF::TRACE(1,"Entering GFF::GeneFeatureSet::annotate_overlaps()\n") ;

    my (@gf,$i,$j,$gf);
    @gf=$self->order_gf();
    Outer: for($i=0;$i<scalar(@gf);$i++){
	GFF::TRACE(1,"\n") if !($i%600) ;
        GFF::TRACE(1,'.')  if !($i%10) ;
	my $gf1=$gf[$i] ;
	for($j=$i+1;$j<scalar(@gf);$j++){
	    my $gf2=$gf[$j] ;

	    next Outer if $gf1->end() < $gf2->start() ; # short circuit non-overlap test!

	    # The overlap test?
	    my @result = $gf1->match($gf2,-1) ; # any overlap?
	    if( $result[0] ) {
                # annotate overlap
		my $f5=$result[1];  
		my $f3=$result[2];
		unless($ignore) {
		    unless($f5) {
			carp "GFF::GeneFeatureSet::annotate_overlaps(): 5' ends of features are aligned?\n" ;
			$gf1->dump(\*STDERR) ;
			print STDERR "\n\tand\n" ;
			$gf2->dump(\*STDERR) ;
			next ;
		    } elsif($f3<=0) {
			carp "GFF::GeneFeatureSet::annotate_overlaps(): embedded feature? 3' end of feature 2 is <= 3' end of feature 1?\n" ;
			$gf1->dump(\*STDERR) ;
			print STDERR "\n\tand\n" ;
			$gf2->dump(\*STDERR) ;
			next ;
		    }
		}
		# normal overlap (or annotating all overlaps anyhow...)?
		$gf1->tagValueList("$tag\_right",[$gf2->Sequence(),$gf2->start()],1) ;
		$gf2->tagValueList("$tag\_left",[$gf1->Sequence(),$gf1->end()],1) ;

                # short circuit analysis if only first 
                # (closest) overlapping feature is to be annotated
		next Outer if $single ; 
	    }
	}
    }    
    GFF::TRACE(1,"\nExiting GFF::GeneFeatureSet::annotate_overlaps\n\n") ;
}

=pod

=over 4

=item self_overlap_merge( $tolerance, $strand, $group_tag, $addscores )

Method to generate a new GeneFeatureSet object
based on a merging of overlapping (but otherwise
similar) GeneFeature records using the
GFF::GeneFeature::overlap_merge() method i.e. if
two GeneFeatureSet's, 1-10 and 8-16 then replace
with one new one, 1-16. See caveat note above
(under self_overlap) about feature merger
semantics.

The optional $tolerance value is passed onto the 
GFF::Genefeature overlap_merge() method.

Defined, non-null '$strand' boolean flag forces
the merging to be strand sensitive.

The optional $group_tag argument is passed to the
GFF::GeneFeature::overlap_merge() function (which see).
This may be a simple [group] tag identifier or a
modifier function.

The optional $addscores argument stipulates that all merged 
feature scores are added to the merged version of the feature.

=back

=cut

sub self_overlap_merge {
    my $self      = shift ;
    my $tolerance = shift ; 
    my $strand    = shift ;
    my $group_tag = shift ;
    my $addscores = shift ;

    $group_tag = 0 if !defined($group_tag) ;

    GFF::TRACE(1,"Entering GFF::GeneFeatureSet::self_overlap_merge()\nSorting the input feature set first\n") ;

    # setup a sorted array of GFF::GeneFeature objects
    my @gf = $self->order_gf();

    GFF::TRACE(1,"Then, traverse the set...\n") ;
    # loop over each
    Outer: for(my $i = 0; $i < scalar(@gf) ; $i++){
      GFF::TRACE(1,'.')  if !($i%10) ;
      GFF::TRACE(1,"\n") if !($i%800) ;
	my $gf1 = $gf[$i];
	# skip any that have already been deleted
	next Outer unless $gf1;
	Inner: for(my $j = $i+1 ; $j < scalar(@gf) ;$j++){
	    my $gf2 = $gf[$j];
	    next Inner unless $gf2;
	    next Outer if $gf1->end() < $gf2->start() ; # short circuit non-overlap test!
#	    print "$i,$j\n";
	    if($gf1->overlap_merge($gf2,$tolerance,$strand,$group_tag,$addscores)){
		$gf[$j] = 0;
	        GFF::TRACE(2,"\nMerging($i,$j): ".$gf1->start." - ".$gf1->end."\n");
	    }
	}
    }

    # build new object from remaining list
    my $newSet = $self->new($self->version(),$self->region()) ;

    foreach my $gf (@gf){
	next unless $gf;
	$newSet->addGeneFeature($gf);
    }
    GFF::TRACE(1,"\nExiting GFF::GeneFeatureSet::self_overlap_merge()\n\n") ;

    return $newSet;
}

=pod

=over 4

=item intersect_range( $start, $end, $exact, $strand )

Method to generate a new GeneFeatureSet object
from all GeneFeatures in the invoking
GeneFeatureSet object which overlap the specified
$start and $end range By default, GeneFeatures
are included if they overlap at all with the
range. 

Use of the optional $exact boolean flag
(default: 'false') specifies that both the start
and end of GeneFeatures must lie completely 
within the specified range. Either, any overlap
with the range is a merge hit.

The optional $strand ('+','-' or '.') argument forces the
intersections to be strand sensitive (default: ignore
the strand value).

=back

=cut

sub intersect_range {
    my $self   = shift ;
    my $rstart = shift ;
    my $rend   = shift ;
    my $exact  = shift ;
    my $strand = shift ;
  
    croak "GeneFeatureSet::intersect_range(): undefined \$rstart?" if !defined $rstart ;
    croak "GeneFeatureSet::intersect_range(): undefined \$rend?"   if !defined $rend ;

    if( $rstart >= $rend ) {
	carp "GeneFeatureSet::intersect_range(): $rstart >= $rend?" ;
	return undef ;
    }

    $exact  = 0 if !defined $exact ;
    $strand = 0 if !defined $strand ;

    my ($hname,$hstart,$hend) = $self->region() ;
    my $newSet = $self->new($self->version(),$hname,$rstart,$rend) ;

    foreach my $gf (@{$self->{'feature'}}) {
	next if ($strand && ($gf->strand ne $strand)) ; # fail early on strand
        my $gStart = $gf->start() ;
        my $gEnd   = $gf->end() ;
        if( (!$exact && 
               ($gStart >= $rstart && $gStart <= $rend) ||
               ($gEnd   >= $rstart && $gEnd   <= $rend)
             ) or
            (  ($gStart >= $rstart && $gStart <= $rend) &&
               ($gEnd   >= $rstart && $gEnd   <= $rend)
             )
           ) 
	{ $newSet->addGeneFeature($gf) ; }
    }
    return $newSet ;
}

=pod

=over 4

=item intersect_overlap( $GeneFeatureSet2, $tolerance, $single, $strand, $soft, $tag )

Method to generate a new GeneFeatureSet object
based on intersection with second GeneFeatureSet
object ("GeneFeatureSet2"), on a GeneFeature by
GeneFeature basis (using the GeneFeature "match"
method). Refer to the GFF::GeneFeature::match()
method for information on $tolerance, $single and
$strand method arguments. 

The optional argument '$soft', when 'true' (non-null) 
specifies that intersection match copies of GeneFeatures 
from both the invoking and $GeneFeatureSet2 sets are
kept in the intersection set (i.e. the
intersection is based upon coordinate matches,
but a 'union' of the matching GeneFeatures). This
is useful for situations where one needs to use
the geometric intersection set in some context
where the feature identity must be preserved
(e.g. when taking a difference set of features
from a given which are not in the geometric
intersection set).

GFF Version 2: The optional '$tag' argument is a
[group] tag which is used to mark up the gene 
feature matches, using the match description.

=back

=cut

sub intersect_overlap {
    my $self      = shift ;
    my $other     = shift ;
    my $tolerance = shift ;
    my $single    = shift ;
    my $strand    = shift ;
    my $soft      = shift ;
    my $tag       = shift ; # optional match markup tag
    $tag = 0 if !defined($tag) ;

    GFF::TRACE(1,"Entering GFF::GeneFeatureSet::intersect_overlap()\n") ;

# new object from intersection
    my $newSet = $self->new($self->version(),$self->region()) ;
    my (%saved);
    my @gf1=$self->order_gf(0,1); # ascending by end
    my @gf2=$other->order_gf();   # ascending by start
    my $n=0 ;
# loop over each
    Outer: foreach my $gf1 (@gf1){
	++$n ;
        GFF::TRACE(1,'.')  if !($n++%10) ;
	GFF::TRACE(1,"\n") if !($n%600) ;
# loop over other 
	foreach my $gf2 (@gf2){
	    next Outer if $gf1->end() < $gf2->start() ; # short circuit non-overlap test!

	    my @string=$gf1->match($gf2,$tolerance,$single,$strand);
	    if($string[0]){
		unless($saved{$gf1}){
		    $gf1->tagValue($tag,0,$string[0]) if $tag ;
		    $newSet->addGeneFeature($gf1);
		    $saved{$gf1}=1;
		}
		unless(!$soft or $saved{$gf2}){
		    $gf2->tagValue($tag,0,$string[0]) if $tag ;
		    $newSet->addGeneFeature($gf2);
		    $saved{$gf2}=1;
		}
	    }
	}
    }
    GFF::TRACE(1,"Leaving GFF::GeneFeatureSet::intersect_overlap()\n") ;
    return $newSet ;
}

=pod

=over 4

=item intersect_overlap_merge( $GeneFeatureSet2, $tolerance, $strand, $tag )

Method to generate a new GeneFeatureSet object
based on intersection merger of features in a
second GeneFeatureSet object ("GeneFeatureSet2"), 
on a GeneFeature by GeneFeature basis (using the 
GeneFeature "overlap_merge" method; refer to the
GFF::GeneFeature::overlap_merge() method for 
information on $tolerance and $strand method 
arguments. Note: the method returns 'copied'
versions of the overlapping invoking features. 

GFF Version 2: The optional '$tag' argument is a
[group] tag which is used to mark up the gene 
feature matches, using the match description.
This value is also used to rename the <feature>
field of the overlap merged gene feature.

=back

=cut

sub intersect_overlap_merge {
    my $self      = shift ;
    my $other     = shift ;
    my $tolerance = shift ;
    my $strand    = shift ;
    my $tag       = shift ; # optional match markup/feature rename tag

    $tolerance = 0 if !defined($tolerance) ;
    $tag       = 0 if !defined($tag) ;
    $strand    = 0 if !defined($strand) ;

    GFF::TRACE(1,"Entering GFF::GeneFeatureSet::intersect_overlap_merge()\n") ;

# new object from intersection
    my $newSet = $self->new($self->version(),$self->region()) ;

    my (%saved);
    my @gf1=$self->order_gf(0,1); # ascending by end
    my @gf2=$other->order_gf();   # ascending by start
    my $n=0 ;
# loop over each
    Outer: foreach my $gf1 (@gf1){
	++$n;
        GFF::TRACE(1,'.')  if !($n%10) ;
	GFF::TRACE(1,"\n") if !($n%600) ;
	my $savecopy = 1 ;
# loop over other 
	foreach my $gf2 (@gf2){
            next Outer if $gf1->end() < $gf2->start() ; # quick shortcut non-overlap test

	    my $gfm = $gf1->overlap_merge($gf2,$tolerance,$strand,0,0,$savecopy) ;
	    if( $gfm ){
		if($tag) {
		    my $label = [$gf1->feature,$gf1->start,$gf1->end,
				 $gf2->feature,$gf2->start,$gf2->end] ;
		    $gfm->tagValueList($tag, $label, 1) ;
		    $gfm->feature($tag) ;
		}
		$newSet->addGeneFeature($gfm)if($savecopy);
		$gf1 = $gfm ; 
		$savecopy = 0 ; # no need for further copies of $gf1...
	    }
	}
    }
    GFF::TRACE(1,"Leaving GFF::GeneFeatureSet::intersect_overlap_merge()\n") ;
    return $newSet ;
}

=pod

=over 4

=item intersect_overlap_set( $GeneFeatureSet2, $tolerance, $single, $strand )

Check for exact match in position between two
GeneFeatureSet objects Refer to the
GFF::GeneFeature::match() method for information
on $tolerance, $single and $strand method
arguments.

=back

=cut

sub intersect_overlap_set {
    my $self = shift;
    my $other = shift;
    my $tolerance = shift;
    my $single = shift;
    my $strand = shift;
    if($self->count!=$other->count){
	return 0;
    }
    my($gff);
    $gff = $self->intersect_overlap($other,$tolerance,$single,$strand);
    if($gff->count!=$other->count){
	return 0;
    }
    return 1;
}

=pod

=over 4

=item intersect_overlap_matches( $GeneFeatureSet2, $tolerance, $single, $strand, $verbose, $file )

Method to generate a new GeneFeatureSet object
based on 'overlap' intersection with a second
GeneFeatureSet object. The new object contains
GeneFeatures from the first object, which have
matches to GeneFeatures in the second
GeneFeatureSet object. 

These matches are recorded
in the new GeneFeatureSet object and may be
retrieved by GeneFeature *Match*() methods. 

Refer to the GFF::GeneFeature::match() method for
information on $tolerance, $single and $strand
method arguments. If the $verbose flag is present
and equal to 1, then errors are reported to $file
(or STDOUT if $file is omitted). 

Note: a side effect of this method is to 
symmetrically record all GeneFeature matches 
in both the invoking GeneFeatureSet object
and the second $GeneFeatureSet2 object.

=back

=cut

sub intersect_overlap_matches {
    my $self      = shift;
    my $other     = shift;
    my $tolerance = shift;
    my $single    = shift;
    my $strand    = shift;
    my $verbose   = shift;
    my $file      = shift;
    my $annotate  = shift ;

    my (@string,$gf1,$gf2);
    $tolerance = 0 if !defined $tolerance ;
    $single    = 0 if !defined $single ;
    $strand    = 0 if !defined $strand ;
    if( defined $verbose ){
	$file = \*STDOUT if !defined $file ;
    }else{
        $verbose = 0 ;
    }
    $annotate    = 0 if !defined $annotate ;

# new object from the matches
    my $newSet = $self->new($self->version(),$self->region()) ;
    my (%saved);
    my @gf1=$self->order_gf(0,1); # ascending by end
    my @gf2=$other->order_gf();   # ascending by start
    my $n=0 ;
    Outer: foreach my $gf1 (@gf1){
#	print STDERR "\n### GF1:\n" ;
#	$gf1->dump(\*STDERR) ;
	++$n ;
	GFF::TRACE(1,"\n") if !($n%600) ;
        GFF::TRACE(1,'.')  if !($n%10) ;
        my $match_found = 0 ;
	foreach my $gf2 (@gf2){
#	    print STDERR "\n### GF2:\n" ;
#	    $gf2->dump(\*STDERR) ;
	    if($gf1->end() < $gf2->start()) {
#		print STDERR "Short circuited!\n" ;
                # short circuit non-overlap test!
		next Outer ;
            } 
	    @string = $gf1->match($gf2,$tolerance,$single,$strand) ;
	    if( $string[0] ){ # a match?
#
#               Record match found
#
#                print STDERR "\n\t\t*** gf1:", $gf1->Sequence, "(", $gf1->start(), ",", $gf1->end(),
#		               ") and gf2:",$gf2->Sequence,"(",$gf2->start(),",", $gf2->end(),
#		               ")\n" ;
#
		unless($saved{$gf1}){
		    $newSet->addGeneFeature($gf1) ;
		    $saved{$gf1}=1;
		}
		$gf1->addMatch( $gf2,  $string[1],  $string[2] ) ;
		$gf2->addMatch( $gf1, -$string[1], -$string[2] ) ; # reverse offset signs...
	    } elsif ($strand && $string[3]==3 && $verbose) {
#
#               Report strand orientation errors?
#
		print $file "  ERROR: DIRECTION? ".$string[4]."\n";
	    }
            # otherwise, just no match...ignore
	}
    }
    return $newSet;
}

=pod

=over 4

=item order_gf( $order, $by_end )

Method returns sorted array of GeneFeature
objects in the GeneFeatureSet. By default, the
array is sorted (low to high) by segment 'start'
coordinates. 

If the optional '$order' is defined and non-null,
then coordinate sorting is descending ('high to low'). 

If the optional '$by_end' argument is 
defined and non-null, then the feature 'end'
coordinate is used for the sort, instead
of the start coordinate.

=back

=cut

# Order low to high by start...
#
sub _by_start {
    $a->start <=> $b->start;
}
# Order low to high by end...
#
sub _by_end {
    $a->end <=> $b->end;
}

# Order high to low by start
#
sub _by_start_descend {
    $b->start <=> $a->start;
}

# Order high to low by end
#
sub _by_end_descend {
    $b->end <=> $a->end;
}

sub order_gf {
    my $self   = shift;
    my $order  = shift ;
    my $by_end = shift ;
    my @array;
    if(defined $order and $order) {
        if(defined $by_end and $by_end) {
	    @array = sort _by_end_descend   @{$self->{'feature'}} ;
	} else {
	    @array = sort _by_start_descend @{$self->{'feature'}} ;
	}
    } else {
	if(defined $by_end and $by_end) {
	    @array = sort _by_end   @{$self->{'feature'}} ;
	} else {
	    @array = sort _by_start @{$self->{'feature'}} ;
	}
    }
    return @array;
}

=pod

=over 4

=item order_gf_by_size( $descending )

Method returns sorted array of GeneFeature
objects in the GeneFeatureSet, sorted from
the smallest gene feature segment length to
the largest gene segment length.

If the optional boolean '$descending'
flag is given, the order is reversed, from
largest to smallest.

=back

=cut

sub _small_to_large {
    ( $a->end - $a->start ) <=> ( $b->end - $b->start ) ;
}

sub _large_to_small {
    ( $b->end - $b->start ) <=> ( $a->end - $a->start ) ;
}

sub order_gf_by_size {
    my $self       = shift;
    my $descending = shift ;
    my @array;
    if(defined $descending and $descending) {
	@array = sort _large_to_small @{$self->{'feature'}} ;
    } else {
	@array = sort _small_to_large @{$self->{'feature'}} ;
    }
    return @array;
}

=pod

=head1 Methods to do calculations on GFF::GeneFeatureSet objects

=over 4

=item count()

Method to count the number of GeneFeature objects
in a GeneFeatureSet object

=back

=cut

sub count {
    my $self = shift;
    my $gf;
    my $count;

    $count=0;
    foreach $gf (@{$self->{'feature'}}){
	$count++;
    }
    $count;
}

=pod

=over 4

=item strands()

Method return strands of the GeneFeature objects
which were found in the GeneFeatureSet object, in
the form of a string.

=back

=cut

sub strands {
    my $self = shift;
    my ($strand,%strand,$strandall,$gf);
# count occurances of each strand type
    foreach $gf (@{$self->{'feature'}}){
	$strand{$gf->strand}++;
    }
# return string of all strands that were found
    foreach $strand (sort keys %strand){
	$strandall.=$strand;
    }
    return $strandall;
}

=pod

=over 4

=item frames()

Method return frames of the GeneFeature objects
which were found in the GeneFeatureSet object, in
the form of a string.

=back

=cut

sub frames {
    my $self = shift;
    my ($frame,%frame,$frameall,$gf);
# count occurances of each frame type
    foreach $gf (@{$self->{'feature'}}){
	$frame{$gf->frame}++;
    }
# return string of all frames that were found
    foreach $frame (sort keys %frame){
	$frameall.=$frame;
    }
    return $frameall;
}

=pod

=over 4

=item minScore( $matched )

Method to find minimum score of the GeneFeature
members of a GeneFeatureSet object. If $matched
is defined and equal to 1, then the method only
evaluates objects which are matched to another
GeneFeature.

=back

=cut

sub minScore {
    my $self = shift;
    my $matched = shift;

    my ($gf,$minscore);
    $minscore = 0;
    foreach $gf (@{$self->{'feature'}}){
	if(!$matched || $gf->getMatches_logical){
	    $minscore = &GFF::GeneFeature::min($minscore,$gf->score);
	}
    }
    return $minscore;
}

=pod

=over 4

=item maxScore( $matched )

Method to find maximum score of the GeneFeature
members of a GeneFeatureSet object. If $matched
is defined and equal to 1, then the method only
evaluates objects which are matched to another
GeneFeature.

=back

=cut

sub maxScore {
    my $self = shift;
    my $matched = shift;

    my ($gf,$maxscore);
    $maxscore=0;
    foreach $gf (@{$self->{'feature'}}){
	if(!$matched || $gf->getMatches_logical){
	    $maxscore = &GFF::GeneFeature::max($maxscore,$gf->score);
	}
    }
    return $maxscore;
}

=pod

=over 4

=item sumScore( $matched, $unweighted)

Method to find (sequence weighted) sum of 
all scores of members of a GeneFeatureSet object.
If $matched argument is set, then only evaluates
for objects which are matched to other GFF::GeneFeatures
In a scalar context, only the sum is returned. In a list context,
the sum plus number of items added is returned. If no items were counted
then 0 is returned. The '$unweighted' argument, when non-NULL, stipulates that
a simple sum of scores rather than a feature length weighted sum of 
scores should be returned.

=back

=cut

sub sumScore {
    my $self   = shift ;
    my $matched = shift ;
    my $unweighted = shift ;

    my ($gf,$sscore,$nscore,$n);
    $sscore=$nscore=0;
    foreach $gf (@{$self->{'feature'}}){
	if(!$matched || $gf->getMatches_logical){
	    if($unweighted) {
		$n = 1 ;
	    } else {
	        $n = (($gf->end)-($gf->start))+1 ;
	    }
	    $nscore += $n;
	    $sscore += ($gf->score) * $n;
	}
    }
    if($nscore){
        my $context = wantarray ;
        if(!defined $context) {
	    croak "GFF::GeneFeatureSet::sumScore: invalid void context!" ; # void context
	} elsif($context){
	    return ($sscore,$nscore); # array context
	} else {
	    return $sscore ;          # scalar context
	}
    } else {
	return 0 ;
    }
}

=pod

=over 4

=item avScore( $matched, $unweighted)

Method to find average score of the GeneFeature
members of a GeneFeatureSet object. 

If $matched is defined and equal to 1, then the 
method only evaluates objects which are aligned to 
another GFF::GeneFeature.

The '$unweighted' argument, when non-NULL, stipulates that
a simple sum of scores rather than a feature length weighted sum of 
scores should be taken for the computation of the average.

=back

=cut

sub avScore {
    my $self     = shift ;
    my $matched  = shift ;
    my $unweighted = shift ;

    my ( $sscore, $nscore ) = $self->sumScore($matched,$unweighted) ;

    if( defined $nscore && $nscore ){
	my $n = sprintf("%10.1f",$sscore/$nscore) ;
	return $n ;
    }else{
	return 0 ;
    }
}

=pod

=over 4

=item scoreRange( $matched )

Method to return the minimum and maximum score
information for members of a GeneFeatureSet
object. Returns a list of minimum and maximum
scores in an array context. Returns the spread
(difference) between maximum and minimum in a
scalar context. If $matched is defined and equal
to 1, then the method only evaluates objects
which are matched to another GeneFeature.

=back

=cut

=pod

=over 4

=item maxScoreGf( \&filter )

Method to return the GeneFeature object having
the maximum score in a GeneFeatureSet object. If
the optional, user-defined &filter predicate
function (taking a GeneFeature object reference
as its argument; returning a value of 1 for
inclusion or 0 for exclusion) is provided, then
only the &filter defined subset of GeneFeatures
is assessed for scores.

=back

=cut

=pod

=over 4

=item max_min_range( $text )  -- Tim's old version (deprecated)

Method to find minimum and maximum range of all
members in a GeneFeatureSet object. The method 
returns (max_range, min_range) list pair unless
the optional non-null '$text' argument is
provided, in which case, the method returns the
string "min_range-max_range".

=back

=cut

sub max_min_range {
    my $self = shift;
    my $text = shift;
    my ($gf,$minrng,$maxrng);
    $maxrng=0;
    $minrng=1000000000;
    
    # EB fixes. Thrown warning on usage and changed name space on 
    # function call for min.max range.

    print STDERR "Warning: using old max_min_range. Should convert to min_max_range as soon as possible!";
    
    foreach $gf (@{$self->{'feature'}}){
	$maxrng=&GFF::GeneFeature::max($maxrng,$gf->end);
	$minrng=&GFF::GeneFeature::min($minrng,$gf->start);
    }
    if($text){
	return "$minrng-$maxrng";
    }else{
	return ($maxrng,$minrng);
    }
}

=pod

=over 4

=item min_max_range( $filter )

Method to find minimum and maximum range of all
members in a GeneFeatureSet object. In a list
context, the method returns (min_range,
max_range) list pair; in a scalar context, the
range difference value, max_range minus
min_range, is returned. 

An optional boolean predicate $filter function 
(taking a GeneFeature reference as its argument 
and returning 1 == inclusion, 0 == exclusion) 
may be used to limit range calculation to 
a specific subset of GeneFeatures in the 
invoking GeneFeatureSet object.

=back

=cut

sub min_max_range {
    my $self     = shift ;
    my $function = shift ;

    my ($gf,$minrng,$maxrng);
    $maxrng = 0 ;
    $minrng = 1000000000 ;

    foreach $gf (@{$self->{'feature'}}) {
        next if($function && !&$function($gf)) ;
	$maxrng = &GFF::GeneFeature::max($maxrng,$gf->end);
	$minrng = &GFF::GeneFeature::min($minrng,$gf->start);
    }
    my $context = wantarray ;
    if(!defined $context) {
	croak "GFF::min_max_range: invalid void context!" ; # void context
    }elsif($context){
	return ($minrng, $maxrng); # array context
    }else{
	return $maxrng - $minrng ; # scalar context
    }
}

=pod

=over 4

=item min_max_range_homol()

Method to find maximum and minimum range, of HomolGeneFeature members
only, in a GeneFeatureSet object. In a list context, the method
returns (min_range, max_range) list pair; in a scalar context, the
range difference value, max_range minus min_range, is returned.

=back

=cut

sub min_max_range_homol{
    my $self = shift;
    my ($gf,$minrng,$maxrng);
    $maxrng=0;
    $minrng=1000000000;
    foreach $gf (@{$self->{'feature'}}){
	if(ref($gf) ne 'GFF::HomolGeneFeature'){
	    next ; # rbsk: just ignore non-HomolGeneFeatures rather than return ''?
	}else{
	    # start2 and end2 can be in wrong order
	    my $start2=$gf->start2;
	    my $end2=$gf->end2;
	    if($start2>$end2){
		if($maxrng<$start2){$maxrng=$start2;}
		if($minrng>$end2){$minrng=$end2;}
	    }else{
		if($maxrng<$end2){$maxrng=$end2;}
		if($minrng>$start2){$minrng=$start2;}
	    }
	}
    }
    my $context = wantarray ;
    if(!defined $context) {
	croak "GFF::min_max_range: invalid void context!" ; # void context
    }elsif($context){
	return ($minrng, $maxrng); # array context
    }else{
	return $maxrng - $minrng ; # scalar context
    }
}

=pod

=over 4

=item direction_homol()

Method to find direction of HomolGeneFeature objects in a
GeneFeatureSet object.  If start2<end2 for all then returns 1, if
start2>end2 for all then returns -1, if mixture then returns 0.

=back

=cut

sub direction_homol{
    my $self = shift;
    my ($gf,$dir);
    foreach $gf (@{$self->{'feature'}}){
	my $start2=$gf->start2;
	my $end2=$gf->end2;
	if($start2>$end2){
	    if(!$dir){
		$dir=-1;
	    }elsif($dir==1){
		$dir=0;
		last;
	    }
	}else{
	    if(!$dir){
		$dir=1;
	    }elsif($dir==-1){
		$dir=0;
		last;
	    }
	}
    }
    return $dir;
}

=pod

=over 4

=item start_range()

Method to find the minimum start coordinate among
members of a GeneFeatureSet object.

=back

=cut

sub start_range {
    my $self = shift;
    my ($gf,$minrng);
    $minrng=1000000000;
    foreach $gf (@{$self->{'feature'}}){
	$minrng = &GFF::GeneFeature::min($minrng,$gf->start);
    }
    return $minrng;
}

=pod

=over 4

=item end_range()

Method to find the maximum end coordinate among
members of a GeneFeatureSet object.

=back

=cut

sub end_range {
    my $self = shift;
    my ($gf,$maxrng);
    $maxrng=0;
    foreach $gf (@{$self->{'feature'}}){
	$maxrng = &GFF::GeneFeature::max($maxrng,$gf->end);
    }
    return $maxrng;
}

=pod

=over 4

=item remap( $offset, $rcend )

Method to add an $offset amount to the start and end coordinates of
every GeneFeature in a GeneFeatureSet object (uses
GFF::GeneFeature->remap()). Note: the start and end of the original
gene features (not copies thereof) of the the gene feature set are
changed.

If the optional $rcend argument is defined 
and non-NULL, then the remapping is presumed to be a 
reverse complementation operation on all the gene features
and the '##sequence-region' bounds of the set.

If $rcend > 0 then this $rcend value is such that it is assumed to 
be the host set coordinate directly remapped onto the coordinate
position $offset+1. 

If $rcend < 0, then the $rcend is set to equal the invoking 
GFF::GeneFeatureSet's '##sequence-region' end coordinate.

=back

=cut

sub remap {
    my $self   = shift ;
    my $offset = shift ; 
    my $rcend  = shift ;

    if( !defined $offset ) {
	confess("GFF::GeneFeatureSet::remap() undefined offset!") ;
    } elsif( $offset !~ /^([\+-]?\d+)$/ ) {
	confess("GFF::GeneFeatureSet::remap() non-numeric offset: $offset!") ;
    } elsif (!$offset and !defined($rcend)) {
        return ; # 0 offset is a do nothing, unless reverse complementing
    } else {
       # Remap the sequence-region
       my($seqname,$start,$end) = $self->region() ;
       if(defined($rcend)) {
	   $rcend = $end if $rcend<0 ;
	   $end   = $rcend-$start+1 ;
	   $start = $rcend-$end  +1 ;
       }
       if(defined($start) and defined($end)) {
	   $self->region( $seqname, $offset+$start, $offset+$end) ;
       }
       # Remap the objects
       my $gf;
       foreach $gf ( @{$self->{'feature'}} ) {
	   $gf->remap($offset, $rcend);
       }
    } 
}

=pod

=over 4

=item mask_length()

Calculates the sum of all segment lengths
(end-start+1) of all the GeneFeatures in the
GeneFeatureSet object. Method assumes
non-overlapping features (i.e. that a
self_overlap_merge() method call has been done
first).

=back

=cut

sub mask_length {
    my $self = shift;
    my ($gf,$sumrng);
    $sumrng=0;
    foreach $gf (@{$self->{'feature'}}){
	$sumrng+=($gf->end)-($gf->start)+1;
    }
    return $sumrng;
}


=pod

=over 4

=item mask_length_true($max,$min)

Calculates true coverage of all the GeneFeatures in the GeneFeatureSet
object.  Need to supply maximum range ($max,$min). $min defaults to 1
if not specified.

=back

=cut

sub mask_length_true {
    my $self = shift;
    my $max = shift;
    my $min = shift;
    $min = 1 if !defined($min) ;
    my $string=' ' x ($max-$min+1);
    foreach my $gf (@{$self->{'feature'}}){
	my $end=($gf->end)-$min+1;
	my $start=($gf->start)-$min+1;
	my $len=$end-$start+1;
	substr($string,$start,$len)='x' x $len;
    }
    $string=~s/ //g;
    return length($string);
}


=pod

=over 4

=item lengthStats($filter)

Calculates the mean length 
of segment lengths (end-start+1) of all the 
GeneFeatures in the GFF::GeneFeatureSet object.

In a scalar context, this method simply returns
the mean value of feature lengths. 

In a list context, the mean, the minimum length,
the maximum length, the standard
deviation and underlying variable values of
total number of features, total feature lengths,
and sum of feature lengths squared, plus
a reference to a (sparse) hash of tallies for each 
non-zero length class (keyed by lengths) are 
also returned, in that order respectively.

An optional boolean predicate $filter function 
(taking a GeneFeature reference as its argument 
and returning 1 == inclusion, 0 == exclusion) 
may be used to limit the calculation to a 
specific subset of GeneFeatures in the 
invoking GeneFeatureSet object.

The method returns 'undef' if no features are
found upon which to compute statistics.

=back

=cut

sub lengthStats {
    my $self   = shift;
    my $filter = shift ;

    my ($gf,$n,$minx,$maxx,$x,$x2,%len) ;
    $minx=100000000;
    $maxx=0 ;
    my $context = wantarray ;
    my $full_stats = (defined($context) && $context) ;
    $n = $x = $x2 = 0;
    foreach $gf (@{$self->{'feature'}}){
	next if($filter && !&$filter($gf)) ;
	$n++ ;
	my $xn =($gf->end)-($gf->start)+1 ;
	$minx=$xn if ($xn<$minx) ;
	$maxx=$xn if ($xn>$maxx) ;
	$x  += $xn ;
	if($full_stats) {
	    $x2  += $xn * $xn ;
	    if(!exists($len{"$xn"})) {
		$len{"$xn"} = 1 ;
	    } else {
		$len{"$xn"}++ ;
	    }
	}
    }
    return undef if !$n ;

    my $mean = sprintf("%.2f",$x/$n) ; # probably adequate precision!
    if(!defined $context) {
	croak "GFF::GeneFeatureSet::lengthStats(): invalid void context!" ; # void context
    } elsif ($context){
	my $SD = sprintf("%.2f",( $x2/$n - $mean**2 )**0.5) ;
	return ($mean, $SD, $n, $minx, $maxx, $x, $x2, \%len); # array context
    } else {
	return ($mean); # scalar context
    }
}

=pod

=over 4

=item shared_matches( $GeneFeatureSet2 )

Method to return a number of (overlap) matching
GeneFeature objects between the invoking
GeneFeatureSet object and a second GeneFeatureSet
object ("$GeneFeatureSet2").

=back

=cut

#
# Method to dump out a GFF object, via method to dump_msp GFF::HomolGeneFeature object.
#
    
sub dump_msp {
    my $self    = shift ;
    my $file    = shift ;

    foreach my $gf (@{$self->{'feature'}}){
	if(ref($gf) eq 'GFF::HomolGeneFeature'){
	    $gf->dump_msp($file);
	}
    }
    return 1;
}


#
# Deletes version 2 tags...
#
sub deleteTag {
    my $self  = shift ;
    my $tag   = shift ;
    foreach my $gf (@{$self->{'feature'}}) {
	$gf->deleteTag($tag) ;
    }
}



#
# Method to extract gene features lying or 
# overlapping a subrange of the calling object
# returns a new GeneFeatureSet object containing these features
# if undef or null, $start defaults to the start of the sequence
# if undef or null, $end   defaults to the end of the invoking sequence
#
sub subrange() {
    my $self  = shift ;
    my $start = shift ;
    my $end   = shift ;

    $start = $self->start_range() if !(defined $start and $start) ;
    $end = $self->end_range() if !(defined $end and $end) ;

    my $gffs = $self->new($self->version(),$self->region()) ; 
}


#
# Method to test if a GFF::GeneFeature object is a member of an existing GeneFeatureSet object
# (note - comparison is by reference only - does not check if 2 GFF::GeneFeature objects
#  contain the same data!)
#

sub member {
    my $self    = shift ;
    my $element = shift ;

    my ($gf,$flag);
    $flag=0;
    foreach $gf (@{$self->{'feature'}}){
	if($gf == $element){
	   $flag = 1;
	   last;
       }
    }
    return $flag;
}

#
# Method to return the minimum and maximum score information 
# for members of a GeneFeatureSet object.
# Returns a list of minimum and maximum scores in an array context.
# Returns the spread (difference) between maximum and minimum in a scalar context
# If $matched argument is set, then only evaluates
# for objects which are matched to other GFF::GeneFeatures
#

sub scoreRange {
    my $self = shift;
    my $matched = shift;

    my ($gf,$minscore,$maxscore);
    $maxscore = 0 ;
    $minscore = 1000000; # impossibly high number...
    foreach $gf (@{$self->{'feature'}}){
	if(!$matched || $gf->getMatches_logical){
	    $minscore = &GFF::GeneFeature::min($minscore,$gf->score);
	    $maxscore = &GFF::GeneFeature::max($maxscore,$gf->score);
	}
    }
    my $context = wantarray ;
    if(!defined $context) {
	croak "GFF::GeneFeatureSet::scoreRange: invalid void context!" ; # void context
    }elsif($context){
	return ($minscore,$maxscore); # array context
    }else{
	return $maxscore - $minscore ; # scalar context
    }
}

#
# Method to return GF object with max score in a GeneFeatureSet object.
#

sub maxScoreGf {
    my $self = shift;
    my $function = shift;
    my ($gf,$mxscore,$mxgf);
    $mxscore=0;
    foreach $gf (@{$self->{'feature'}}){
	if(!$function || &$function($gf)){
	    if($gf->score>$mxscore){
		$mxscore=$gf->score;
		$mxgf=$gf;
	    }
	}
    }
    return $mxgf;
}


#
# Version 2 GFF: adds specified tag-value pair 
# to every GFF::GeneFeature in the invoking GeneFeatureSet object
#
sub addGroupValue {
    my $self = shift ;
    my $tag  = shift ;
    my $value = shift ;
}

#
# Method to add a label to each GFF::GeneFeature in this GeneFeatureSet point back to this GeneFeatureSet
#

sub label {
    my $self = shift;
    my $membergroup = shift;
    my ($gf);
    foreach $gf (@{$self->{'feature'}}){
	$gf->addMember($self,$membergroup);
    }
}

#
# Method to make a list of the corresponding list of GeneFeatureSet's 
# that the GFF::GeneFeatures in this GeneFeatureSet are matched with
# by virtue of a common $membergroup label (e.g. true matching predicted 'Gene')
#

sub label_pair {
    my $self = shift;
    my $membergroup = shift;
    my ($gf1,$gf2,%gff,$gff,@gff );

    foreach $gf1 (@{$self->{'feature'}}) {
        next if !$gf1->getMatches_logical() ;  # check for matches...

        my $match = $gf1->getMatches() ; # is a reference to a hash of matching GFF::GeneFeatures
	foreach $gf2 (keys %{$match}) {
	     $gff = $match->{$gf2}->[0]->getMember($membergroup) ; # retrieve the corresponding GFF::GeneFeatureSet
	     if($gff){
		 unless($gff{$gff}){
		     push(@gff,$gff); # keep a (non-duplicate) list of such objects to return
		 }
		 $gff{$gff}++;
	    }
	}
    }
    return @gff;
}

#
# Method to return a number of matched GFF::GeneFeature objects between 2 GeneFeatureSet's
#

sub shared_matches {
    my $self  = shift ;
    my $other = shift ;

    my ($gf,%gf,$match,$n) ;
#
#   Hash remember all the 'other' GFF::GeneFeature references
#
    foreach $gf (@{$other->{'feature'}}){
	$gf{$gf} = 1 ;
    }

#
#   Then compare them to all of my $self GFF::GeneFeatures
#
    $n = 0 ;
    foreach $gf (@{$self->{'feature'}}) {
        my $match = $gf->getMatches() ;  # is a reference to a hash...
	foreach my $gfm (keys %{ $match }){
	    if( $gf{ $gfm } ){
		$n++;
	    }
	}
    }
    $n;
}

=pod

=over 4

=item score( $tolerance, $nep, $net, \*VERBOSE )

Method returns an (4 x 3) array of scores for a
set of GeneFeature objects which contain
(overlap) matches to other GeneFeatures. This
scoring provides the number of features (N) and
fractional accuracy (specificity) and coverage
(sensitivity) for exact matches, overlaps, 5' and
3' alignments of the gene feature matches,
respectively. The $tolerance indicates how
precise the boundaries must match to count as a
match to be scored. The $nep is the number of
predicted exons and the $net is the number of
true exons. If a score is "infinite" due to
division by zero, a -1 is returned for that
value. The optional '\*VERBOSE' argument, if
defined, is a file device for dumping of detailed
(signed) match offset statistics. 

=back

=cut

sub score {
    my $self      = shift ;
    my $tolerance = shift ;
    my $nep       = shift ;
    my $net       = shift ;
    my $Offset5   = shift ;
    my $Offset3   = shift ;

    $Offset5 = 0 if !defined $Offset5 ;
    $Offset3 = 0 if !defined $Offset3 ;

    my($no,$ne,$n5,$n3) ;
    my($nenep,$nonep,$n5nep,$n3nep,
       $nenet,$nonet,$n5net,$n3net);
    my(%Offset5, %Offset3) ;
#
#   For each GFF::GeneFeature in the set, 
#
    foreach my $gf (@{$self->{'feature'}}){
        next if !($gf->getMatches_logical()) ; # skip if does not have matches

	$no++ ;  # count any match at all as an overlap hit

        my ($m3,$m5,$me) ;
        $m3 = $m5 = $me = 0 ;
        my $match = $gf->getMatches() ; # is a reference to a hash 
        foreach my $gfm (keys %{$match}) {
	   my $exact = 0 ;
           if($Offset5) {
	       my $hash = $match->{$gfm}->[1] ;
	       if(!defined $Offset5->{"$hash"}){
		   $Offset5->{"$hash"} = 1 ;
	       } else {
		   $Offset5->{"$hash"}++ ;
	       }
	   }
           if($Offset3) {
	       my $hash = $match->{$gfm}->[2] ;
	       if(!defined $Offset3->{"$hash"}){
		   $Offset3->{"$hash"} = 1 ;
	       } else {
		   $Offset3->{"$hash"}++ ;
	       }
	   }
	   if(abs($match->{$gfm}->[1]) <= $tolerance){ # $match->{$gfm}->[1] == e5'
		$m5++; 
		$exact++ ;
	   }
	   if(abs($match->{$gfm}->[2]) <= $tolerance){ # $match->{$gfm}->[2] == e3'
		$m3++;
		$exact++ ;
	   }
	   if($exact == 2){
	       $me++;
	   }
        }
#
#       Count a GFF::GeneFeature overlaps only once...
#
        $n3++ if $m3 ;
        $n5++ if $m5 ; 
        $ne++ if $me ;
    }
    if($nep==0){
	$nenep= -1;
	$nonep= -1;
	$n5nep= -1;
	$n3nep= -1;
    }else{
	$nenep=$ne/$nep;
	$nonep=$no/$nep;
	$n5nep=$n5/$nep;
	$n3nep=$n3/$nep;
    }
    if($net==0){
	$nenet= -1;
	$nonet= -1;
	$n5net= -1;
	$n3net= -1;
    }else{
	$nenet=$ne/$net;
	$nonet=$no/$net;
	$n5net=$n5/$net;
	$n3net=$n3/$net;
    }
    $ne,$nenep,$nenet,
    $no,$nonep,$nonet,
    $n5,$n5nep,$n5net,
    $n3,$n3nep,$n3net;
}

1;  # says use was ok

=pod

=head1 Revision History

2.101 (18/05/2000) rbsk: in 'complement' method, modified $tag management to
                        ...the complementary features in the form [Between_$tag $feature1 $feature2]...

2.100 (15/05/2000) rbsk: annotate_overlaps() method created and refined

2.099 (04/04/2000) rbsk: repaired bug in Perl regex for comment parsing

2.098 (31/03/2000) rbsk: featureLengthStats() now also returns min and max length stats.

2.097 14/03/2000 - rbsk: - copy method versioning fixed
2.096 02/2000 - rbsk:    - minor bug repairs to region() meta-comment parsing in

2.095 21/01/2000 - rbsk: - added reverse complementation to remap() method
                           dump() returns '#' of features dumped.
                           added source_version() meta-comment method

2.094 (29/12/99) - rbsk: - intersect_range() method should return local start,end range in GFF header

2.093 (24/11/99)- rbsk: - fixed logical errors and performance issues in intersect_overlap*() methods
                          Thanks to Alessandro Guffante for bringing these to my attention ;-)
                        - order_gf() argument semantics changes: the single $reverse argument is replaced by
                          $descending, to designate 'high to low' sorting by coordinates,
                          $end argument added to force usage of the gene feature end coordinate instead
                          of the start coordinate. The old '$reverse' argument is thus replaced by
                          $descending == $end == non-null; 
                          the method still defaults to 'start coordinate, ascending sort';

2.092 (19/11/99)- rbsk: - cluster() method, for matches, the comparison function can now
                          return a non-null $name string for unique labelling of the cluster
                          which becomes the ##sequence-region $name for the cluster
                        - dump_matches(): $show_nomatches argument added; Normally, 
                          only features with matches are dumped. The optional boolean 
                          flag '$show_nomatches' when defined and non-null,
                          directs that 'no match' records are reported too.

2.091 (10/11/99)- rbsk: - major algorithmic performance enhancement of pairwise overlap methods!
2.090 (5/11/99) - rbsk: - added the 'read_header()' method
2.089 (30/10/99)- rbsk: - added the 'pipe()' method
2.088 (21/10/99)- rbsk: - added GFF complement() method
2.087 (13/10/99)- rbsk: - added lengthStats() method
2.086 (12/10/99)- rbsk: - added method nextGeneFeature() (FIFO inverse of addGeneFeature()).

2.085 (3/10/99) - rbsk: - $group_tag in self_overlap_merge() access CODE ref
                          and subsumed into GFF::GeneFeature::overlap_merge() method.
2.084 (30/9/99) - rbsk: - $tag argument in dump()
                        - $group_tag argument moved over in self_overlap_merge()
                        - using order_gf() in *overlap*() methods (and run tracing)
2.083 (27/9/99) - rbsk: - $strand argument in self_overlap_merge() && intersect_range()
                        - make $strict mode in self_overlap() strand sensitive
2.082 (21/9/99) - rbsk: - added optional '$tag' argument to intersect_overlap()
                        - created intersect_overlap_merge() method
                        - created the deleteTag() method
2.081  (9/9/99) - rbsk: order_by_gf_size() method added
    3/9/99 - rbsk: $exclude to addGeneFeature() method; exclude() method added
   31/8/99 - rbsk: self_overlap_merge() method: $group_tag specification
                   allows for recording of merged features (by $group_tag)
                   optional '$tolerance' value 
                   provides for overlap merge where the two features lie
                   within $tolerance base pairs of each other
   21/7/99 - rbsk: rewriteField() $target matches made case insensitive
                   and framed by /^...$/
   12/7/99 - rbsk: creation from miscellaneous GFF analysis code; transferred
                   methods makeGenes(), constructGene() and mRNA from 
                   the GFF::GeneFeatureSet to new GFF::Analysis module
    7/7/99 - rbsk: read() $echo argument added to provide user feedback during
                   reading in of large GFF files...
    5/7/99 - rbsk: transcript() method renamed to mRNA() - biologically more accurate :-)
   28/6/99 - rbsk: transcript() method added
   14/6/99 - rbsk: recoded self_overlap() slightly to account for undefined group_values
   24/5/99 - rbsk: $tag argument to self_overlap()
                   new() region() setting bug fixed
   21/5/99 - rbsk: revised the rewriteField() method (see above) ;  Note that the  
                   $group_tag/value & $saveold arguments were removed from this method
   14/5/99 - rbsk: created the 'copy()' method
    7/5/99 - rbsk: reinserted Tim's old max_min_range() method, 
                   for backwards compatibility (deprecated?)
    6/5/99 - rbsk: added '$verbose' argument to score method()
   28/4/99 - rbsk: renamed GFF.pm to GFF::GeneFeatureSet.pm
   27/4/99 - rbsk: added '$filter' argument to min_max_range()
                   added minScore() and scoreRange() methods
   26/4/99 - rbsk: renamed max_min_range() to min_max_range()
   23/4/99 - rbsk: addGeneFeatureSet() documentation fixed...
   21/4/99 - rbsk: intersect_range(), containsMembers(), and getAllMembers() methods added
                   rewriteField(): added '$saveold' argument.
                   $copy argument added to group(), addGeneFeature(), addGeneFeatureSet() methods
   19/4/99 - rbsk: deletion bug in self_overlap() fixed...
                   GeneFeatureSetPair functionality (score() method et al.) merged with GeneFeatureSet.pm
   16/4/99 - rbsk: added $exact flag to self_overlap() method; overlap method debugged too
                   added $soft flag to intersect_overlap()
   1/4/99  - rbsk: order_gf: optional $reverse argument to reverse sorting order to
                   high to low by segment end coordinates.
   26/3/99 - rbsk: new method 'rewriteField()' added.
   23/3/99 - rbsk: $strict flag added to self_overlap() method
   16/3/99 - rbsk: GeneFeatureSet objects now subclassed from GeneFeatureSetObject class;
                   moved "version()" method from GeneFeatureSet.pm into
                   the new base class GeneFeatureSetObject.pm
   25/2/99 - rbsk: Extensively revised and improved the documentation
                   Added Version 2 GeneFeatureSet code, including &version class function
       Also:
             Converted all "constructor" type methods in all GeneFeatureSet libraries to class methods
             (i.e. must now be invoked as class->new*(args) or "new class args")
             Standardize all file glob arguments to \*FILE references
             Use "croak" instead of "die" in theFeature
             Added GeneFeatureSet "intersection" set operation;
             Renamed "intersect_not" to "difference"
             Default $type for read_parse is now "GeneFeature"
             Added "frames" method; rename "strand" to "strands"
             max_min_range() modified to sense array v/s scalar context
             max_min_range_homol() just ignores non-HomolGen

=cut

