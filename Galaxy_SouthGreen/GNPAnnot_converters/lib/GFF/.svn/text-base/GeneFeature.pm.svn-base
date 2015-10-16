=head1 NAME

GFF::GeneFeature.pm - Perl object module extension for GFF Gene Features

=head1 SYNOPSIS

use GFF ;  # which contains an implicit 'use GFF::GeneFeature'

=cut

package GFF::GeneFeature;

=head1 AUTHORS

Copyright (c) 1999 ,2000
Created by Tim Hubbard th@sanger.ac.uk
Augmented by Richard Bruskiewich rbsk@sanger.ac.uk

Sanger Centre, Wellcome Trust Genome Campus, Cambs, UK
All rights reserved.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation.

=head1 DESCRIPTION

GFF::GeneFeature (derived from the GFF.pm http://www.sanger.ac.uk/Software/GFF/GFF.shtml base class) 
is a class of Perl object describing a single
gene feature record ("line") in the Gene Finding Feature ("GFF") format. A GFF::GeneFeatureSet 
http://www.sanger.ac.uk/Software/GFF/GeneFeatureSet.shtml
is  the container object for a set of GFF::GeneFeature objects. 
A GFF::HomolGeneFeature http://www.sanger.ac.uk/Software/GFF/HomolGeneFeature.shtml 
is a homology-specific gene feature object class derived from the GFF::GeneFeature object class.

=head2 How to Read Method Protocols

Normal Perl data type notations are used for argument declarations in the method protocols. A
backslash denotes argument passing by reference. Arguments shown in italics are optional
parameters to the method calls. Class methods are invoked using the 'class->method(args)' or
'method class args' Perl call formats.

=head1 SOURCE CODE

The most current release of the Perl source code for this module is available at
 ftp://ftp.sanger.ac.uk/pub/gff/ .
Bug reports may be submitted to Richard Bruskiewich at rbsk@sanger.ac.uk . 

Future releases may be issued as Bioperl archived modules.

=cut

use Carp;
use strict;
use vars qw($AUTOLOAD @ISA @EXPORT @EXPORT_OK $VERSION);
require Exporter ;
use GFF::HomolGeneFeature ;

$VERSION = '3.005' ;
#
# @ISA has our inheritance.
#
@ISA = qw( Exporter GFF ); 
#
# Place functions/variables you want to *export*, ie be visible from the caller package into @EXPORT_OK
#
@EXPORT    = qw();
@EXPORT_OK = qw();

=pod

### old implementation ###

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
);

=cut

################### 02/02/2000 attempt at increasing memory efficiency #####################

# Concept: 1. indirection for string storage: 
#
#		Hash all seqname, source, feature, and tag strings to a unique index number
#		which the gene feature fields will index.  Treat the 'Sequence' tag like a special field tag.
#
#           2. convert genefeature objects to a simple array too

{
	my %SymTab ;

	# returns a 'unique identifier' reference 
	# to a given scalar string
	# used in seqname, source, feature, or tag fields
	sub uid {
	    my $string = shift ;
	    if(!exists($SymTab{$string})) {
	        $SymTab{$string}=\$string ;
	    }
	    return $SymTab{$string} ;
	}
	    

}
    
# Fixed length array indices             
use constant VERSION  => 0 ;
use constant SEQNAME  => 1 ;
use constant SOURCE   => 2 ;
use constant FEATURE  => 3 ;
use constant GFSTART  => 4 ;
use constant GFEND    => 5 ;
use constant SCORE    => 6 ;
use constant STRAND   => 7 ;
use constant FRAME    => 8 ;
use constant GROUP    => 9 ; # GFF Version 1
use constant SEQUENCE => 9 ; # GFF Version 2
use constant TAGS     => 10 ;
use constant MATCHES  => 11 ;
use constant MEMBER   => 12 ;
use constant COMMENT  => 13 ;

{
    sub _debug {
	my $self = shift ;
	my $pos  = shift ;
	print STDERR "GFF::GeneFeature Dump" ;
	print STDERR "\@$pos" if $pos ;
	print STDERR ":\n" ; 
	for(my $i=VERSION;$i<=COMMENT;++$i) {
	    print STDERR "\t$i: ",(defined($self->[$i])?$self->[$i]:'undefined'),"\n" ;
	}
    }
}

################### 02/02/2000 attempt at increasing memory efficiency #####################

=head1 GFF::GeneFeature Construction/Input Methods

For all construction methods, an optional
"$version" argument may be given which sets the
created object to a specified GFF specification
version. If this argument is not given, then the
GFF class (package) default version is used (see
GFF->version()).

=over 4

=item new( $version )

Class method to construct a new, empty
GFF::GeneFeature object.

=back

=cut

sub new {
    my $ref = shift;
    my $class = ref($ref) || $ref;
    my $version = shift ;
    
    if(!defined($version)) {
	if(ref($ref) && GFF::GeneFeature->verify($ref)) {
	   $version = $ref->version ;
        } else {
	   $version = GFF->version ;
	}
    }

    my $self=[] ; # 02/2000 data representation: ARRAY
    $self->[VERSION] = $version ;

    # set some things to their default GFF null values
    $self->[SCORE]   = "." ; 
    $self->[STRAND]  = "." ; 
    $self->[FRAME]   = "." ; 

    # defer allocation of remaining data structure

    bless $self, $class ;
    return $self;
}

=pod

=over 4

=item new_from_line( $line, $version )

Class method to parse a gene feature $line string
into a GFF::GeneFeature object (creates and
returns the object reference).

=back

=cut

sub new_from_line {
    my $ref     = shift ;
    my $string  = shift ;
    my $version = shift ;

    my $self ;
    my @array ;

    chomp($string);

    $self = $ref->new($version) ;

    # Tabs must delimit fields in all GFF Versions
    # Note: double quoted free text string values in the optional [attributes] field
    #       can have tabs encoded as C/UNIX/Perl style backslashed escaped '\t' but these 
    #       should be ignored by split? 
    @array = split(/\t/,$string);

    $self->seqname(shift(@array));
    $self->source(shift(@array));
    $self->feature(shift(@array));
    $self->[GFSTART]=shift(@array);
    $self->[GFEND]=shift(@array);
    $self->[GFEND] >= $self->[GFSTART]
       || carp "Illegal GeneFeature: end(",
	       $self->[GFEND],") < start(",$self->[GFSTART],
               ") from '$string'\n" ;

    $self->[SCORE]=shift(@array);
    $self->[STRAND]=shift(@array);    
    $self->[FRAME]=shift(@array);

    my $attributes = shift @array ;
    $self->parseTags($attributes) if $attributes ;

    # Flatten remainder of line beyond any additional tabs as a comment
    my $remark = join ' ', @array ; 
    $remark =~ s/^\s+// ;  # trim leading blanks... 
    $remark =~ s/^#\s*// ; # ... and hash
    if($remark) { 
        # check for a set comment as a precaution:
        # parseTags() may have already found a '#' comment?
        if(my $comment=$self->comment()) { 
	    $self->comment("$comment $remark") ; # append extra remark
        } else {
	    $self->comment($remark) ;
        }
    }
    return $self;
}

=pod

=over 4

=item new_from_parse( $line, \&parser, $version )

Class method to parse the $line string into the
GFF::GeneFeature object using a user-defined
&parser function (creates and returns the object
reference). The &parser should expect the (empty)
GFF::GeneFeature object reference as its first
argument and the input line (string) as its
second argument. So given, the function should
perform the appropriate parsing of the input
$line to load the GFF::GeneFeature object with
data. A typical use for this method is to parse
non-GFF native formats into GFF format.

=back

=cut

# for reading a line parsed by a function
sub new_from_parse {
    my $ref      = shift;
    my $string   = shift;
    my $function = shift;
    my $version  = shift ;
    my $self=$ref->new($version) ;
    if($function){
	if(&$function($self,$string)){
	    return $self;
	}else{
	    return '';
	}
    } else {
	confess "GFF::GeneFeature::new_from_parse(): parse function argument required!";
    }
}

=pod

=over 4

=item parseTags( $attributes )

Method to parse the $attributes string into the
invoking GFF::GeneFeature object [attribute] field.
The method assumes that the invoking object
already knows what GFF version it is.

For GFF Version 1, the entire $attributes
field is taken as the value of the [group] field.

=back

=cut

sub parseTags {
    my $self = shift ;
    my $attributes = shift ;

    confess "GFF::GeneFeature::parseTags(): undefined $attributes?" if !defined($attributes) ;

    if( $self->version() == 1 ) {
	$self->Sequence($attributes) ; # GFF V.1 GROUP == GFF V.2 Sequence field
    } else { 
#      Version 2 or greater: parse $attributes for semicolon ';' delimited 
#      tag-value structures, taking special care about double quoted string context
#      Code adapted from AceParse.pm, courtesy of James Gilbert
       my $tag = '' ;
       my $valuelist = [] ;
       my $end_tag = 0 ;
       my $G_String = $attributes ;
       while (length($G_String)) {
	  my $q = 0; # Set inside quoted text
	  my $field; # Element to add to tags array
	  for (;;) {
	    # Are inside single " ?
	    if ( $q ) {
                # Are at end of quoted text
                if ($G_String =~ s/^"+(\s*|;|\Z)//) {
                    $q = 0;
		    if(!length($G_String) or $1 eq ";") {
			$end_tag = 1 ; # EOL simulated end of tag-value
		    }
                    last;
                # Save " from \" at start of string
		} elsif ($G_String =~ s/^\\(")//) {  
		    $field .= $1;
                # Or everything until next " (could be empty string?)
                } elsif (length($G_String) and $G_String =~ s/^([^"]*)//) {
                    $field .= $1;
                # Or something is wrong
                } else {
	            carp "GFF::GeneFeature::parseTags($self): unbalanced \" in tag field: '$attributes'\n" ;
                    return ;
                }
	    } else {
		# Are at start of new quoted string
		if ( $G_String =~ s/^"// ) {
		    $q = 1;
                # Or see a tag-value delimiter... clear the tag?
		} elsif ($G_String =~ s/^;//) {
		    $end_tag = 1; 
		    $field = undef ;
		    last;
		# or have the start of a '#' 
		} elsif ( $G_String =~ s/^#\s*(.*)// ) {
		    my $comment = $1 ;
		    $self->comment($comment) ; 
		    $end_tag = 1 ; # EOL simulated end of tag-value
		    last ;
                # Or have an unquoted string
		} elsif ($G_String =~ s/^([^\;\s]+)\s*//) {
                    $field = $1;
		    if(!length($G_String)) {
			$end_tag = 1 ; # EOL simulated end of tag-value
		    }
		    last;
                # Or have spaces on start of line
		} else {
                    $G_String =~ s/^\s+//;
                    # Exit infinite loop if nothing left
	            last unless length($G_String);
		}
	     }
	 } # end for(;;)

	 if (defined $field) { 
            if( $tag ) {
               push(@{$valuelist}, $field) ;
            } else {
               $tag = $field ;
            }
	    $field = undef ;
         }
         if($end_tag) {
	     if( $tag ) {
		 if($tag eq 'Sequence') {
		     $self->Sequence($valuelist->[0]) ; # simple label only assumed
		 } else {
		     $self->tagValueList($tag,$valuelist) ;
		 }
	     }
	     $end_tag = 0 ; 
	     $tag = '' ;
	     $valuelist = [] ;
	 }
      } # end while(length($G_String))
   } # else Version 2
}

=pod

=over 4

=item copy( $version )

Method to duplicate the invoking GFF::GeneFeature
object. 

If the optional '$version' argument is
specified (and greater than 0) then the new copy
is cast into the specified version.  This allows
for GFF version casting of GeneFeatureSets.
Note, however, that the GFF Version 1 [group] field
is taken to be equivalent to the GFF Version 2
'Sequence' [attribute] tag-value.

For the 'matches' and member' data,
this just copies 2nd order object references 
(not replicating the objects) from
the old object to the new.

=back

=cut

sub copy {
    my $self    = shift ;
    my $version = shift ;  

    $version = $self->version() if !(defined($version) and $version > 0) ;
    my $gfcopy = GFF::GeneFeature->new($version) ;

    $gfcopy->[SEQNAME] = $self->[SEQNAME];
    $gfcopy->[SOURCE]  = $self->[SOURCE];
    $gfcopy->[FEATURE] = $self->[FEATURE];
    $gfcopy->[GFSTART] = $self->[GFSTART];
    $gfcopy->[GFEND]   = $self->[GFEND];
    $gfcopy->[SCORE]   = $self->[SCORE];
    $gfcopy->[STRAND]  = $self->[STRAND];
    $gfcopy->[FRAME]   = $self->[FRAME];
    $gfcopy->[COMMENT] = $self->[COMMENT] if(defined($self->[COMMENT]));

    if(defined($self->[MATCHES])) {
	$gfcopy->[MATCHES]={} ;
	%{$gfcopy->[MATCHES]} = ( %{$self->[MATCHES]} ) ; 
    }
    if(defined($self->[MEMBER])) {
	$gfcopy->[MEMBER]={} ;
	%{$gfcopy->[MEMBER]}  = ( %{$self->[MEMBER]} ) ;
    } 
    if( $self->version() == 1 and $version == 1) {
        $gfcopy->[GROUP] = $self->[GROUP] if(defined($self->[GROUP])) ;
    } elsif( $self->version() == 1 and $version > 1 ) {
        $gfcopy->[SEQUENCE] = $self->[GROUP] if(defined($self->[GROUP])) ;
    } elsif( $self->version() > 1 and $version == 1 ) {
        $gfcopy->[GROUP] = $self->[SEQUENCE] if(defined($self->[SEQUENCE])) ;
    } else { 
        # Version 2 or greater copied to Version 2 or greater object, 
        # has tag-value matches to copy over (including the 'Sequence' name)?
	$gfcopy->[SEQUENCE] = $self->[SEQUENCE] if(defined($self->[SEQUENCE])) ;
        my $tag ;
	if(defined($self->[TAGS])) {
	    $gfcopy->[TAGS]={} ;
	    foreach $tag (keys %{$self->[TAGS]}) {
		my $valuelist = [] ;
		# $self->[TAGS]->{$tag} should be a reference to an anonymous array
		push @{$valuelist}, @{$self->[TAGS]->{$tag}} if ref($self->[TAGS]->{$tag}) ;
		$gfcopy->[TAGS]->{$tag} = \@{$valuelist} ;
	    }
	}
    }
    return $gfcopy;
}

=pod

=over 4

=item fromHomolGeneFeature( $version )

A "type cast" method to convert an 
HomolGeneFeature object into a GFF::GeneFeature 
object. A new GFF version number may also 
be specified. All optional ([group/attribute], 
comments) and HomolGeneFeature specific data is 
undefined in the new object.

=back

=cut

sub fromHomolGeneFeature {
    my $self    = shift ;
    my $version = shift ;
    my @gf=splice(@{$self},VERSION,FRAME+1) ; # GROUP+ data discarded
    $self=\@gf ;
    bless $self,'GFF::GeneFeature';
    if( defined $version and $version and $version =~/\d+/ ) {
        $self->[VERSION]=$version ;
    }
    return $self ;
}

=pod

=head1 GFF::GeneFeature Output Methods

=over 4

=item dump( \*OUTPUT, $tab, $newline, $flen, $tag )

Method uses dump_string() to write a formatted
output of a GFF::GeneFeature object to a
filehandle, OUTPUT. If \*OUTPUT is not given,
\*STDOUT is used. When the optional $tab, $flen
and $newlines arguments are omitted, this method
is guaranteed to dump well-formed GFF records
meeting GFF version standards. The use of the
optional parameters provide alternate, non-GFF
compliant tabular text formats for output. 

=over 4

=item C<tab>

The "$tab" argument is a boolean flag, where a "true"
(non-null) value directs the use of tab as the
field delimiter in the output line; otherwise,
use blank space (flag is assumed "true"
(non-null) if not specified). 

=item C<newline>

The "$newline" argument is passed to dumpTags 
(see below) via dump_string. 

=item C<newline>

The "$flen" argument is a boolean 
flag (assumed "false" if not specified), where a
non-null value stipulates that the length of the
current output line should be printed as an extra
field at the end of the output line. (Note: the
extra length of this field is *not* added to the
displayed line size, but the extra field is tab
delimited, if $tab is set).

=item C<tag>

Version 2 GFF: the optional $tag argument, passed
to the dump_string() method, restricts [attributes] field
dumping of the GFF to the specified tag. If $tag is
undefined, all [attributes] tag values are dumped. If
$tag is defined but null (empty string), then
no [attributes] tag-values are dumped. Otherwise, a
defined and non-empty $tag value is assumed to be
a simple identifier or Perl regex matching tags
whose values are to be included in the dump.

=back

=cut

sub dump {
    my $self    = shift ;
    my $file    = shift ;
    my $tab     = shift ; 
    my $newline = shift ;
    my $flen    = shift ;
    my $tag     = shift ;

    my $len;
    
    if( !(defined $file and $file) ){
	$file = \*STDOUT;
    }
    print $file $self->dump_string($tab,$newline,$tag);
    if($flen){
	$len=$self->length;
	if(!defined $tab or $tab){
	    print $file "\t$len";
	}else{
	    print $file " $len";
	}
    }
    print $file "\n";
}

=pod

=over 4

=item dump_string( $tab, $newline, $tag )

Method to return a formatted output of all the
fields (including [attributes]) of a GFF::GeneFeature
object to a string. 

The optional "$tab" argument
is a boolean flag, where a "true" (non-null)
value directs the use of tab as the field
delimiter in the output line; otherwise, use
blank space (flag is assumed "true" (non-null) if
not specified). 

The "$newline" and $tag arguments
are passed to dumpTags (see below). See dump()
method above and the dumpTags() method below 
for explanation about $tag argument.

=back

=cut

sub dump_string {
    my $self    = shift;
    my $tab     = shift ;
    my $newline = shift ;
    my $tag     = shift ;
    $tab = 1 if !defined $tab ; # $tab is assumed true if unspecified (undef)
    my $comment = $self->comment() ;
    if( defined($comment) and $comment ) {
	$comment = " # ".$comment ;  # '#' used to delimit comment...
    } else {
	$comment = '' ;
    }
    my $string ;
    if($tab){
	$string = sprintf("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t",
		$self->seqname(),
		$self->source(),
		$self->feature(),
		$self->[GFSTART],
		$self->[GFEND],
		$self->[SCORE], # version 2 GFF can have '.'
		$self->[STRAND],
		$self->[FRAME]) ;
    }else{
	$string = sprintf("%s %s %s %d %d %s %s %s ",
		$self->seqname(),
		$self->source(),
		$self->feature(),
		$self->[GFSTART],
		$self->[GFEND],
		$self->[SCORE], # version 2 GFF can have '.'
		$self->[STRAND],
		$self->[FRAME]) ;
    }
    if(defined($tag) && !$tag) { 
	$string .= $comment ; # completely suppress [attributes] field dumping?
    } else {
	$string .= ($self->dumpTags($tab, $newline, length($string), $tag).$comment) ;
    }
}

=pod

=over 4

=item dumpTags( $tab, $newline, $len, $tag )

Method to return a formatted output of a
GFF::GeneFeature object [group/attributes] field to a
string. 

For GFF Version 1, the string of the whole [group] 
field is returned.

For GFF Version 2 object dumps, tag-value sets consisting
of semicolon delimited tag-value structures (see
GFF specification) are returned.

The optional "$tab" argument
is a boolean flag, where a "true" (non-null)
value directs the use of tab as the field
delimiter in the output line; otherwise, use
blank space (flag is assumed "true" (non-null) if
not specified).

The "$newline" argument is a boolean flag
(assumed "false" if omitted), where a "true"
(non-null) value directs both that [attributes] field
tag-value pairs are printed one per line (with
the semi-colon omitted), and that any '\n'
characters in normally double quoted free text
[attributes] value strings are converted to "real"
newlines which wrap the text into multiline
printed format. Both newline effects are printed
within the restricted context of the [attributes]
field column. In other words, the $newline flag
is used for some semblance of "pretty printing"
TAG fields.

If $tab and/or $newline are specified, then 
the $len argument should contain the length of the
dump line preceeding the [attributes] field.

If the $tag argument is undefined, then all 
tag-value pairs are dumped. Otherwise, $tag is 
assumed to be a simple tag or Perl regex expression
matching tag value fields which are to be dumped.

Note: tags are dumped in Perl unordered (by 'each')

=back

=cut

sub dumpTags {
    my $self    = shift ;
    my $tab     = shift ;
    my $newline = shift ; # use '\n' to split tag-value matches and embedded newlines; otherwise GFF-like
    my $len     = shift ; # size of line before the group/attribute field

    # optional: controls what tag fields are dumped, 
    # undef => dump all tags, otherwise, 
    # specified tag-value field(s) only are dumped 
    my $theTag = shift ; 

    $tab = 1     if !defined $tab ;
    $newline = 0 if !defined $newline ;
    $len = 64    if !defined $len or !$len ; # defaults to 8x8 sp/tab width

    if( $self->version() == 1 ) {
        return $self->Sequence() ; # GFF V.1 GROUP == GFF V.2 SEQUENCE
    } else { # GFF V.2+
        my $string = '' ;

	# maybe a 'Sequence' tag-value?
	$string .= ('Sequence "'.$self->Sequence().'"') 
	    if(defined($self->[SEQUENCE])) ;

	# other tag-value?
	return $string if !defined($self->[TAGS]) ; # null string if empty
	confess "GFF::GeneFeature::dumpTags() invalid tag-value hash?" 
	    if(ref($self->[TAGS])!~/HASH/) ;

        # Dump all other attribute tag-values,
        # sorted by tag; not necessarily the original order :-)
        while( my ($tag,$values) = each(%{$self->[TAGS]}) ) {
	    next if (defined($theTag) and ($tag !~ /^$theTag$/)) ;

	    # offset to [attributes] for next line of tag-value
            my $offstr ;
            if($tab) {
	       $offstr = "\n".(" "x72) ; # estimate only: 8 regular tabs + one overflow field, @ 8 sp/tab
	    } else {
	       $offstr = "\n".(" "x($len+length($tag)+3)) ;
	    }

	    unless (! $string ) { # except for the first tag-value pair
	       if( !$newline ) {
		   $string .= "; " ; # delimit off the previous tag-value pair first
	       } else {
		   # offset to [attributes] for next line of tag-value
		   $string .= $offstr ;
	       }
	    }
            $string .= "$tag " ;
	    if(defined($values) && ref($values)=~/ARRAY/) {
		foreach my $avalue (@{$values}) {
		    if( !defined($avalue) ) {
			carp "### Warning: Undefined value in tagValueList($tag) encountered?\n" ;
			next ;
		    }
		    if(!$avalue) { # NULL values here will be interpreted as numeric zeros
			$string .= " 0" ; 
		    } elsif( $avalue =~ /^\s*(\+|\-)?\d+(\.\d+)?$/ ) { # if numeric
			$string .= " $avalue" ; # then leave unquoted?
		    } else {
			if( !$newline ) {
			    $string .= " \"$avalue\"" ; # simply return a double quoted text value
			} else {
			    # convert embedded newlines in $value to a real newlines+TAG field offset
			    $avalue =~ s|\\n|$offstr|g ;
			    $string .= "  \"avalue\"" ;
			}
		    }
		}
	    }
	    $string =~ s/\s+$// ; # trim any trailing blanks, just in case...
	}
        return $string ;
    } # else GFF V.2+
}

=pod

=over 4

=item dump_matches( \*OUTPUT, $tab )

Method uses dump_string() to write a formatted
output of a GFF::GeneFeature object to a
filehandle. Includes information about any
(overlap) matching GFF::GeneFeatures (note: match
output for each record is multi-line, the matches
designated by an indented '=>' bullet). If
\*OUTPUT is not given, \*STDOUT is used. The
"$tab" argument is a boolean flag, where a "true"
(non-null) value directs the use tab as the field
delimiter in the output line; otherwise, blank
space is used as the delimiter (assumed "true"
(non-null) if not specified).

=back

=cut

sub dump_matches {
    my $self = shift;
    my $file = shift;
    my $tab = shift;
    
    if( !defined $file ){
	$file = \*STDOUT;
    }
#
#   Print self...
#   Only dump feature w/ 'Sequence' tag-value if seen
#
    print $file $self->dump_string($tab,0,'Sequence') ; 
#
#   ... then get references to objects matched with self
#
    if( $self->getMatches_logical() ){
        print $file ": matches\n" ;
#
#       Retrieve reference to a hash keyed on GeneFeature references
#
        my $match = $self->getMatches() ; 
        foreach my $gf (keys %{$match}) {
#           
#           $tab is assumed "true" if omitted (undef)
#
            my $out;
	    if(!defined $tab or $tab) {
	        $out = "\t";
	    } else {
	        $out = " ";
	    }
	    $out .= sprintf("Match Offsets: 5'(%d), 3'(%d)",
			       $match->{$gf}->[1], # e5'
			       $match->{$gf}->[2], # e3'
			   );

            print $file "=>\t",$match->{$gf}->[0]->dump_string($tab),"$out\n" ;
        }
    } else {
        print $file ": no matches\n" ;
    }


}

=pod

=head1 GFF::GeneFeature Access Methods

The various GFF record fields may be set or
queried by the following access methods. All the
methods can take a single string argument to set
the variable. With or without an argument, the
methods return the current (or newly set) value,
as a string, except as specifically noted below:

=over 4

=item version()

GFF version.

=back

=cut

sub version {
    my $self = shift;
    my $version = shift ;
    my $class = ref($self) || carp "\'$self\' is not an object of mine!";
    
    if ($version) {
	return $self->[VERSION] = $version+0;
    } else {
	return $self->[VERSION];
    }
}

=pod

=over 8

=item seqname()

Name of the host sequence.

=back

=cut

sub seqname {
    my $self  = shift;
    my $class = ref($self) || carp "\'$self\' is not an object of mine!";
    my $label = shift ;
    
    if(defined($label)) {
	return ${$self->[SEQNAME] = uid($label)};
    } elsif(defined($self->[SEQNAME]) and 
	    ref($self->[SEQNAME])) {
	return ${$self->[SEQNAME]};
    } else {
        return undef ;
    }
}

=pod

=over 8

=item source()

Source of the sequence.

=back

=cut

sub source {
    my $self = shift;
    my $class = ref($self) || carp "\'$self\' is not an object of mine!";
    my $label = shift ;
    
    if (defined($label)) {
	return ${$self->[SOURCE] = uid($label)};
    } elsif(defined($self->[SOURCE]) and 
	    ref($self->[SOURCE])) {
	return ${$self->[SOURCE]};
    } else {
        return undef ;
    }
}

=pod

=over 8

=item feature()
eneFeature
Feature type name.

=back

=cut

sub feature {
    my $self = shift;
    my $class = ref($self) || carp "\'$self\' is not an object of mine!";
    my $label = shift ;
    
    if (defined($label)) {
	return ${$self->[FEATURE] = uid($label)};
    } elsif(defined($self->[FEATURE]) and 
	    ref($self->[FEATURE])) {
	return ${$self->[FEATURE]};
    } else {
        return undef ;
    }
}

=pod

=over 8

=item start()

Start coordinate of feature.

=back

=cut

sub start {
    my $self = shift;
    my $class = ref($self) || carp "\'$self\' is not an object of mine!";
    my $start = shift ;
    if (defined($start)) {
	if($start !~ /^[\+-]?\d+$/) {
	    carp "GFF::GeneFeature::start must be an integer\n" ;
	    return undef ;
	} else {
	    return $self->[GFSTART] = $start;
	}
    } else {
	return $self->[GFSTART];
    }
}

=pod

=over 8

=item end()

End coordinate of feature.

=back

=cut

sub end {
    my $self = shift;
    my $class = ref($self) || carp "\'$self\' is not an object of mine!";
    my $end = shift ;
    if (defined($end)) {
	if($end !~ /^[\+-]?\d+$/) {
	    carp "GFF::GeneFeature::end must be an integer\n" ;
	    return undef ;
	} else {
	    return $self->[GFEND] = $end;
	}
    } else {
	return $self->[GFEND];
    }
}

=pod

=over 8

=item score()

Source score of feature (by method). Returned as a string even if a float number.

=back

=cut

sub score {
    my $self = shift;
    my $class = ref($self) || carp "\'$self\' is not an object of mine!";
    my $score=shift ;
    if (defined($score)) {
	if($score !~ /^[\+-]?(\d+(\.\d+)?|\.)$/) {
	    carp "GFF::GeneFeature::score must be a integer/float number or '.'\n" ;
	    return undef ;
	} else {
	    return $self->[SCORE] = $score;
	}
    } else {
	return $self->[SCORE];
    }
}

=pod

=over 8

=item strand()

"+" (forward), "-" (reverse) or "." (n/a).

=back

=cut

sub strand {
    my $self = shift;
    my $class = ref($self) || carp "\'$self\' is not an object of mine!";
    my $strand = shift ;
    if(defined($strand)) {
	if($strand !~ /^[\+-\.]$/) {
	    carp "GFF::GeneFeature::frame must be one of '+','-' or '.'\n" ;
	    return undef ;
	} else {
	    return $self->[STRAND] = $strand;
	}
    } else {
	return $self->[STRAND];
    }
}

=pod

=over 8

=item frame()

"0", "1", "2" or "." (n/a)

=back

=cut

sub frame {
    my $self  = shift;
    my $class = ref($self) || carp "\'$self\' is not an object of mine!";
    my $frame = shift ;
    if (defined($frame)) {
	if($frame !~ /^[012\.]$/) {
	    carp "GFF::GeneFeature::frame must be one of '0','1','2' or '.'\n" ;
	    return undef ;
	} else {
	    return $self->[FRAME] = $frame;
	}
    } else {
	return $self->[FRAME] ;
    }
}


# Subsumed into GFF::GeneFeature by GFF Version 2
sub start2 {
    my $self = shift;
    my $class = ref($self) || carp "\'$self\' is not an object of mine!";
    my $start2 = shift ;
    if (defined($start2)) {
	if($start2 !~ /^[\+-]?\d+$/) {
	    carp "GFF::GeneFeature::start2 must be an integer\n" ;
	    return undef ;
	} else {
	    return $self->group_value('Target',1,$start2);
	}
    } else { 
	return $self->group_value('Target',1);
    }
}

sub end2 {
    my $self = shift;
    my $class = ref($self) || carp "\'$self\' is not an object of mine!";
    my $end2 = shift ;
    if (defined($end2)) {
	if($end2 !~ /^[\+-]?\d+$/) {
	    carp "GFF::GeneFeature::end2 must be an integer\n" ;
	    return undef ;
	} else {
	    return $self->group_value('Target',2,$end2);
	}
    } else {
	return $self->group_value('Target',2);
    }
}

sub percentid {
    my $self = shift;
    my $class = ref($self) || carp "\'$self\' is not an object of mine!";
    my $pcid=shift ;
    if (defined($pcid)) {
	if($pcid !~ /^[\+-]?(\d+(\.\d+)?|\.)$/) {
	    carp "GFF::GeneFeature::percentid must be a integer/float number or '.'\n" ;
	    return undef ;
	} else {
	    return $self->group_value('Percent_Id',0,$pcid);
	}
    } else {
	return $self->group_value('Percent_Id');
    }
}

sub exp_value {
    my $self = shift;
    my $class = ref($self) || carp "\'$self\' is not an object of mine!";
    my $evalue=shift ;
    if (defined($evalue)) {
	if($evalue !~ /^[\+-]?(\d+(\.\d+)?|\.)$/) {
	    carp "GFF::GeneFeature::percentid must be a integer/float number or '.'\n" ;
	    return undef ;
	} else {
	    return $self->group_value('E_Value',0,$evalue);
	}
    } else {
	return $self->group_value('E_Value');
    }
}

=pod

=over 8

=item group($value) - GFF Version 1

Returns or sets the optional [group] field 
$value, assumed to be a scalar string.

=item group($tag,@values) - GFF Version 2

Calls the attribute() method below, unless the
given $tag is equal to 'Sequence', which is treated
differently, in that it is assumed to be a 
single valued string, hence is returned this way.
(In fact, use of the group() method in this context
is deprecated: use the 'Sequence' method instead).

=item attribute($tag,@values) - Version 2+ GFF

Under the Version 2 specification, the 
(optional) [attribute] field of a GFF record 
must be structured as an ACEDB B<.ace> style 
B<tag-value> set, flattened to one line by 
using semicolon delimiters instead of newlines. 

The attribute() method can set or return these
tag-values:

- If no arguments are given, then the method
  returns the (possibly undefined if empty)
  reference to the anonymous Perl hash
  which indexes references anonymous arrays 
  of tag-values indexed by the tag names.

- If the $tag argument only is given, then
  the associated (possibly undefined) 
  anonymous array of tag-values is returned.
  (Note: valueless tags are created in the 
   in the GeneFeature object as a side effect 
   of the attribute call, but return 'undef')
 
- If a $tag and @values are provided in the 
  call, the $tag is set to be an anonymous
  array containing the @values (overwriting 
  any previous values - see also tagValueList()).

Note: to simply test the existence of a $tag
in a given GeneFeature record, use the $tag
name as a method invocation (which is AUTOLOAD
tested in the [attribute] tag array).

=back

=cut

# GFF Version 2+ only
sub attribute {
    my $self = shift;
    my $class = ref($self) || carp "\'$self\' is not an object of mine!";
    if(my $tag=shift) {
	if($tag eq 'Sequence') { # V.2+
            # single scalar string 'Sequence' tag-value assumed
	    return $self->Sequence(shift) ; 
	} else {
	    $self->[TAGS]={} if !(defined($self->[TAGS])) ;
	    confess "GFF::GeneFeatureSet::attribute(): invalid tag-value hash?!" 
		if(!ref($self->[TAGS])=~/HASH/) ;
	    my @values = @_ ;
	    if(scalar(@values)) {
		return ($self->[TAGS]->{$tag} = \@values) ;
	    } else {
		unless(exists($self->[TAGS]->{$tag})) {
		    return ($self->[TAGS]->{$tag} = undef) ; # but tag now exists...
		} else {
		    return ($self->[TAGS]->{$tag}) ;
		}
	    }
	}
    } else {
	return($self->[TAGS]) ;
    }
}

sub group {
    my $self = shift;
    my $class = ref($self) || carp "\'$self\' is not an object of mine!";
    
    if( $self->version() == 1 ) {
	$self->Sequence(shift) ; # V.1 GROUP value == V.2 SEQUENCE attribute
    } else {
	$self->attribute(@_) ;
    }
}

=pod

=over 8

=item group_value_list

Same as tagValueList() -- deprecated, use tagValueList

=item tagValueList( $tag, \@values, $append ):

Version 2 GFF: method sets and/or returns the
reference to the array of values associated with
a given $tag of a (Version 2) tag-value
structure. The (optional) specification of a
reference to such an array of values (C<\@values>)
resets the $tag hash value to the new array
reference (and returns it).

Caution: modifying the values of the array using the
returned reference, *does* modify the specific 
tag values of the given gene feature! Make a copy
of the array if that is not your intent!

If $tag is undefined or null, then the function 
just returns the value of the 'Sequence' tag as
as a single valued list reference
(note: this may be empty). 

If the $append argument is defined and non-NULL, 
then the given @values are appended to any existing 
value list (default: 0) 

Again, the 'Sequence' tag is a special case, 
for which only the first value in the $valuelist 
is kept, so the '$append' argument is ignored here. 

Version 1 GFF: just returns the group string
value embedded in a single member list. 
All the other arguments are ignored.

=back

=cut

sub group_value_list {
    tagValueList(@_) ;
}

sub tagValueList {
    my $self      = shift ;
    my $tag       = shift ; # a string tag name
    my $valuelist = shift ; # (optional) Version 2 $values should be a reference to an array of values?
    my $append    = shift ; # append these values rather than overwrite given tag values

    confess "*** Error: GFF::GeneFeature::tagValueList(): a defined valuelist argument *must* be an ARRAY reference\n"
	if(defined($valuelist) and $valuelist and ref($valuelist) !~ /ARRAY/) ;

    if( $self->version() == 1 ) {
	if(defined($self->Sequence)) {
	    return [$self->Sequence] ; # GFF V.1 GROUP == V.2 SEQUENCE
	} else {
	    return undef ;
	}
    } else { # Version 2+ GFF
         if(!(defined($tag) and $tag)) {
	     return [$self->Sequence] ;
         } else {
	    if($tag eq 'Sequence') {
		my $value ;
	        if( defined ($valuelist) and $valuelist) { 
		    return [$self->Sequence($valuelist->[0])] ;
		} else {
		    if(defined($self->Sequence)) {
			return [$self->Sequence] ; # GFF V.1 GROUP == V.2 SEQUENCE
		    } else {
			return undef ;
		    }
		}
            } else { # general tag
		my $attributes ;
	        if( defined ($valuelist) and ref($valuelist) =~ /ARRAY/ ) { # values to set?
		    $self->[TAGS]={} if !(defined($self->[TAGS])) ; # incoming new values? allocate hash if still undefined
		    $attributes = $self->[TAGS] ;
		    if(defined($attributes->{$tag}) && $append) {
		        push(@{$attributes->{$tag}}, @{$valuelist}) ; # append the new array of values to any old
		    } else {
		        $attributes->{$tag} = $valuelist ;   # overwrite with reference to the new array of values
		    }
                } else {
		    $attributes = $self->[TAGS] ; 
		    return undef if !(defined($attributes) and defined($attributes->{$tag})) ; # no attributes to report?
		}
		return $attributes->{$tag} ;
            }
        }
    }
}

=pod

=over 8

=item group_value( $tag, $index, $value )

Same as tagValue() -- deprecated (use tagValue() directly)

=item tagValue( $tag, $index, $value )

(Version 2 GFF) method sets and/or returns the
element at the (zero based)"index" position in
the value list associated with a given $tag of a
(Version 2) tag-value structure. 

If $tag is undefined, then the value of the
'Sequence' tag value (if any) is returned.

If an $index is not given, the first value of 
the value list is returned. 

If $value is specified, the ith element
is set to it (note: if the $tag associated value
array is undefined when this method is called,
then it is created and its value is set to a
single element list containing $value). 

(Version 1 GFF) just returns the [group/Sequence] string, 
ignoring any $tag or $index provided. If $value is
provided, then the [attributes] name is reset.

=back

=cut

sub group_value {
    tagValue(@_) ;
}
    
sub tagValue {
    my $self  = shift ;
    my $tag   = shift ;
    my $index = shift ; 
    my $value = shift ; # optional; generally a scalar, but could be a list in version 2

    if( $self->version() == 1 or
	!defined($tag) or
        $tag eq 'Sequence' ) {
        return $self->Sequence($value) ;
    } else {
        # default to 0th value, if index is not given
        $index = 0 if !defined $index ;

        my $valuelist = $self->tagValueList($tag) ; # retrieves an array reference?

        if( defined($value) ) {  # resetting the value?
            unless( $valuelist && ref ($valuelist) =~ /ARRAY/ ) { # is the value list array defined?
		my $new_values = [] ;
		$new_values->[$index] = $value ; 
                $self->tagValueList($tag,$new_values) ; 
	        return $value ;
	    }
            return $valuelist->[$index] = $value ; # set and return the ith element of the list (zero based)
	} else {
            return undef unless ref($valuelist);
            return $valuelist->[$index]; # return the ith element of the list (zero based)
        } 
   }
}

=pod

=over 8

=item TAG(\@values) - Version 2 GFF only

[TAGS] field tags given as method names are now
Perl AUTOLOAD recognized as accesses to the
tag-value field. An optional reference to an
array of tag-values may be given to associate 
with the tag.

With or without values, the current
value list array reference is returned. 

Valueless tags simply return a boolean '1' 
if they exist. If the tag does NOT exist 
in the given object, then '0' is returned.

=item Sequence( $value )

A special 'tag-value' case, corresponding 
to the 'Sequence' tag. Usually a single valued
tag-value list (and stored this way).
Use directly instead of tagValue().  
In Version 1 GFF, synonymous with group().

=back

=cut

sub AUTOLOAD {
    my $self = shift;
    my $class = ref($self) || carp "\'$self\' is not an object of mine!";
    my $name = $AUTOLOAD;

    # don't propagate DESTROY messages...
    $name =~ /::DESTROY/ && return;
    $name =~ s/.*://; #get only the bit we want

    if( $self->version() >= 2 && 
	defined($self->[TAGS])) {
	if( exists $self->[TAGS]->{$name} ) { 
	    if (@_) {
		return $self->[TAGS]->{$name} = shift ; # should ideally be a reference to an array...
	    } elsif(defined($self->[TAGS]->{$name})) {
		return $self->[TAGS]->{$name} ;
	    } else {
		return 1 ; # valueless but tag now exists, so return a boolean
	    }
	} else {
	    return 0 ; # silently fails...
	}
    } else { # version 1 is simplistic
        confess "In type $class, can't access $name - probably passed a wrong variable into $class\n";
    }
}

=pod

=over 8

=item Sequence( $label )

A special 'tag-value' case, corresponding 
to the 'Sequence' tag, which is usually a single 
(string scalar) valued tag-value list 
(and is stored as such).
Use directly instead of tagValue(). 

In Version 1 GFF, synonymous with group().

=back

=cut

sub Sequence {
    my $self  = shift;
    my $class = ref($self) || carp "\'$self\' is not an object of mine!";
    my $label = shift ; 

    if(defined($label)) {
	return ${$self->[SEQUENCE] = uid($label)};
    } elsif(defined($self->[SEQUENCE])) {
	my $ref=ref($self->[SEQUENCE]) ;
	if($ref=~/ARRAY/) {
	    return $self->[SEQUENCE]->[0];
        } else {
	    return ${$self->[SEQUENCE]};
	}
    } else {
        return undef ;
    }
}

=pod

=over 8

=item grepTag($tpat,$vpat)

This boolean predicate method tests the [attribute] field
for a tag matching a specified Perl regex pattern '$tpat'.
If the optional $vpat pattern is also provided, then
this is applied to the value lists of the $tpat hits.
If $tpat is undefined, then the $vpat is applied to all
values in the [TAGS] tag-values.

Returns 0 if unsuccessful; Otherwise, the set of hits 
is returned in an array context, or the count of hits
in a scalar context. If only the $tpat pattern is provided, then
the set of hits consist of a list of references to
matching tag-value arrays. If the $vpat pattern is also
provided, then only the matching values, of the tag-value
pairs, are returned. 

Version 1 GFF - Just applies the $tpat pattern to 
the whole [TAGS] field. The $vpat is ignored.

Version 2 GFF - Applies the pattern(s) to all
tag-values of the [TAGS] field.

=back

=cut

sub grepTag {
    my $self = shift;
    GFF::GeneFeature->verify($self);
    my $tpat = shift ;

    my $context = wantarray ;
    if( $self->version() == 1 ) {
	if( !defined($tpat) or ($self->[GROUP] =~ /$tpat/)) {
	    my $hit = $& ;
	    if(!defined $context) {
		return ;
	    }elsif($context){
		return ($hit); # array context
	    }else{
		1 ; # scalar context
	    }	    
	} else {
	    return 0 ;
	}
    } else { # version 2+ tag-values
        my @hits ;
	my $vpat = shift ;
	foreach my $tag (keys %{$self->[TAGS]}) {
	    if( !defined($tpat) or ($tag =~ /$tpat/)) {
		my $values =$self->[TAGS]->{$tag} ; 
		if(defined($vpat)) {
		    if(defined($values)) {
			my @vhits = grep $vpat, @{$values} ;
			push(@hits,@vhits) ;
		    }
		} else {
		    if(defined($values)) {
			push(@hits,[$tag,@{$self->[TAGS]->{$tag}}]) ;
		    } else {
			push(@hits,[$tag]) ;
		    }
		}
	    }
	}
	if(!defined $context) {
	    return ;
	}elsif($context){
	    return @hits ; # array context
	}else{
	    return scalar(@hits) ; # scalar context
	}
    }
}

=pod

=over 8

=item renameTag($oldtag,$newtag)

This (Version 2 GFF) method renames [attribute] tags.
Note: this operation directly modifies the 
gene feature objects concerned.

=back

=cut

sub renameTag {
    my $self = shift;
    GFF::GeneFeature->verify($self);

    if( $self->version() >=2 ) {
	my $oldtag = shift ;
	my $newtag = shift ;
	confess "GFF::GeneFeature::renameTag(): both old and new tag arguments must be defined\n" 
	    unless(defined($oldtag) and $oldtag and defined($newtag) and $newtag) ;
	my $values ;
	if(exists($self->[TAGS]->{$oldtag})) {
	    $values = $self->[TAGS]->{$oldtag} ;
	    delete($self->[TAGS]->{$oldtag})
	} else { # no old tag?
	    $values = [] ;
	    carp "GFF::GeneFeature::renameTag() warning: tag '$oldtag' does not exist but empty tag '$newtag' will be created\n" ;
	}
	$self->[TAGS]->{$newtag} = $values ;
    }
}

=pod

=over 8

=item deleteTag($tag)

This method deletes [attribute] values and tags.
Note: this operation directly modifies the 
gene feature objects concerned.

Version 1 GFF - $tag argument is ignored; the method
undefines the [attribute] field value.

Version 2 GFF - if no $tag argument is provided,
the entire [attribute] tag-value array is cleared of
tags and values. If a $tag argument is provided, 
only the indicated tag (if it exists) is deleted
from the [attribute] field, along with any associated
values.

=back

=cut

sub deleteTag {
    my $self = shift;
    GFF::GeneFeature->verify($self);

    if( $self->version() == 1 ) {
	$self->[TAGS]=undef ;
    } else {
	my $tag = shift ;
	if(!defined($tag)) {
	    # clear the whole hash
	    undef($self->[TAGS]) ;
	} else { # defined tag
	    if($tag eq 'Sequence') {
		undef($self->[SEQUENCE]) ;
	    } else {
		delete($self->[TAGS]->{$tag}) ;
	    }
	}
    }
}

=pod

=over 8

=item comment()

trailing line comments associated with a GFF
record. Under Version 2, such comments, starting
after the [attribute] field, must be delimited with a
'#' character.

=back

=cut

sub comment {
    my $self = shift;
    GFF::GeneFeature->verify($self);
    
    if (@_) {
	return $self->[COMMENT] = shift;
    } else {
	return $self->[COMMENT];
    }
}


=pod

=head1 GFF::GeneFeature Analysis Methods

The following methods compute upon GFF features.

=over 4

=item length()

Calculate the segment length (end-start+1) of a
given feature.

=back

=cut

sub length {
    my $self = shift;
    my $len = $self->[GFEND]-$self->[GFSTART]+1;
    return $len ;
}

=pod

=over 4

=item remap( $offset, $rcend )

Method to add an $offset amount to the start and
end coordinates of a GFF::GeneFeature. Note: the
start and end of the original object, not a copy,
are changed.

If the optional $rcend argument is defined 
and non-NULL, then the remapping is presumed to be a 
reverse complemention operation in which the $rcend
point remaps to the $offset+1 coordinate position.
This value is generally the end coordinate of the
'##sequence-region' of the source GFF::GeneFeatureSet.

Note that if reverse complementation is done,
then the sign of the <strand> field, if known,
is toggled.

=back

=cut

sub remap {
    my $self    = shift ;
    my $offset  = shift ;
    my $rcend   = shift ;

    unless($offset or defined($rcend)) {
        return ; # do nothing
    } else {
	my ($start,$end) ;
	unless(defined($rcend)) {
	    $start=$self->[GFSTART];
	    $end  =$self->[GFEND];
	} else { # reverse complementation
	    $end   = $rcend-$self->[GFSTART]+1 ;
	    $start = $rcend-$self->[GFEND]+1 ;
	    if($self->[STRAND] =~ /^([-\+])$/) {
		$self->[STRAND]=(($1 eq '+')?'-':'+') ;# toggle strand
	    }
	}
	$offset=0 if!defined($offset) ;
	$self->[GFSTART]=$offset+$start ;
	$self->[GFEND]  =$offset+$end ;
    }        
}

=pod

=over 4

=item match( $GF2, $tolerance, $single, $strand )

Method to compare two GFF::GeneFeature objects to
look for an overlap returns 5 scalars as an
array. The first GFF::GeneFeature object invokes
the method giving the second GFF::GeneFeature
object as the first argument ("GF2" above), plus
optional method arguments "$tolerance", "$single"
and "$strand" to guide the analysis (each are
assumed to be 0 when not explicitly given).

Based upon their specified start and end
coordinates, two GFF::GeneFeatures will either
overlap perfectly, partially or not at all (are
"misses"). The $tolerance value specified
controls the match decision for each category of
overlap as follows:

=over 4

=item 1

Specifying a $tolerance value of 0 dictates that
an exact match is required, that is, that the
corresponding 5' and 3' coordinate ends of both
GFF::GeneFeature objects must be equal to one
another.

B<Release 2.106 revision>: optionally, the $tolerance 
argument can now be a reference to an array  5' and 3'
end specific tolerance pairs [t55, t53, t35, t33],
where t55 => the second gene feature lies within
t55 basepairs 5' of the 5' end of the first gene feature,
t53 == 3' of the 5' end of the gene feature, etc.
(Note: 5' is a function of strandedness, if any, or
simply 'start' for '.' strand objects),

=item 2

For partial overlaps of GFF::GeneFeature objects,
if the "$tolerance" is set to a negative number
then, ceteris paribus, then two overlapping
GFF::GeneFeature objects are matched
unconditionally.

=item 3

If two GFF::GeneFeature object segments overlap
imperfectly but a positive, non-zero $tolerance.
is specified, then a match is successful if the
(absolute value of the) extent of both the 5' and
3' mismatch in coordinates is less than or equal
to the tolerance value (see also the effect of
the "$single" method argument below).

=item 4

If two GFF::GeneFeature object segments do not
overlap at all, but if a negative $tolerance
value less than -1 is specified, then a match is
declared if the misses are within $tolerance,
that is, if the difference in coordinates of the
closest segment ends of the two GFF::GeneFeature
objects is less than the $tolerance.

=back

If a value of 1 is given for the "$single"
argument to the method, then the above $tolerance
conditions, for positive $tolerance values and
imperfect GFF::GeneFeature overlaps, are relaxed
such that either a 5' or a 3' end mismatch within
tolerance results in a positive match.

If "$strand" is given (i.e. not 0), then the
match fails unless the two GFF::GeneFeature lie
on the same strand. If $strand is zero, then
strandedness of features is completely ignored
in the match comparison (i.e. can also be '.' == unknown)

The 5 scalar fields returned in the match @result
array have the following values:

=over 4

=item $result[0]

1 ("true") or 0 ("false") indicates if objects match by criteria supplied

=item $result[1]

For overlaps: Number indicating signed difference in 5' end of 2 objects
For misses within tolerance: 3' segment end "closest" coordinate of the 5'-most
GFF::GeneFeature

=item $result[2]

For overlaps: Number indicating signed difference in 3' end of 2 objects
For misses within tolerance: 5' segment end "closest" coordinate of the 3'-most
GFF::GeneFeature

=item $result[3]

Detailed match information: 0 if no overlap; 1 if match; 2 if overlap, but rejected on
tolerance+single; 3 if overlap, but rejected on strand; 4 if no overlap but accepted on
tolerance

=item $result[4]

Descriptive string about overlap, if there is one

=back

=back

=cut

sub match {
    my $self      = shift;
    my $other     = shift;
    my $tolerance = shift;
    my $single    = shift;
    my $strand    = shift;

    my($a1s,$a1e,$a2s,$a2e,$str,
       $maxs,$mine,$f3,$f5,
       $match,$string,$type,
       $len,$slen);

    $tolerance = 0 if !defined $tolerance ;
    $single    = 0 if !defined $single ;
    $strand    = 0 if !defined $strand ;
#
# Check if Genefeature objects overlap
# in absolute terms, along the coordinate axis
#
    $a1s  = $self->start();
    $a1e  = $self->end();
    $a2s  = $other->start();
    $a2e  = $other->end();
    $str  = $self->strand();
    $maxs = &max( $a1s, $a2s );
    $mine = &min( $a1e, $a2e );

    if($maxs<=$mine){
#############################################
# Yes! the two features absolutely overlap! #     
#############################################
#
# define string
#
	$len  = $mine - $maxs + 1 ;
	$slen = &min($a1e-$a1s,$a2e-$a2s) + 1 ;

        my $label1 = $self->tagValue() ;
        $label1 = "Unknown" if !defined($label1) ;

        my $label2 = $other->tagValue() ;
        $label2 = "Unknown" if !defined($label2) ;

	$string = $len.'/'.$slen.' '.$label1.": $a1s-$a1e <=> $a2s-$a2e :".$label2;

	GFF::TRACE(2,"GeneFeature::match() => $string\n");
#
#       Reject on strand early, iff $strand argument is set
#       If $strand not set, ignore strand in the comparison (even if '.')
#
	if( $strand and ( 
	    ($str ne $other->strand) or
             $str eq '.' )
           ) {
	     $type = 3;
	     return 0,0,0,$type,$string;
	}
#	
#       Absolute offset of $other GeneFeature coordinates relative to $self
#       (i.e. positive means that second one lies downstream of first)
#
	my $fw = ($str eq '+' or $str eq '.') ;
	if($fw){ 
	    $f5 = $a2s - $a1s ;
	    $f3 = $a2e - $a1e ;
	} else {
	    $f3 = $a2s - $a1s ;
	    $f5 = $a2e - $a1e ;
	}

#
# test overlap for 'match' criteria
#
        # check for an array of end specific tolerances
	my ($t55, $t53, $t35, $t33 ) ;
	if( ref($tolerance) =~ /ARRAY/) {
	    ($t55, $t53, $t35, $t33 ) = @{$tolerance} ;
	    # if any $tnm is undefined, then that tolerance criterion is ignored?
	    if( ( defined($t55) and (($fw && $a2s<=$a1s && $a2s+$t55 >= $a1s) || (!$fw && $a1e<=$a2e && $a1e+$t55 >= $a2e))) or
		( defined($t53) and (($fw && $a1s<=$a2s && $a1s+$t53 >= $a2s) || (!$fw && $a2e<=$a1e && $a2e+$t53 >= $a1e))) or
		( defined($t35) and (($fw && $a2e<=$a1e && $a2e+$t35 >= $a1e) || (!$fw && $a1s<=$a2s && $a1s+$t35 >= $a2s))) or
		( defined($t33) and (($fw && $a1e<=$a2e && $a1e+$t33 >= $a2e) || (!$fw && $a2s<=$a1s && $a2s+$t33 >= $a1s)))
		) {
		$match = 1;
		$type  = 1;

	    } else { # failure

		$match = 0;
		$type  = 2;

	    }
	} else { # simple scalar tolerance
	    if( # a negative tolerance value accepts all matches
	       $tolerance < 0 or 

	       # accept a match in which both sides are within tolerance
	       (abs($f3) <= $tolerance && abs($f5) <= $tolerance) or

	       # $single => accept matchif either side within tolerance
	       (  $single && 
	         ( abs($f3) <= $tolerance || abs($f5) <= $tolerance ))
               ) {
		$match = 1;
		$type  = 1;

	    } else { # failure

		$match = 0;
		$type  = 2;

	    }
	}
	$match,$f5,$f3,$type,$string;

    }
#####################################################
# No! The two features *do not* absolutely overlap! #
#####################################################
#
# but a tolerance < -1 allows linking of non-overlapping features
# Q: how can I pass back the extent of the mismatch?
#
   elsif($tolerance < -1 && ($mine-$maxs) > $tolerance ){
#       rbsk: changed type from 3 to 4; return $mine, $maxs
	1,$mine,$maxs,4; 

#
# otherwise, match fails!
#
    }
    else{
	0,0,0,0;
    }
}

=pod

=over 4

=item overlap_logical( $GF2, $verbose )

Method that uses GFF::GeneFeature::match() method
to return a non-null ("true") or 0 ("false")
answer on whether or not there is any detectable
'overlap' (on either strand) between pair of
features: the invoking GFF::GeneFeature object
and a second object specified as the first
argument ("$GF2"). If the $verbose flag is
defined and not null, then a detailed match
description is returned.

=back

=cut

#
# Method that uses overlap method to return a simple y/n answer on there being an 'overlap'
# if verbose=1, returns overlap string if there's a match
# $tolerance defaults to -1 unless otherwise specified
#

sub overlap_logical {
    my $self = shift;
    my $other = shift;
    my $verbose = shift;

    my $tolerance = -1 ;
    my $single = 0;
    my $strand = 0;
    my @ret=$self->match($other,$tolerance,$single,$strand);
    if($verbose && $ret[0]){
	return $ret[4];
    }else{
	return $ret[0];
    }
}

=pod

=over 4

=item match_logical( $GF2, $verbose )

Method that uses GFF::GeneFeature::match() method
to return a non-null ("true") or 0 ("false")
answer on whether or not there is an exact
strand/coordinate (overlap) match between pair of
features: the invoking GFF::GeneFeature object
and a second object specified as the first
argument ("$GF2"). If the $verbose flag is
defined and not null, then a detailed match
description is returned.

=back

=cut

sub match_logical {
    my $self = shift;
    my $other = shift;
    my $verbose = shift;

    my $tolerance = 0;
    my $single = 0;
    my $strand = 1;
    my @ret=$self->match($other,$tolerance,$single,$strand);
    if($verbose && $ret[0]){
	return $ret[4];
    }else{
	return $ret[0];
    }
}

=pod

=over 4

=item overlap_merge( $GF2, $tolerance, $strand, $label, $addscores, $copy )

Merge the start and end coordinates of the
invoking GFF::GeneFeature, and a second
GFF::GeneFeature object (given as the "$GF2
object reference argument) to itself, if the two
GFF::GeneFeatures overlap within $tolerance
positions (default 0). 

Optional '$strand' argument forces strand 
sensitivity in merging.

The optional $label argument is overloaded to
two forms:

1. a simple [attributes] field tag may be given, under which  
   the method records merge info.  This includes: 
   'Sequence' tag value (if available), <source>,
   <feature>, <start>, <end> values.

2. a reference to a Perl function may be given, which
   expects to be invoked with the two overlapping 
   gene features in the set. This function would generally
   modify the first gene feature object in some
   customized way.

The optional $addscores argument stipulates that
the scores of merged objects are to be added.

If the $copy argument is set,
then the method makes returns a copy constructed
from the merged object, instead of the original
invoking object (Note: separate copies are made
for each merger event, so a 'self_overlap_merge'
may be useful after this method is called).

Returns all merged (copy of) invoking GeneFeature
iff an overlap_merge occurred, otherwise 0.

=back

=cut

sub overlap_merge {
    my $self      = shift ;
    my $other     = shift ;
    my $tolerance = shift ;
    my $strand    = shift ; # '1' => strand sensitive merge
    my $label     = shift ; 
    my $addscores = shift ;
    my $copy      = shift ; # make a copy of the gene feature?

    $tolerance = 0 if !defined($tolerance) ;
    $strand    = 0 if !defined($strand) ;
    $addscores = 0 if !defined($addscores) ;
    $copy      = 0 if !defined($copy) ;

    my ($s1,$s2,$e1,$e2);
    $s1 = $self->start;
    $s2 = $other->start;
    $e1 = $self->end;
    $e2 = $other->end;

#    print "$s1,$s2,$e1,$e2\n";
    if( (!$strand || ($self->strand eq $other->strand)) 
                     and
        (($s2 >= $s1 && $s2 <= $e1+$tolerance) ||
         ($e2 <= $e1 && $e2+$tolerance >= $s1))
      ) {
	my $gf ;
	if($copy) {
	    $gf = $self->copy() ;
	} else {
	    $gf = $self ;
	}
	if(defined($label)) { 
	    if(ref($label) =~ /CODE/) {
		# ISA reference to a gf1 modifying function?
		&{$label}($self,$other) ;
	    } else { 
		# then simple attribute tag: record some standard data?
		my @newlist ;
		my $sequence = $other->tagValueList('Sequence') ;
		push( @newlist, @{$sequence}) if(defined($sequence)) ;
		push( @newlist, ($other->source(), $other->feature(), $other->start(), $other->end())) ;
		my $oldlist  = $other->tagValueList($label) ;
		push( @newlist, @{$oldlist}) if(defined($oldlist)) ;
		$self->tagValueList($label,
				       \@newlist,
				       1) ; # append rather than overwrite
	    }
	}

	$gf->start(&min($s1,$s2));
	$gf->end(&max($e1,$e2));

	if($addscores) {
	    my $sscore = $self->score() ;
	    my $oscore = $other->score() ;
	    $gf->score($sscore+$oscore) if($sscore ne '.' && $oscore ne '.')  ;
	}

#	print $gf->start." - ".$gf->end."\n";

	return $gf ; # will act boolean in a scalar context?
    } else {
	return 0;
    }
}

=pod

=over 4

=item addMatch( $gfm, $e5, $e3 )

Method to record, for the invoking
GFF::GeneFeature, an overlap match of this
GFF::GeneFeature with another GFF::GeneFeature
object ('$gfm'). The associated 5' and 3' overlap
offsets (see GFF::GeneFeature::match())) are also
recorded.

If the GFF is version 2 or higher, then
the hits are also annotated as [attribute] tag-values
under the tag "'Matches_'.gfm->source".

=back

=cut

sub addMatch {
    my $self  = shift ;
    my $gfm   = shift ;
    my $e5    = shift ; # 5' overlap offset to match
    my $e3    = shift ; # 3' overlap offset to match

    $self->[MATCHES]={} if(!defined($self->[MATCHES])) ;

    $self->[MATCHES]->{$gfm} = [$gfm, $e5, $e3] ; 
    if( $self->version() > 1 ) { # also annotate with tag-values
	my $name = $gfm->Sequence() ;
	$name='?' if !defined($name) ;
	$self->tagValueList('Matches_'.$gfm->source,
			    [$name,$gfm->feature,$e5,$e3],1) ;
    }
}

=pod

=over 4

=item getMatches($source,$feature)

Method to get GFF::GeneFeature match records for
the invoking GFF::GeneFeature object. Returned as
a reference to a hash, keyed on GFF::GeneFeature
object references ('$gfm'), with values equal to
(references to an anonymous array [$gfm, $e5,
$e3] ) of the matched GFF::GeneFeature and
associated offsets (see addMatch above).

Specification of the optional '$source' and
'$feature' Perl regex arguments constrains 
the list returned to matches of a given 
<source> and/or <feature> match respectively.
Note: both $source and $feature are assumed
to be $pattern's intending to be a case insensitive
match the whole field (ie. /^$pattern$/i)

=back

=cut

sub getMatches {
    my $self = shift ;
    my $source = shift ;
    my $feature = shift ;

    return {} if !defined($self->[MATCHES]);

    # matches are returned as a reference 
    # to a hash keyed on GeneFeature references
    my ($Imatches,$Smatches,$Fmatches) ;
    $Imatches = $self->[MATCHES] ; # reference to a hash
    if(defined($source)) {
	$Smatches={} ;
	foreach my $gfm (keys %{$Imatches}) {
	    if( $Imatches->{$gfm}->[0]->source() =~ /^$source$/i) {
		$Smatches->{$gfm} = $Imatches->{$gfm} ;
	    }
	}
    } else { 
	$Smatches = $Imatches ; 
    }
    if(defined($feature)) {
	$Fmatches={} ;
	foreach my $gfm (keys %{$Smatches}) {
	    if( $Smatches->{$gfm}->[0]->feature() =~ /^$feature$/i) {
		$Fmatches->{$gfm} = $Smatches->{$gfm} ;
	    }
	}
    } else { 
	$Fmatches = $Smatches ; 
    }
    return $Fmatches ;
}

=pod

=over 4

=item getMatches_logical()

Method to return non-null ("true") or 0 ("false")
depending upon whether the GFF::GeneFeature has
matches or not.

=back

=cut

sub getMatches_logical {
    my $self = shift;
    return scalar(%{$self->getMatches()});
}

=pod

=over 4

=item addMember( $gffset, $member )

Method to add a "member" record to indicate the
parent GFF $gffset object for a specified grouping
(subset "$member'ship") of GFF::GeneFeatures

=back

=cut

sub addMember {
    my $self = shift;
    my $parent = shift;
    my $member = shift;
    $self->[MEMBER]={} if(!defined($self->[MEMBER])) ;
    $self->[MEMBER]->{$member} = $parent;
}

=pod

=over 4

=item getMember( $memberTAG )

Method to retrieve the reference to the parent
GFF object of GFF::GeneFeature object indexed
under this particular TAGing ("$memberTAG").

=back

=cut

sub getMember {
    my $self = shift;
    my $member = shift;
    return undef if(!defined($self->[MEMBER])) ;
    $self->[MEMBER]->{$member};
}

# misc subroutines:

# max:
# D: return max of i and j
sub max{
    my($i,$j)=@_;
    if($i>$j){
	$i;
    }else{
	$j;
    }
}

# min:
# D: return min of i and j
sub min{
    my($i,$j)=@_;
    if($i<$j){
	$i;
    }else{
	$j;
    }
}

=pod

=head1 Revision History

 3.006 3/7/2000   - rbsk: - fixed bug in which numerically zero (null) tag values are 
                            not properly printed in dumpTags?!?? Don't know how this crept in!

 3.005 22/06/2000 - rbsk: - tagValue() bug: did not properly index a new value into 
                            the value list of a tag if the value list is undefined
                            only set the valuelist->[0] position... Of course,
                            it is bad programming not to lower valuelist->[n] positions 
                            first, before setting the higher positions...

 3.004 05/04/2000 - rbsk: - renameTag() method created
 3.003 14/03/2000 - rbsk: - copy method versioning fixed
 3.002 21/2/2000  - rbsk: - match(): fixed logical bug in ARRAY tolerance processing
                            tagValueList() returns 'undef' for undefined values

 3.001 2/2/2000   - rbsk: - memory consumption was excessive for large feature sets,
                            so I devised an novel implementation of GFF::GeneFeatures
                            to conserve storage space.
			  - Generally use 'attribute' and 'tag' nomenclature now
			  - Treat the 'Sequence' GFF Version 2 attribute as a special field entry
                            (since it is so pervasive)
                          - 'remap' and 'length' made 'zero-less' coordinate axis sensitive

 2.109 21/01/2000 - rbsk: - added reverse complementation to remap() method

 2.108 (31/12/99) - rbsk: - grepTag() and grepTagValue() methods added

 2.107 (08/12/99) - rbsk: - addMatch() method now directly annotates matches
                           as a descriptive Version 2 tag-value set

 2.106 (19/10/99) - rbsk: - match() method can take a reference to an array for
                           the $tolerance value, consisting of pairs of 5' and 3'
                           end specific tolerances.
 2.105 (3/10/99) - rbsk: - $group_tag in overlap_merge() and associated code
                           inherited from GFF::GeneFeatureSet::self_overlap_merge().
                         - group() method bug: fixed Version 1 GFF crash bug
 2.104 30/9/99 - rbsk: - $tag argument in dump(), dump_string(), dump_group()
 2.103 27/9/99 - rbsk: - $strand argument in overlap_merge()
 2.102 21/9/99 - rbsk: - added $copy argument to overlap_merge() method; now 
                         also returns $self instead of simple '1' iff overlap occurs
                       - created the deleteTag() method
   16/9/99 - rbsk: match() method should totally ignore <strand> if $strand not set?
   8/9/99 -  rbsk: overlap_merge() $addscores argument.
   8/9/99 -  rbsk: fixed GeneFeatureSet.pm handling of valueless [group] tags:
                   - group() method can now be used to set tags without values
                     or tags with values, but $gf->group('tag') or 
                     $gf->group('tag','value0','value1',...,'valueN') ; Seems
                     redundant to methods group_value() and group_value_list()
                     which are similar but slight different in their operation.
                   - For Version 2, AUTOLOAD now returns a boolean '1' for any [group] 
                     value tag which exists but has no values. Note that for Version 2,
                     AUTOLOAD names not recognized ALL default to tag methods and
                     fail silently by returning NULL (rather than throwing an exception,
		     as in Version 1 GFF).
  31/8/99 -  rbsk: $append argument added to group_value_list() method;
                   overlap_merge() method: optional '$tolerance' value 
                   provides for overlap merge where the two features lie
                   within $tolerance base pairs of each other
  27/8/99  - rbsk: sub parse_group() bug fix: couldn't parse some
                   instances of 'end' double quotes...
  18/8/99  - rbsk: coded explicit primary field access methods, rather than 
                   relying upon AUTOLOAD (i.e. to gain efficiency - George Hartzwell suggestion :-)
  12/7/99  - rbsk: custom AUTOLOAD recognizes GFF Version 2 [group] tags as access functions
                   i.e. $gf->tagname() == $gf->group_value_list('tagname') 
                   and  $gf->tagname(\@VALUES) == $gf->group_value_list('tagname',\@VALUES) 
   4/5/99  - rbsk: renamed GeneFeature.pm => GFF::GeneFeature.pm
   20/4/99 - rbsk: Instead of a list, the getMatches() method now returns a reference to a hash,
                   keyed on GeneFeature references (getMatches_logical is now strictly boolean)
   19/4/99 - rbsk: GeneFeaturePair.pm functionality merged with GeneFeature.pm (semantic change:
                   Changed 'pairs' to 'matches' in object data structure
                   in order to capture 'one-to-many' semantics:
                       $fields{'Pairs'} => $fields{'matches} 
                       addPair          => addMatch           # adds a matching GeneFeature
                       getPair          => getMatches         # returns array of matches
                       getPair_logical  => getMatches_logical # returns number of matches (0=> false)
                       dump_pairs       => dump_matches       # dumps all GF's with matches (and the matches too)
                    intersect_overlap_pairs() => intersect_overlap_matches:
                    bug found: use 'absolute' coordinate offsets of second relative to first GF
   17/3/99 - rbsk: added $newline arg to dump(), dump_string() & dump_group()
   16/3/99 - rbsk: GeneFeature objects now subclassed from GFFObject class;
   13/3/99 - rbsk: bug fix in group_value() function, dump_group(), et al.
   3/3/99  - rbsk: extensively revised and improved the documentation
                   added Version 2 GFF code, especially, group() field management methods
                   the GFF V.1 GROUP is generally taken to be equal to the GFF V.2 Sequence tag-value

=head1 SEE ALSO

http://www.sanger.ac.uk/Software/GFF/GFF.shtml,
http://www.sanger.ac.uk/Software/GFF/HomolGeneFeature.shtml and
http://www.sanger.ac.uk/Software/GFF/GeneFeatureSet.shtml

=cut

1;  # says use was ok
