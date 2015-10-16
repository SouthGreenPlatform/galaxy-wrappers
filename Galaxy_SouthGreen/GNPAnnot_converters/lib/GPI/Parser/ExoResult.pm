package GPI::Parser::ExoResult;

use strict;
use Carp;

=head1 NAME

 ExoResult.pm

=head1 SYNOPSIS

 use ExoResult;
 $res = GPI::Parser::Exonerate::ExoResult->new();

=head1 DECRIPTION

ExoResult object can used from the parsing of Exonerate result file. The module doesn't parse the file result, but you can create object with data that you parse in the result file.

=head1 VERSION

$Id: ExoResult.pm,v 1.3 2007/03/01 09:28:03 equevill Exp $
$Log: ExoResult.pm,v $
Revision 1.3  2007/03/01 09:28:03  equevill
Added CVS log in the header of the file


=head1 METHODS

=over 1

=item B<new>


 Description : Create a new ExoResult object.
 Arguments    : -
 Return      : An ExoResult object

=item B<set_values>


 Decription  : Save information from the parsing of exonerate result file
               tags are describes in the hash table %TAGS
 Arguments   : the tag of the value, and the value itself
 Return      : -
 Usage       : $ExoResul_object->set_values(-id => "queryID",
                                            -desc => "function");

=item B<id>


 Decription  : Get query identifier

=item B<desc>


 Decription  : Get query decription

=item B<target>


 Decription  : Get target id

=item B<model>


 Decription  : Get used model

=item B<rawscore>


 Decription  : Get rawscore

=item B<query_start>


 Decription  : Get query start

=item B<query_stop>


 Decription  : Get query stop

=item B<target_start>


 Decription  : Get target start

=item B<target_stop>


 Decription  : Get target stop

=item B<alignString>

 Description : Get Alignement for this result in human readable format
 Arguments   : -
 Return      : Return a String with alignement

=item B<get_formatted_alignment>

 Description : Get alignement information in short format
 Arguments   : 'cigar', 'sugar' or 'vulgar'
 Return      : a String that contains the alignement in requested format

=cut
my %TAGS = (
	    "id"          => undef,
	    "desc"        => undef,
	    "target"      => undef,
	    "model"       => undef,
	    "score"       => undef,
	    "query_start" => undef,
	    "query_stop"  => undef,
	    "target_start"=> undef,
	    "target_stop" => undef,
	    "align"       => undef,
	    "sugar"       => undef,
	    "cigar"       => undef,
	    "vulgar"      => undef,
	    "ident"       => undef,
	    "cons"        => undef,
	    "query_length"=> undef,
	    "subject_length" => undef,
	    "cover"       => undef
);

sub new {
  my $class = shift;
  my $self  = {};
  bless($self,$class);
  return $self;
}

sub set_values {
  my $self = shift;
  my $arg;
  my $val;

  while($arg = shift) {
    $val = shift;
    $arg =~ s/^-//;
    if(exists $TAGS{$arg} ) {
      $self->{$arg} = $val;
    }
    else {
      warn "$arg is not a recognized tag of ExoResult class\n";
    }
  }
}

sub id {
  my $self = shift;
  if(exists $self->{id}) {
    return $self->{id};
  }
  else {
    return $TAGS{id};
  }
}
sub desc{
  my $self = shift;
  if(exists $self->{desc}) {
    return $self->{desc};
  }
  else {
    return $TAGS{desc};
  }
}

sub target{
  my $self = shift;
  if(exists $self->{target}) {
    return $self->{target};
  }
  else {
    return $TAGS{target};
  }
}

sub model{
  my $self = shift;
  if(exists $self->{model}) {
    return $self->{model};
  }
  else {
    return $TAGS{model};
  }
}

sub score{
  my $self = shift;
  if(exists $self->{score}) {
    return $self->{score};
  }
  else {
    return $TAGS{score};
  }
}

sub query_start{
  my $self = shift;
  if(exists $self->{query_start}) {
    return $self->{query_start};
  }
  else {
    return $TAGS{query_start};
  }
}

sub query_stop {
  my $self = shift;
  if(exists $self->{query_stop}) {
    return $self->{query_stop};
  }
  else {
    return $TAGS{query_stop};
  }
}

sub target_start{
  my $self = shift;
  if(exists $self->{target_start}) {
    return $self->{target_start};
  }
  else {
    return $TAGS{target_start};
  }
}

sub target_stop {
  my $self = shift;
  if(exists $self->{target_stop}) {
    return $self->{target_stop};
  }
  else {
    return $TAGS{target_stop};
  }
}

sub pident {
    my $self = shift;
  if(exists $self->{ident}) {
    return $self->{ident};
  }
  else {
    return $TAGS{ident};
  }
}

sub pcons {
    my $self = shift;
  if(exists $self->{cons}) {
    return $self->{cons};
  }
  else {
    return $TAGS{cons};
  }
}

sub query_length {
    my $self = shift;
  if(exists $self->{query_length}) {
    return $self->{query_length};
  }
  else {
    return $TAGS{query_length};
  }
}

sub subject_length {
    my $self = shift;
  if(exists $self->{subject_length}) {
    return $self->{subject_length};
  }
  else {
    return $TAGS{subject_length};
  }
}

sub pcover {
    my $self = shift;
  if(exists $self->{cover}) {
    return $self->{cover};
  }
  else {
    return $TAGS{cover};
  }
}

sub _align {
  my $self = shift;
  my $raw = shift;

  $self->{align} .= $raw."\n";
}

sub alignString {
  my $self = shift;

  if(exists $self->{align}) {
    return $self->{align};
  }
  else {
    return $TAGS{align}
  }
}

sub get_formatted_alignment {
  my $self = shift;
  my $type = shift;

  if($type !~ /^(cigar|vulgar|sugar)$/) {
    warn "alignment format you try to get in get_formatted_alignment() doesn't exists\nSelect one of this format : cigar, sugar or cigar\n";
    return undef;
  }
  else {
    return $self->{$type};
  }
}
1;
