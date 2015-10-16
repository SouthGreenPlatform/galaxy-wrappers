=head1 NAME

Garnier - Module to Garnier (Proteins Secondary structure identifier) results.

=head1 SYNOPSIS

#example here
use GPI::Parser::Garnier;

=head1 DESCRIPTION

It parses Garnier (from EMBOSS package) output resutls.

=head1 VERSIONS

$Id: Garnier.pm,v 1.6 2007/03/01 09:28:03 equevill Exp $
$Log: Garnier.pm,v $
Revision 1.6  2007/03/01 09:28:03  equevill
Added CVS log in the header of the file


Copyright (c) INRA/URGI 2005

=head1 AUTHORS / ACKNOWLEDGEMENTS

Emmanuel Quevillon <emmanuel.quevillon@versailles.inra.fr>

=cut

package GPI::Parser::Garnier;
use lib '/apps/GnpAnnot/lib';

use strict;
use English;
use Data::Dumper;
use vars qw(@ISA);
use GPI::Parser;
@ISA = qw( GPI::Parser );


=head1 new

Descrption: Create a new GPI::Parser::Garnier object.
 Arguments: $params reference to hash table with parameters.
   Returns: Reference to GPI::Parser::Garnier object or undef.

=cut

sub new ($;$$){

    my($class, @params) = @_;

    my $self = bless { }, $class;

    my %params =  $self->_init(\@params);

    $self->setParams(\%params);

    if($self->getParam('debug')){
	for my $key (sort keys %{$self->getParams()}){
	    print "[DEBUG] $key => $params{$key}\n";
	}
    }

    return $self;

}



=head1 parse

Description: Parse Garnier output results
             object.
  Arguments: $in Garnier results file.
    Returns: 0, $hres on success
             1, message on error

=cut

sub parse ($){

    my ($self, $ifile) = @_;

    unless($ifile){
	$self->_exitOnError("No input file given to parse");
    }

    #unless($hres && ref($hres) eq 'HASH'){
    #    $self->_exitOnError("Need a reference to a hash table.");
    #}
    my $hres = { };

    unless( -f $ifile ){
	close($self->getParam('fh')) if($self->isFileHandle($self->getParam('fh')));
	$self->_exitOnError("[$ifile] does not exist : $!");
    }

    unless(open(IN, "<$ifile")){

	close($self->getParam('fh')) if($self->isFileHandle($self->getParam('fh')));
	$self->_exitOnError("Could not open $ifile : $!");
    }

    my ($topSeen, $fileName);
    my ($seqLength, $seqID, $seqDesc);
    my $nresults = 0;

    RES: while(<IN>){

	  chomp;

	  next if($_ =~ /^\s*$/);

	  if($_ =~ /garnier plot of (\S+),\s*(\d+) aa;/){
	      ++$nresults;
	      $fileName = $1;
	      $seqLength = $2;
	      $topSeen = 1;
	      next;
	  }
	  elsif($topSeen && $_ =~ /^\s+(\S+)\s+(.+)/){
	      $seqID = $1;
	      $seqDesc = $2;

	      $hres->{$nresults}->{topdesc}->{title} = "##sequence-region";
	      $hres->{$nresults}->{topdesc}->{name}  = $seqID ? $seqID : "Region_$nresults";
	      $hres->{$nresults}->{topdesc}->{name}  =~ s/.+\/([^\/]+)$/$1/;
	      $hres->{$nresults}->{topdesc}->{start} = 1;
	      $hres->{$nresults}->{topdesc}->{end}   = $seqLength;

	      my $nhits = 0;

	    HIT: while(<IN>){

		  chomp;

		  #Given positions
		  my ($seqStart, $seqEnd) = (0,0);

		  #my $hits = { };
		  #$hits->{name} = $seqID;

		  if($_ =~ /percent/){
		      $topSeen = 0;
		      next RES;
		  }
		  elsif($_ =~ /\s+(\d+)$/){ #    .   10    .   20    .   30    .   40    .   50    .   60
		      $seqEnd = $1;
		      $seqStart = (($seqEnd - 60) + 1);
		      my $nhsps = 0;
		      #++$nhits;
		      #my $hits = { };
		      #$hits->{name} = $seqID;
		      #$hits->{id} = $nhits;

		    HSP: while(<IN>){

			  chomp;

			  if(/^\s*$/){
			      #$hres->{$nresults}->{hits}->{$nhits} = $hits;
			      next HIT;
			  }



			  if($_ =~ /(helix)\s(.+)$/){
			      #We define each secondary structure as a hsp
			      my @hsp;
			      ++$nhits;
			      my $hits = { };
			      #$hits->{name} = $seqID;# . "_$1";# . sprintf('%0' . $NAMELEN . 'd', $nhits);
			      $hits->{extra} = $1;
			      $hits->{id} = $nhits;

			      my $extra = {
					   starts => $seqStart,
					   ends   => $seqEnd,
					   name   => $seqID # . '_helix' #We don't need name because there is no Target, no a seq/seq comparison.
					  };
			      my $string = $2;
			      my($res, $msg) = $self->getPositions($string, 'H', \@hsp, \$nhsps, $extra);
			      $self->_exitOnError($msg) if($res);
			      push(@{$hits->{hsps}}, @hsp);
			      $hres->{$nresults}->{hits}->{$nhits} = $hits;
			  }
			  elsif($_ =~ /(sheet)\s(.+)$/){
			      #We define each secondary structure as a hsp
			      my @hsp;
			      ++$nhits;
			      my $hits = { };
			      #$hits->{name} = $seqID;# . "_$1" ;# . sprintf('%0' . $NAMELEN . 'd', $nhits);;
			      $hits->{extra} = $1;
			      $hits->{id} = $nhits;

			      my $extra = {
					   starts => $seqStart,
					   ends   => $seqEnd,
					   name   => $seqID # . '_sheet'
					  };
			      my $string = $2;
			      my($res, $msg) = $self->getPositions($string, 'E', \@hsp, \$nhsps, $extra);
			      $self->_exitOnError($msg) if($res);
			      push(@{$hits->{hsps}}, @hsp);
			      $hres->{$nresults}->{hits}->{$nhits} = $hits;
			  }
			  elsif($_ =~ /(turns)\s(.+)$/){
			      #We define each secondary structure as a hsp
			      my @hsp;
			      ++$nhits;
			      my $hits = { };
			      #$hits->{name} = $seqID;# . "_$1";# . sprintf('%0' . $NAMELEN . 'd', $nhits);
			      $hits->{extra} = $1;
			      $hits->{id} = $nhits;

			      my $extra = {
					   starts => $seqStart,
					   ends   => $seqEnd,
					   name   => $seqID # . '_turns'
					  };
			      my $string = $2;
			      my($res, $msg) = $self->getPositions($string, 'T', \@hsp, \$nhsps, $extra);
			      $self->_exitOnError($msg) if($res);
			      push(@{$hits->{hsps}}, @hsp);
			      $hres->{$nresults}->{hits}->{$nhits} = $hits;
			  }
			  elsif($_ =~ /(coil)\s{2}(.+)$/){ #coil is one letter long less than sheet/helix or turns
			      #We define each secondary structure as a hsp
			      my @hsp;
			      ++$nhits;
			      my $hits = { };
			      #$hits->{name} = $seqID;# . "_$1";# . sprintf('%0' . $NAMELEN . 'd', $nhits);
			      $hits->{extra} = $1;
			      $hits->{id} = $nhits;

			      my $extra = {
					   starts => $seqStart,
					   ends   => $seqEnd,
					   name   => $seqID # . '_coil'
					  };
			      my $string = $2;
			      my($res, $msg) = $self->getPositions($string, 'C', \@hsp, \$nhsps, $extra);
			      $self->_exitOnError($msg) if($res);
			      push(@{$hits->{hsps}}, @hsp);
			      $hres->{$nresults}->{hits}->{$nhits} = $hits;
			  }

		      } #end while(<IN>) #3

		  } # if($_ =~ /^\s+\.

	      } #end while #2

	  } #end elsif topSeen

	  else{
	      $self->_exitOnError("Could not find the beginning of the file?");
	  }
      } #end while #1

    $hres->{nresults} = $nresults - 1;

#    print Dumper $hres;
 #   exit;

    return($hres);

}

=head1 getPositions

Descrpition: Get start and stop positions for a specific substring in another string
  Arguments: $string String to parse
             $sub    Substring to search
             $hsp    Reference to an array
             $nhsps  Current hsps number
             $extra  Reference to a hash table
    Returns: 0, ''  on success
             1, msg on errors

=cut

sub getPositions ($$$;$$$$){

    my($self, $string, $sub, $hsp, $nhsps, $extra) = @_;

    return (0, '') unless(length($string));
    return (1, "No substring given to search.") unless($sub && length($sub));
    return (1, "Need a reference to an array.") unless($hsp && ref($hsp) eq 'ARRAY');

    #my $NameLen = $NAMELEN;
    my $curPos = 1;
    my $ok = 0;
    my $starts = undef;
    my $ends   = undef;

    foreach my $char (split('', $string)){
	$ok = $char eq $sub ? 1 : 0;

	if($ok){
	    unless(defined($starts)){
		$starts = $curPos;
	    }
	}else{
	    $ends = $curPos;
	    if(defined($starts)){
		my $hash = { };
		++${$nhsps};
		%$hash = %$extra if($extra);

		$hash->{startq} = $starts;
		$hash->{startq} += $extra->{starts} -1 if($extra);
		$hash->{endq}   = $ends -1;
		$hash->{endq}   += $extra->{starts} -1 if($extra);
		$hash->{id}     = $$nhsps;

		$hash->{starts} = $hash->{startq};
		$hash->{ends}   = $hash->{endq};
		#$hash->{name}   .= sprintf('%0'.$NameLen.'d', $$nhsps);

		push(@$hsp, $hash) if(defined($starts));
		$starts = undef;
	    }
	}
	$curPos++;
    }

    #In case the last letter is a substring.
    if(defined($starts)){
	my $hash = { };
	++${$nhsps};
	%$hash = %$extra if($extra);

	$hash->{startq} = $starts;
	$hash->{startq} += $extra->{starts} -1 if($extra);
	$hash->{endq}   = $curPos -1;
	$hash->{endq}   += $extra->{starts} -1 if($extra);
	$hash->{id}     = $$nhsps;

	$hash->{starts} = $hash->{startq};
	$hash->{ends}   = $hash->{endq};
	#$hash->{name}   .= sprintf('%0'.$NameLen.'d', $$nhsps);

	push(@$hsp, $hash) if(defined($starts));
	$starts = undef;
    }

    return(0, '');
}

1;
