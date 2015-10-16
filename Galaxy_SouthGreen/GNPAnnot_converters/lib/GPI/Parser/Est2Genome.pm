=head1 NAME

Est2Genome - Est2Genome parser utility

=head1 SYNOPSIS

#example here
use GPI::Parser::Est2Genome;

=head1 DESCRIPTION

This module parses a est2genome (EMBOSS) output.
Either the aliagnment is present or not, it is not important.

=head1 VERSIONS

$Id: Est2Genome.pm,v 1.3 2007/03/01 09:28:03 equevill Exp $
$Log: Est2Genome.pm,v $
Revision 1.3  2007/03/01 09:28:03  equevill
Added CVS log in the header of the file


Copyright (c) INRA/URGI 2005

=head1 AUTHORS / ACKNOWLEDGEMENTS

Emmanuel Quevillon <emmanuel.quevillon@versailles.inra.fr>

=cut

package GPI::Parser::Est2Genome;
use lib '/apps/GnpAnnot/lib';
use strict;
use Data::Dumper;
use vars qw(@ISA);
use GPI::Parser;
use GPI::GFF;
use Bio::SeqFeature::Generic;
@ISA = qw( GPI::Parser );

my $TYPES = {
	     'Span'    => 'match_set',
	     'Segment' => 'match_part',
	    };


=head1 new

Descrption: Create a new GPI::Parser::Est2Genome object.
 Arguments: $params reference to hash table with parameters.
   Returns: Reference to GPI::Parser::Est2Genome object or undef.

=cut

sub new ($;$$){

    my($class, @params) = @_;

    my $self = bless { }, $class;

    my %params =  $self->_init(\@params);

    $self->setParams(\%params);

    if($self->debug()){
	for my $key (sort keys %{$self->getParams()}){
	    print STDERR "[DEBUG] $key => $params{$key}\n";
	}
    }

    return $self;
}



=head1 parse

Description: Parse est2genome output
  Arguments: $in Est2genome results file.
    Returns: 0, $features on success
             1, message on error

=cut

sub parse ($){

    my ($self, $ifile) = @_;

    unless($ifile){
	$self->_exitOnError("No input file given to parse");
    }

    unless($self->_checkApplication($self->appl())){
	my $list = join(', ', sort keys %${GPI::ApplInfo2GFFFeat});
	$self->_exitOnError("Application name not set or not supported!: \n\tSupported: $list\n");
    }

    my $writer;
    my $features = [ ];
    my $ofh;

    $self->_print_verbose("Initilalising reader Bio::Root::IO object with file [$ifile] ...");
    #A generic reader using Bio::Root::IO to access _readline and _pushback methods
    my $reader = Bio::Root::IO->new(
				    -file => $ifile,
				   );

    $reader || $self->_exitOnError("Could not open $ifile : $!\n");

    if($self->flush()){
	$self->_print_verbose("Initialising writer GPI::GFF ...");
	$writer = GPI::GFF->new(
				       -gff_version => 3,
				       -fh          => $self->fh() || \*STDOUT
				      );
	$ofh = $self->fh() || \*STDOUT;
    }

    my $span = 0;
    my $seg  = 0;
    my $ID;
    my $PI = $self->pidentity() || 0;
    my $currest = undef;
    my $strand  = undef;
    my $ok      = 0;
    #FIXME : Pass this value as an argument n the command line.(Length of the EST).
    my $eLength = 60;

    $self->_print_verbose("Reading input file ...");

    while(defined(my $line = $reader->_readline())){

	next if $line =~ /^\s*$|^Note Best/;

	#If we have some Intron defined, we can get the strand.
	if($line =~ /^(.)Intron/){
	    #$strand = $1;
	    next;
	}

	my($type, $score, $ident, $start, $end, $seq_id, $tstart, $tend, $tid, @tname) = split(/\s+/, $line);

	if($currest eq $tid){

	    if($line =~ /^Span/){

		if($ident >= $PI){

		    my $ESTLENGTH = $tend - $tstart + 1;
		    next if($ESTLENGTH != $eLength);

		    $span++;
		    $seg = 0;
		    $ID = join("_", $seq_id, $tid, $TYPES->{$type}.sprintf("%04d", $span));
		    my $tags = {
				'target_start' => $tstart,
				'target_end'   => $tend,
				'target_id'    => $tid,
				'target_description' => join(" ", @tname),
				'ID'           => $ID,
				'Name'         => join("_", $seq_id, $tid),
				'Target'       => join("+", $tid, $tstart, $tend),
				'program'      => $self->appl() || 'est2genome',
				'target_pident'=> $ident,
			       };

		    my $feat = Bio::SeqFeature::Generic->new(
							     -seq_id      => $seq_id,
							     -source      => $self->appl() || 'est2genome',
							     -primary_tag => $TYPES->{$type},
							     -start       => $start,
							     -end         => $end,
							     -score       => $score,
							     -strand      => $strand || '.',
							     -tag         => $tags
							    );

		    push @$features, $feat;
		    $ok = 1;
		}
	    }
	    elsif($line =~ /^Segment/){

		next unless $ok;

		$seg++;
		my($type, $score, $ident, $start, $end, $seq_id, $tstart, $tend, $tid, @tname) = split(/\s+/, $line);
		my $id = join("_", $ID, $TYPES->{$type}.sprintf("%04d", $seg));
		
		my $tags = {
			    'ID'           => $id,
			    'Name'         => join("_", $seq_id, $tid),
			    'Target'       => join("+", $tid, $tstart, $tend),
			    'program'      => $self->appl() || 'est2genome',
			    'Parent'       => $ID,
			   };

		my $feat = Bio::SeqFeature::Generic->new(
							 -seq_id      => $seq_id,
							 -source      => $self->appl() || 'est2genome',
							 -primary_tag => $TYPES->{$type},
							 -start       => $start,
							 -end         => $end,
							 -score       => $score,
							 -strand      => $strand || '.',
							 -tag         => $tags
							);
		push @$features, $feat;
	    }
	}
	else{
	    if(defined($currest) && $self->flush()){
		if(@$features){
		    $writer->write_feature(@$features);
		    print $ofh "###\n";
		    $features = [ ];
		}
	    }
	    $currest = $tid;
	    $ok      = 0;
	    $strand  = undef;
	    $reader->_pushback($line);
	}	
    }

    $reader->close();
    $writer->close() if($self->flush());

    return($features);

}
