=head1 NAME

HmmTop - Module to parse and create GFF3 output from HmmTop results

=head1 SYNOPSIS

#example here
use GPI::Parser::HmmTop;

=head1 DESCRIPTION

It uses BioPerl modules to parse different output results.

=head1 VERSIONS

$Id: HmmTop.pm,v 1.4 2007/03/01 09:28:03 equevill Exp $
$Log: HmmTop.pm,v $
Revision 1.4  2007/03/01 09:28:03  equevill
Added CVS log in the header of the file


Copyright (c) INRA/URGI 2005

=head1 AUTHORS / ACKNOWLEDGEMENTS

Emmanuel Quevillon <emmanuel.quevillon@versailles.inra.fr>

=cut

package GPI::Parser::HmmTop;

### FIXME : This module is not yet installed in the perl lib, so check were we are. ###


use strict;
use English;
use Data::Dumper;
use vars qw(@ISA);
use GPI::Parser;
@ISA = qw( GPI::Parser );
use lib '/apps/GnpAnnot/lib';
use GPI::Bio::Tools::HmmTop;
use GPI::GFF;

=head1 new

Descrption: Create a new GPI::Parser::HmmTop object.
 Arguments: $params reference to hash table with parameters.
   Returns: Reference to GPI::Parser::HmmTop object or undef.

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

    $self->usebioperl(1);

    return $self;

}



=head1 parse

Description: Read a HmmTop result using GPI::Bio::Tools::HmmTop parsing object.
  Arguments: $in HmmTop results file.
    Returns: 0, '' on success
             1, message on error

=cut

sub parse ($){

    require Bio::SeqIO;
    require GPI::Bio::Tools::HmmTop;

    my ($self, $in) = @_;

    $self->_exitOnError("Fitlers not supported yet for HmmTop.") if($self->haveFilters());

    unless($self->_checkApplication($self->appl())){
	my $list = join(', ', sort keys %${GPI::ApplInfo2GFFFeat});
	$self->_exitOnError("Application name not set or not supported!: \n\tSupported: $list\n");
    }

    my $writer;
    my $features;
    my $bioperl = $self->usebioperl();

    #Array reference for Bio::Seq objects
    my $seqs = [ ];

    if($bioperl){

	require GPI::GFF;

	$features = [ ];

	if($self->flush()){
	    $writer = GPI::GFF->new(
						-gff_version => 3,
						-fh          => $self->fh() || \*STDOUT
					       );
	}
    }else{
	$features = { };
    }

    unless($in){
	$self->_exitOnError("No input file given to parse");
    }

    my $hmmtop = GPI::Bio::Tools::HmmTop->new(-file => $in);

    while(my $prediction = $hmmtop->next_prediction()){

	#### TODO Paly with results and predictions;
	next if($prediction->get_Annotations('transmembraneHelixNumber')->value() == 0);

	#### TODO Check if locations are present, if so, add _part features
	if($self->flush()){
	    $writer->write_feature($prediction);
	}
	else{
	    push @$features, $prediction;
	}
    }


    return($features);
}
