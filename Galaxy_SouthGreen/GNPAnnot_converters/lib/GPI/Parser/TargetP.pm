=head1 NAME

TargetP - Module to Parse and create GFF3 output from TargetP results

=head1 SYNOPSIS

#example here
use GPI::Parser::TargetP;

=head1 DESCRIPTION

It uses BioPerl modules to parse different output results.
(Bio::SeqIO)

=head1 VERSIONS

$Id: TargetP.pm,v 1.4 2007/03/01 09:28:04 equevill Exp $
$Log: TargetP.pm,v $
Revision 1.4  2007/03/01 09:28:04  equevill
Added CVS log in the header of the file


Copyright (c) INRA/URGI 2005

=head1 AUTHORS / ACKNOWLEDGEMENTS

Emmanuel Quevillon <emmanuel.quevillon@versailles.inra.fr>

=cut

package GPI::Parser::TargetP;

### FIXME : This module is not yet installed in the perl lib, so check were we are. ###

my $lib;

BEGIN {
	if($ENV{USER} eq 'tuco'){
		$lib = '/home/tuco/src/perl/lib/GPI-Bio/lib';
	}else{
		$lib = '/home/equevill/src/perl/lib/GPI-Bio/lib';
	}	

}


use strict;
use English;
use Data::Dumper;
use vars qw(@ISA);
use GPI::Parser;
@ISA = qw( GPI::Parser );
use lib $lib;
use GPI::Bio::Tools::TargetP;
use GPI::GFF;

=head1 new

Descrption: Create a new GPI::Parser::TargetP object.
 Arguments: $params reference to hash table with parameters.
   Returns: Reference to GPI::Parser::TargetP object or undef.

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

Description: Read a TargetP results file using a GPI::Bio::Tools::TargetP object
  Arguments: $in TargetP output results file
    Returns: 0, '' on success
             1, message on error

=cut

sub parse ($){

    require GPI::Bio::Tools::TargetP;
    require Bio::Annotation::SimpleValue;

    my ($self, $in) = @_;

    $self->_exitOnError("Fitlers not supported yet for TargetP.") if($self->haveFilters());

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
    
    my $targetp = GPI::Bio::Tools::TargetP->new(-file => $in);
    my $predcnt = 1;

    while(my $prediction = $targetp->next_prediction()){
	
	#Shall we ignore 'Any other location' and keep only real detected location?
	my $location = $prediction->get_Annotations('location');
	next if($location =~ /Any other|Unknown/);

	if($location =~ /Mitochondrion/){
	    $prediction->score($prediction->get_Annotations('mitochondrionCutOff'));
	}
	elsif($location =~ /Secretory/){
	    $prediction->score($prediction->get_Annotations('signalPeptideCutOff'));
	}
	elsif($location =~ /Chloroplast/){
	    $prediction->score($prediction->get_Annotaions('chloroplastCutOff'));
	}
	else{
	    return (1, "Location [$location] is not supported or unknown.");
	}

	my $ID = join("_", $self->appl(), $prediction->seq_id(), $self->getApplInfo('part') . sprintf('%0'.$GPI::Parser::NAMELEN.'d', $predcnt));

	#$prediction->add_Annotation('ID', new Bio::Annotation::SimpleValue(-value => $ID));
	$prediction->add_tag_value('ID', $ID);
	$prediction->add_tag_value('programm', 'targetp');
	$prediction->add_tag_value('lib', 'matrix');
	$prediction->add_tag_value('Dbxref', 'TargetP:targetp');
	#We set the end of the prediction to the end of the signal
	my $signalLength = $prediction->get_Annotations('signalPeptideLength');
	$prediction->end($signalLength);
	
        $predcnt++;

	if($self->flush()){
	    print $writer->write_feature($prediction);
	}
	else{
	    push @$features, $prediction;
	}
    }

    return($features);
}
