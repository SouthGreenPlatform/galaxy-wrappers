=head1 NAME

Fasta - Module to create Fasta GFF3 formated output.

=head1 SYNOPSIS

#example here
use GPI::Parser::Fasta;

=head1 DESCRIPTION

It uses BioPerl modules to parse different output results.
(Bio::SeqIO)

=head1 VERSIONS

$Id: Fasta.pm,v 1.7 2007/03/01 09:28:03 equevill Exp $
$Log: Fasta.pm,v $
Revision 1.7  2007/03/01 09:28:03  equevill
Added CVS log in the header of the file


Copyright (c) INRA/URGI 2005

=head1 AUTHORS / ACKNOWLEDGEMENTS

Emmanuel Quevillon <emmanuel.quevillon@versailles.inra.fr>

=cut

package GPI::Parser::Fasta;

use lib '/apps/GnpAnnot/lib';
use strict;
use English;
use Data::Dumper;
use vars qw(@ISA);
use GPI::Parser;
@ISA = qw( GPI::Parser );


=head1 new

Descrption: Create a new GPI::Parser::Fasta object.
 Arguments: $params reference to hash table with parameters.
   Returns: Reference to GPI::Parser::Fasta object or undef.

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

    unless($self->source()){
	print STDERR "You need to provid a source name to create Fasta GFF3.\n";
	exit 1;
    }

    return $self;

}



=head1 parse

Description: Read a Fasta file through a Bio::SeqIO object.
  Arguments: $in Fasta file.
    Returns: 0, '' on success
             1, message on error

=cut

sub parse ($){

    require Bio::SeqIO;

    my ($self, $in) = @_;

    my $seqi = Bio::SeqIO->new(
			       -file   => "$in",
			       -format => "fasta"
			      );

    my $seqo = Bio::SeqIO->new(
			       -format => 'fasta',
			       -fh     => $self->fh(),
			      );

    unless($seqi && $seqo){
	unlink($self->out()) if($self->out());
	close($self->fh()) if($self->isFileHandle($self->fh()));
	$self->_exitOnError("Could not create Bio::SeqIO object.");
    }

    my $ofh = $self->fh();

    unless($self->isFileHandle($ofh)){
        $ofh = \*STDOUT;
    }

    my $delim = "\t";

    #We don't need this as the loader will complain if there is no sequence definition
    print $ofh "##gff-version 3\n";

    #This array will keep all the Sequence object to add them at the end of the file.
    my @seqs;

    while(my $seq = $seqi->next_seq()){
	print $ofh join($delim, $seq->id, $self->source(), $self->getApplInfo('region'), '1', $seq->length(), ".", "+", ".", "ID=".$seq->id()), "\n";
	push @seqs, $seq;
	#$seqo->write_seq($seq), "\n";
    }

    print $ofh "##FASTA\n";

    for(@seqs){
	$seqo->write_seq($_);
    }

    return;
}
