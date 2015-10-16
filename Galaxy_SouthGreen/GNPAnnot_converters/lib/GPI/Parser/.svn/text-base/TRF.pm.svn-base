=head1 NAME

TRF - Module to parse Tandem Repeat Finder results. At the moment this module only supports non HTML results.

=head1 SYNOPSIS

#example here
use GPI::Parser::TRF;

=head1 DESCRIPTION

It uses BioPerl modules to parse different output results.
(GPI::Bio::Tools::TRF)

=head1 VERSIONS

$Id: TRF.pm,v 1.14 2007/03/01 09:28:04 equevill Exp $
$Log: TRF.pm,v $
Revision 1.14  2007/03/01 09:28:04  equevill
Added CVS log in the header of the file



Copyright (c) INRA/URGI 2006

=head1 AUTHORS / ACKNOWLEDGEMENTS

Emmanuel Quevillon <emmanuel.quevillon@versailles.inra.fr>

=cut

### FIXME : This module is not yet installed in the perl lib, so check were we are. ###

my $lib;

BEGIN {
	if($ENV{USER} eq 'tuco'){
		$lib = '/home/tuco/src/perl/lib/GPI-Bio/lib';
	}else{
		$lib = '/home/equevill/src/perl/lib/GPI-Bio/lib';
	}	

}


package GPI::Parser::TRF;

use strict;
use English;
use Data::Dumper;
use vars qw(@ISA);
use GPI::Parser;
use lib $lib;
use GPI::Bio::Tools::TRF;
use GPI::GFF;

@ISA = qw( GPI::Parser );


=head1 new

Descrption: Create a new GPI::Parser::Blast object.
 Arguments: $params reference to hash table with parameters.
   Returns: Reference to GPI::Parser::Blast object or undef.

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

    if(defined($self->evalue())){

	my $evalue = $self->evalue();

	#We just create this object to reformat the evalue.
	my $mathfloat = Math::BigFloat->new($evalue);
	return undef unless(defined($mathfloat));

	#We set the new formatted evalue to access it directly in the future.
	$self->evalue($mathfloat->bstr());

    }

    #We force usage of bioperl to be sure it will use GPI::Bio::Tools::GFF
    #Also flushing will use it, we won't have to modify GPI::Parser.
    $self->usebioperl(1);
    $self->flush(1);

    return $self;

}



=head1 parse

Description: Parse Tandem Repeat Finder results using GPI::Bio::Tools::TRF
             object.
  Arguments: $in TRF results file (.dat), non HTML output.
    Returns: 0, $features on success
             1, message on error

=cut

sub parse ($){

    require GPI::Bio::Tools::TRF;
    require Bio::Seq;

    my ($self, $ifile) = @_;

    $self->_exitOnError("Fitlers not supported yet for Tandem Repeat Finder.") if($self->haveFilters());

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
	require Bio::SeqFeature::Generic;
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

    unless($ifile){
	$self->_exitOnError("No input file given to parse");
    }

    my $trf = GPI::Bio::Tools::TRF->new(-file   => "$ifile");

    unless($trf){
	unlink($self->out()) if($self->out());
	close($self->fh()) if($self->isFileHandle($self->fh()));
	$self->_exitOnError("Could not create GPI::Bio::Tools::TRF object.");
    }

    local $OUTPUT_AUTOFLUSH = 1;

    my $maxhits = $self->maxhits();
    my($nresults, $nfeats) = (0, 0);


    while(my $result = $trf->next_result()){
	$nresults++;

        next unless($result->num_features());

	my ($progress, $next_update, $pbhits) = (0, 0, 0);
	
	my $RefTag = { };
	$RefTag->{ID} = $result->analysis_query() || $self->getApplInfo('region') . sprintf('%0'.$NAMELEN.'d', $nresults);

	my $RefFeat =  Bio::SeqFeature::Generic->new(
						     -seq_id      => $result->analysis_query() || "Tandem_Repeat_Finder_$nresults",
						     -source_tag  => $result->analysis_method(),
						     -primary_tag => $self->qtype() || 'region',
						     -start       => 1,        #$result->_get_features()->[0]->start(),                         #FIXME
						     -end         => $result->query_length(),   #By default this the end of the last feature
						     -score       => '.' ,                      #Right now we cannot retrieve the seq length!
						     -strand      => '+',
						     -frame       => '.',
						     -tag         => $RefTag
						    );

	unless($self->flush()){
	    push(@$features, $RefFeat);
	}else{
	    $writer->write_feature($RefFeat);
	}


	if($self->verbose()){
	    $progress = Term::ProgressBar->new({count => $result->num_features(),
						name  => "Reading result $nresults [" . $result->analysis_query() . "]"}
					      );
	}

	while(my $feat = $result->next_feature()){
	    $nfeats++;

	    $pbhits++;
	    $next_update = $progress->update($pbhits) if($pbhits >= $next_update && $self->verbose());

	    my $tag = { };

	    my $ID             = ($feat->get_tag_values('ID'))[0];
	    my $consensusSeq   = ($feat->get_tag_values('consensusSeq'))[0];
	    my $periodSize     = ($feat->get_tag_values('periodSize'))[0];
	    my $copyNumber     = ($feat->get_tag_values('copyNumber'))[0];
	    my $consensusSize  = ($feat->get_tag_values('consensusSize'))[0];
	    my $percentMatches = ($feat->get_tag_values('percentMatches'))[0];
	    my $percentIndels  = ($feat->get_tag_values('percentIndels'))[0];
	    my $percentA       = ($feat->get_tag_values('percentA'))[0];
	    my $percentC       = ($feat->get_tag_values('percentC'))[0];
	    my $percentG       = ($feat->get_tag_values('percentG'))[0];
	    my $percentT       = ($feat->get_tag_values('percentT'))[0];
	    my $entropy        = ($feat->get_tag_values('entropy'))[0];


	    my $note = "periodsize:$periodSize,copynumber:$copyNumber,percentmatches:$percentMatches,percentindles:$percentIndels,consensussize:$consensusSize";
	    $note   .= ",percentA:$percentA,percentC:$percentC,percentG:$percentG,percentT:$percentT";
	    $note   .= ",entropy:$entropy" if(defined($entropy));
	    #We create a new fasta entry
	    my $targetName = join("_", $feat->seq_id(), "ConsSeq", "TRF".sprintf("%03d", $nfeats));

	    my $seq = Bio::Seq->new();
	    $seq->id($targetName);
	    $seq->seq($consensusSeq);

	    push(@$seqs, $seq);

	    #Note values are set with %5D ...
	    $tag->{'Note'}   = $note;
	    $tag->{'Target'} = $targetName . "+1+" . length($consensusSeq);
	    $tag->{'ID'}     = $ID;

	    $feat->remove_tag();
	    $feat->set_attributes(-tag => $tag);

	    unless($self->flush()){
		push(@$features, $feat);
	    }else{
		print $writer->write_feature($feat);
	    }

	}
	#Progress Bar, to complete to 100% the progress bar
	$progress->update($result->num_features()) if($self->verbose());# if($result->num_hits >= $next_update);

    }

    #Print FASTA sequences
    my($res, $mess) = $self->write_sequences($seqs);
    $self->_exitOnError($mess) if($res);

    #To avoid writing twice the sequences
    $self->addseq(0);

    return ($features);

}
