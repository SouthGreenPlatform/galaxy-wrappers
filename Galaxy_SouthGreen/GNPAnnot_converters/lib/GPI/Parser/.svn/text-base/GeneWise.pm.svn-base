=head1 NAME

GeneWise - Module to parse GeneWiseresults.

=head1 SYNOPSIS

#example here
use GPI::Parser::GeneWise;

=head1 DESCRIPTION

It uses BioPerl modules to parse different output results.
(Bio::Tools::Genewise)

=head1 VERSIONS

$Id: GeneWise.pm,v 1.7 2007/03/01 09:28:03 equevill Exp $
$Log: GeneWise.pm,v $
Revision 1.7  2007/03/01 09:28:03  equevill
Added CVS log in the header of the file


Copyright (c) INRA/URGI 2005

=head1 AUTHORS / ACKNOWLEDGEMENTS

Emmanuel Quevillon <emmanuel.quevillon@versailles.inra.fr>

=cut

### FIXME : This module is not yet installed in the perl lib, so check were we are. ###

my $lib;

package GPI::Parser::GeneWise;

use strict;
use English;
use Data::Dumper;
use lib '/apps/GnpAnnot/lib';
use vars qw(@ISA);
use GPI::Parser;
@ISA = qw( GPI::Parser );
use GPI::GFF;
use Bio::Tools::Genewise;



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

    $self->usebioperl(1);

    return $self;

}



=head1 parse

Description: Parse (t)Blast(n,x) results using Bio::SearchIO
             object.
  Arguments: $in Blast results file.
    Returns: 0, $features on success
             1, message on error

=cut

sub parse($$){

    my ($self, $ifile) = @_;

    unless($self->_checkApplication($self->appl())){
	my $list = join(', ', sort keys %${GPI::ApplInfo2GFFFeat});
	$self->_exitOnError("Application name not set or not supported!: \n\tSupported: $list\n");
    }

    my $writer;
    my $features;
    my $bioperl = $self->usebioperl();

    if($bioperl){

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

    my $seqio = GPI::Bio::Tools::Genewise->new(-file => "$ifile");

    unless($seqio){
	unlink($self->out()) if($self->out());
	close($self->fh()) if($self->isFileHandle($self->fh()));
	$self->_exitOnError("Could not create GPI::Bio::Tools::Genewise object.");
    }

    local $OUTPUT_AUTOFLUSH = 1;

    my $maxhits = $self->maxhits();
    my($nresults) = (0);
    my $nbtranscript = 1;

    while(my $prediction = $seqio->next_prediction()){

	my ($progress, $next_update, $pbhits) = (0, 0, 0);

	$nresults++;
	my $nbtranscripts = scalar($prediction->transcripts());
	my $protid        = $seqio->_prot_id();
	my $target_id     = $seqio->_target_id();
	my $score         = $seqio->_score();
	my $gstart        = $prediction->start();
	my $gend          = $prediction->end();

	my @transcripts = $prediction->transcripts();

	foreach my $transcript (@transcripts){

	    if($self->verbose()){
		$progress = Term::ProgressBar->new({
						    count => $nbtranscripts,
						    name  => "Reading result $nresults [$protid]"
						   });
	    }

	    my @exons = $transcript->exons();
	    foreach my $exon (@exons){
		if($self->debug()){
		    print STDERR "------------------\n";
		    print STDERR "Start: " . $exon->location->start() . " End: " . $exon->location->end() . " Length: " . $exon->location->length(). "\n";
		    print STDERR "Location Type : " . $exon->location->location_type() . " FT: " . $exon->location->to_FTstring() . "\n";
		    print STDERR "------------------\n";
		}
	    }

#	    my @sorted = sort { $a->start() <=> $b->start() } (@polyAsites, @prime3utrs, @prime5utrs, @introns, @exons);
	    my @sorted = sort { $a->start() <=> $b->start() } (@exons);
	    my $item = 1;

	    #Depending on the strand, start and end do not behave the same.
	    my $strand = $transcript->strand == '1' ? '+' : '-';
	    my $start  = $strand eq '+' ? $gstart + 1 : $gstart + 2;
	    my $end    = $strand eq '+' ? $gend : $gend + 1;

	    my $setseqfg = Bio::SeqFeature::Generic->new(
							 -seq_id      => $target_id,
							 -source      => $self->appl(),
							 -primary_tag => $transcript->primary_tag(),
							 -start       => $start,
							 -end         => $end,
							 -score       => $score,
							 -strand      => $transcript->strand == '1' ? '+' : '-',
							 -frame       => '.',
							);
	    my $SetID = join("_", $self->appl(), $protid, $target_id, $transcript->primary_tag().sprintf('%0' . $NAMELEN . 'd', $nbtranscript));
	    $setseqfg->add_tag_value('ID', $SetID);
	    #Target is false as it is a DNA/Protein alignment
	    #$setseqfg->add_tag_value('Target', join('+', $protid, $sorted[0]->start(), $sorted[-1]->end()));
	    push(@$features, $setseqfg);

	    my $types = { };

	    foreach my $part (@sorted){

		if($self->debug()){
		    print STDERR "Type: " . $part->primary_tag() . "\n";
		    print STDERR "Start: " . $part->start() . " End: " . $part->end() . " Length: " . $part->length(). "\n";
		    print STDERR "------------------\n";		
		}

	       $types->{$part->primary_tag()}++;

		my $tags = { };
		my $ID   = join("_",  $SetID, $part->primary_tag() . sprintf('%0' . $NAMELEN . 'd', $types->{$part->primary_tag()}));
		
		#Depending on the strand, start and end positions need to be recalculated
		my $partstart = my $partend = 0;
		if($strand eq '+'){
		    $partstart = $part->start() + $gstart + 1;
		    $partend   = $part->end()   + $gstart;
		}else{
		    $partstart = $gend - $part->end()   + 2;
		    $partend   = $gend - $part->start() + 1;
		}

		my $seqfg = Bio::SeqFeature::Generic->new(
							  -seq_id        => $target_id,
							  -source        => $self->appl(),
							  -primary_tag   => $part->primary_tag(),  #Can be exon, intron, ....
							  -start         => $partstart,            #To be mapped on the ref sequence
							  -end           => $partend,              #To be mapped on the ref sequence
							  -score         => '.',
							  -strand        => $strand,
							  -frame         => $part->frame() || '.'
							 );

		$seqfg->add_tag_value('ID', $ID);
		$seqfg->add_tag_value('Name', $protid);
		#$seqfg->add_tag_value('Note', $Note);

		unless($self->flush()){
		    push(@$features, $seqfg);
		}else{
		    print $writer->write_feature($seqfg);
		}
		$item++;
	    }

	    $progress->update($nbtranscript) if($self->verbose());

	    $nbtranscript++;
	}
	
    }

    return($features);

}
