=head1 NAME

Sim4 - Module to parse Sim4 results.

=head1 SYNOPSIS

#example here
use GPI::Parser::Sim4;

=head1 DESCRIPTION

It uses BioPerl modules to parse different output results.
(Bio::SearchIO::sim4)

=head1 VERSIONS

$Id: Sim4.pm,v 1.35 2007/03/01 09:28:04 equevill Exp $
$Log: Sim4.pm,v $
Revision 1.35  2007/03/01 09:28:04  equevill
Added CVS log in the header of the file


Copyright (c) INRA/URGI 2005

=head1 AUTHORS / ACKNOWLEDGEMENTS

Emmanuel Quevillon <emmanuel.quevillon@versailles.inra.fr>

=cut


package GPI::Parser::Sim4;
use lib '/apps/GnpAnnot/lib';

use strict;
use English;
use Data::Dumper;
use vars qw(@ISA);
use GPI::Parser;
@ISA = qw( GPI::Parser );


=head1 new

Descrption: Create a new GPI::Parser::Sim4 object.
 Arguments: $params reference to hash table with parameters.
   Returns: Reference to GPI::Parser::Sim4 object or undef.

=cut

sub new ($;$$){

    my($class, @params) = @_;

    my $self = bless { }, $class;

    my %params =  $self->_init(\@params);

    $self->setParams(\%params);

    if($self->debug()){
	for my $key (sort keys %{$self->getParams()}){
	    $self->_print_debug("$key => $params{$key}");
	}
    }

    return $self;

}



=head1 parse

Description: Parse Sim4 results using Bio::SearchIO
             object.
             It returns a hash table containing each result,
             hits and hsp in a specific structure.
  Arguments: $in sim4 results file.
    Returns: 0, $features on success
             1, message on error

=cut

sub parse ($){

    require Bio::SearchIO;

    my ($self, $ifile) = @_;

    unless($self->_checkApplication($self->appl())){
	my $list = join(', ', sort keys %${GPI::ApplInfo2GFFFeat});
	$self->_exitOnError("Application name not set or not supported!: \n\tSupported: $list\n");
    }

    my $writer;
    my $features;
    my $bioperl = $self->usebioperl();

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

    my $seqio = Bio::SearchIO->new(-file => "$ifile",
				   -format => 'sim4',
				  );

    unless($seqio){
	unlink($self->out()) if($self->out());
	close($self->fh()) if($self->isFileHandle($self->fh()));
	$self->_exitOnError("Could not create Bio::SearchIO object.");
    }

    local $OUTPUT_AUTOFLUSH = 1;

    my $maxhits = $self->maxhits();
    my($nresults, $nhits, $nhsps) = (0, 0, 0);

    while(my $result = $seqio->next_result()){

	my ($progress, $next_update, $pbhits) = (0, 0, 0);

	$nresults++;

	$self->printDebug($result) if($self->debug());

	if($self->verbose()){
	    $progress = Term::ProgressBar->new({count => $result->num_hits(),
						name  => "Reading result $nresults [" . $result->query_name() . "]"}
					      );
	}

	$self->_print_verbose("\tReading results $nresults [", $result->query_name(), "] ");
	#Info about the query
	my $query_acc  = $result->query_accession();
	my $query_desc = $result->query_description();
	my $query_len  = $result->query_length();
	my $query_name = $result->query_name();

	if($bioperl){

	    my $RefTag = { };
	    $RefTag->{ID} = $query_name || $self->getApplInfo('region') . sprintf('%0'.$GPI::Parser::NAMELEN.'d', $nresults);

	    my $RefFeat =  Bio::SeqFeature::Generic->new(
							 -seq_id        => $query_name || "Region_$nresults",
							 -source        => $self->appl() || 'Generic_source',
							 -primary_tag   => $self->qtype() || $self->getApplInfo('region'),
							 -start         => 1,
							 -end           => $query_len,
							 -score         => '.',
							 -strand        => '+',
							 -frame         => '.',
							 -tag           => $RefTag
							);

	    unless($self->flush()){
		push(@$features, $RefFeat);
	    }else{
		print $writer->write_feature($RefFeat);
	    }
		
	}
	
	while(my $hit = $result->next_hit()){

	    #Check if there is/are some result(s) for this hit. Right place for calling this test otherwise lots of WARNING printed.

	    $pbhits++;
	    $next_update = $progress->update($pbhits) if($pbhits >= $next_update && $self->verbose());

	    next unless($hit->num_hsps() > 0);

	    my @PartFeat = ();

	    $nhits++;

	    $self->printDebug($hit) if($self->debug());

	    my $hits = { };
	    my $msg;

	    if($self->haveFilters()){
		$hits->{sequence} = $self->initArraySeq($hit->length());
		$self->_exitOnError($msg) unless(defined($hits->{sequence}));
	    }

	    $hits->{name}      = $hit->name();
	    $hits->{desc}      = $hit->description();
	    $hits->{qacc}      = $query_acc;
	    $hits->{qdesc}     = $query_desc;
	    $hits->{qname}     = $query_name;
	    #It causes Bio::Tools::GFF to crash
	    $hits->{desc}    =~ s/[;,]+/ /g;
	    $hits->{length}    = $hit->length();
	    $hits->{length_aln}=$hit->length_aln();
	    if($hits->{length_aln} > $hits->{length} && $self->debug()){
		print STDERR "******************************************************************\n";
		print STDERR "You might have some problems (overlap/hsp matching somewhere else)\n";
		print STDERR "with hit(target_id) ".$hit->name()." and query $query_name\n";
		print STDERR "Query length (".$hits->{length}.") < length_aln(".$hits->{length_aln}.")\n";
		print STDERR "Have a look at the end of the parsing or into the Blast report.\n";
		print STDERR "******************************************************************\n";
	    }
	    $hits->{algo}      = $hit->algorithm();
	    $hits->{score}     = defined($hit->raw_score()) || "";
	    $hits->{bits}      = $hit->bits();
	    $hits->{evalue}    = defined($hit->significance()) || "";
	    $hits->{nbhits}    = $hit->n();
	    $hits->{rank}      = $hit->rank();
	    $hits->{strandq}   = $hit->strand('query');
	    $hits->{strands}   = $hit->strand('subject');
	    $hits->{strand}    = $hit->strand('subject');
	    $hits->{matchid}   = $hit->matches('id');
	    $hits->{matchcon}  = $hit->matches('cons');
	    $hits->{frame}     = $hit->frame();
	    $hits->{eachacc}   = $hit->each_accession_number();
	    $hits->{id}        = $nhits;
	    #Get the start of the first hsp and the end of the last one to get the coordinates of the coverage
	    $hits->{hspsstart} = (sort { $a <=> $b } map { $_->start('subject') } $hit->hsps())[0] || '';
	    $hits->{hspsend}   = (sort { $a <=> $b } map { $_->end('subject') } $hit->hsps())[-1] || '';

	    #We fill the sequence array in case of certain filters are requested
	    if($self->haveFilters()){
		for($hit->seq_inds('hit', 'identical')){
		    $hits->{sequence}->[$_-1] = 3;
		}

		for($hit->seq_inds('hit', 'conserved')){
		    $hits->{sequence}->[$_-1] = 2 unless($hits->{sequence}->[$_-1] > 1);
		}
	    }

	    while(my $hsp = $hit->next_hsp()){

		$self->printDebug($hsp) if($self->debug());

		$nhsps++;

		unless($bioperl){

		    my $hash = { };

		    $hash->{evalue}   = $hsp->evalue();
		    $hash->{gapst}    = $hsp->gaps('total');
		    $hash->{gapsq}    = $hsp->gaps('query');
		    $hash->{gapss}    = $hsp->gaps('hit');

		    $hash->{lent}     = $hsp->length('total');
		    $hash->{lenq}     = $hsp->length('query');
		    $hash->{lens}     = $hsp->length('hit');

		    #Specific case depending on the strand of the subject
		    $hash->{strandq}  = $hsp->strand('query');
		    $hash->{strands}  = $hsp->strand('subject');
		    $hash->{strand}   = $hash->{strands} =~ /-/ ? $hash->{strands} : $hash->{strandq};

		    $hash->{startq}   = $hsp->start('query');
		    $hash->{starts}   = $hsp->strand('subject') =~ /-/ ? $hit->length() - $hsp->end('subject') + 1 : $hsp->start('subject');
		    $hash->{start}    = $hsp->start('query');

		    $hash->{endq}     = $hsp->end('query');
		    $hash->{ends}     = $hsp->strand('subject') =~ /-/ ? $hit->length() - $hsp->start('subject') + 1 : $hsp->end('subject');
		    $hash->{end}      = $hsp->end('query');

		    #Real coordinates without changes because of the strand
		    $hash->{rstarts}    = $hsp->start('subject');
		    $hash->{rends}      = $hsp->end('subject');

		    #print STDERR "HSP: Start $hash->{starts}($hash->{rstarts}) End $hash->{ends}($hash->{rends})\n";

		    $hash->{indidss}  = ($hsp->seq_inds('subject', 'identical'))[0];
		    $hash->{indidse}  = ($hsp->seq_inds('subject', 'identical'))[-1];
		    $hash->{indconss} = ($hsp->seq_inds('subject', 'conserved'))[0];
		    $hash->{indconse} = ($hsp->seq_inds('subject', 'conserved'))[-1];

		    $hash->{indidqs}  = ($hsp->seq_inds('query', 'identical'))[0];
		    $hash->{indidqe}  = ($hsp->seq_inds('query', 'identical'))[-1];
		    $hash->{indconqs} = ($hsp->seq_inds('query', 'conserved'))[0];
		    $hash->{indconqe} = ($hsp->seq_inds('query', 'conserved'))[-1];

		    $hash->{rank}     = $hsp->rank();
		    $hash->{matchesq} = $hsp->matches('query');
		    $hash->{matchess} = $hsp->matches('hit');
		    $hash->{rangeq}   = $hsp->range('query');
		    $hash->{ranges}   = $hsp->range('hit');
		    $hash->{id}       = $nhsps;
		    $hash->{parent}   = $nhits;
		    $hash->{name}     = $hits->{name};

		    $hits->{hsps} = [ ] unless(defined($hits->{hsps}));

		    push(@{$hits->{hsps}}, $hash);
		
		}
		elsif($bioperl){
		    unless(defined($self->fh())){
			$self->_exitOnError("No file handle defined for writing.");
		    }

		    my $PartTag = { };
		    my $PartID = $self->annot() ."_" if($self->annot());
		    my $PartID = join("_",
				      $self->appl(),
				      $query_name ? $query_name : "Region_$nresults",
				      $hit->name() ,
				      'match_set' . sprintf('%0' . $GPI::Parser::NAMELEN . 'd', $nhits),
				     );

		    $PartTag->{ID}     = $PartID . '_match_part' . sprintf('%0' . $GPI::Parser::NAMELEN . 'd', $nhsps);
		    $PartTag->{Parent} = $PartID;
		    $PartTag->{Target} = join('+', $hit->name(), $hsp->start('subject'), $hsp->end('subject'));
		    #To be consitstent with E.Gicquelo
		    #$PartTag->{Note}   = "DB:" . $result->database_name() . ",";
		    #$PartTag->{Note}  .= "Align:" . $hsp->length() . ",";
		    #$PartTag->{Note}  .= "Identity:" . $hsp->percent_identity() . ",";
		    #$PartTag->{Note}  .= "Desc:" . $hit->name() . "_" . $hit->description();

		    $PartTag->{target_pcover} = $hsp->length('subject');
		    $PartTag->{target_pcons}  = $hsp->num_conserved();
		    $PartTag->{target_pident} = $hsp->num_identical();
		    $PartTag->{Note}         .= "Gaps:=".$hsp->gaps() if($hsp->gaps());
		    $PartTag->{Mstring}       = $hsp->homology_string();
		    $PartTag->{Starts}        = $hsp->start('subject');
		    $PartTag->{Ends}          = $hsp->end('subject');
		    $PartTag->{target_start}  = $hsp->start('query');
		    $PartTag->{target_end}    = $hsp->end('query');




		    my $Partfeature = Bio::SeqFeature::Generic->new(
								    -seq_id       => $query_name ? $query_name : "Region_$nresults",
								    -source_tag   => $self->appl() || 'Generic_source',
								    -primary_tag  => 'match_part',
								    -start        => $hsp->start('query'),
								    -end          => $hsp->end('query'),
								    -score        => $hsp->evalue(),
								    -strand       => $hsp->strand('subject') =~ /-/ ? $hsp->strand('subject') : $hsp->strand('query'),
								    -frame        => $hit->{frame} ? $hit->{frame} : '.',
								    -tag          => $PartTag
								   );
		    push(@PartFeat, $Partfeature);

		}else{
		    $self->_exitOnError("Problem while parsing results.");
		}
	    } #end next_hsp 

	    $nhsps = 0;

	    #This won't be true if we use bioperl because hits->{hsps} won't never be fill.
	    if(defined($hits->{hsps})){

		if($self->haveFilters()){
		    ($hits->{pcoverage}, $hits->{pidentity}, $hits->{nconserved}) = $self->getStats($hits);
		}

		if($self->debug()){
		    #print STDERR "Type\tCoverage\tIdentity\tConserved\n";
		    #print STDERR "----\t--------\t--------\t---------\n";
		    #print STDERR "\t" . join("\t", $hits->{pcoverage}, $hits->{pidentity}, $hits->{nconserved}). "\n";
		    #print STDERR "======================================================\n\n";
		}

		if($self->applyFilters($hits)){
		    $features->{$nresults} = { } unless($features->{$nresults});;
		    $features->{$nresults}->{topdesc}->{title} = "##sequence-region";
		    $features->{$nresults}->{topdesc}->{name}  = $query_name ? $query_name : "Region_$nresults";
		    $features->{$nresults}->{topdesc}->{name}  =~ s/.+\/([^\/]+)$/$1/;
		    $features->{$nresults}->{topdesc}->{start} = 1;
		    $features->{$nresults}->{topdesc}->{end}   = $query_len;
		    $features->{$nresults}->{hits}->{$nhits}   = $hits;

		    if($self->flush() && !$bioperl){
			my ($res, $msg) = $self->writeGFF3($features, $self->fh());
			$self->_topSeen(1);
			$self->_exitOnError($msg) if($res);
			$features = { };
		    }

		}

	    } #end if(defined($hits->{hsps}))

	    elsif($bioperl){

		#Here we need to store in the right order, results, hits and hsps to be right written in the gff output.
		my @q = @PartFeat;
		my @s = @PartFeat;

		@q = sort { $a->start() <=> $b->start() } @PartFeat;
		@s = map { ($_->get_tag_values('Target'))[0] } @PartFeat;
		
		my $SetStart = $q[0]->start();
		my $SetStop  = $q[-1]->end();
		
		my (@s2, @s3);
		
		for(@s){
		    push(@s2, (split(/\+/, $_))[1]);
		    push(@s3, (split(/\+/, $_))[2]);
		}
		
		@s = sort { $a <=> $b } @s2;
		my $FeatStart = $s[0];
		
		@s = sort { $a <=> $b } @s3;
		my $FeatStop  = $s[-1];
		
		my $SetStrand = $hit->strand('query');
		$SetStrand =~ s/\d+//;
		$SetStrand = $SetStrand eq '' ? '+' : '-';
		my $SetFrame = $hit->frame();
		$SetFrame =~ s/[-\+]//;
		my $SetNum = sprintf('%0'.$GPI::Parser::NAMELEN.'d', $nhits);
		
		my $SetTag = { };
		
		my $SetID;
		$SetID .= $self->annot() ."_" if($self->annot());
		$SetID .= join("_", $self->appl(), $query_name || "Region_$nresults");
		$SetID .= "_" . $hit->name() if($hit->name());
		$SetID .= "_" . join("_", $self->getApplInfo('set') . $SetNum);
		
		my $SetName = join("_", $query_name, $hit->name());
		
		my $SetTarget = join('+', $hit->name(), $FeatStart, $FeatStop);
		
		$SetTag->{ID}     = $SetID;
		$SetTag->{Name}   = $SetName;
		$SetTag->{Target} = $SetTarget;
		my ($gaps);

		for(@PartFeat){

		    $gaps   += ($_->get_tag_values('Gapst'))[0];

		    $_->remove_tag('Starts')  if($_->has_tag('Starts'));
		    $_->remove_tag('Startq')  if($_->has_tag('Startq'));
		    $_->remove_tag('Ends')    if($_->has_tag('Ends'));
		    $_->remove_tag('Endq')    if($_->has_tag('Endq'));
		    $_->remove_tag('Gapst')   if($_->has_tag('Gapst'));
		    $_->remove_tag('Pcover')  if($_->has_tag('Pcover'));
		    $_->remove_tag('Pident')  if($_->has_tag('Pident'));
		    $_->remove_tag('Mstring') if($_->has_tag('Mstring'));
		    $_->remove_tag('Ncons')   if($_->has_tag('Ncons'));
		}

		my ($pcoverage, $pidentity, $nconserved) = $self->getStats($hits,\@PartFeat);

		if($self->addnote() || defined($self->pcoverage()) || defined($self->pidentity()) || defined($self->nconserved())){

		    $SetTag->{target_pcover} = sprintf("%.2f", $pcoverage);
		    $SetTag->{target_pident} = sprintf("%.2f", $pidentity);
		    $SetTag->{target_pcons}  = sprintf("%.2f", $nconserved);
		    $SetTag->{target_length} = $hit->length() if($hit->length());
		    $SetTag->{target_start}  = $FeatStart;
		    $SetTag->{target_end}    = $FeatStop;
		    $SetTag->{target_id}     = $hit->name();
		    $SetTag->{Note}         .= ",Gaps:$gaps" if($gaps);
		    $SetTag->{lib}           = $self->annot() if ($self->annot());
		    $SetTag->{program}       = $self->appl() if ($self->appl());
		}

		my $sourcetag=$self->appl() || "Generic_source" ;
		if($self->annot()){$sourcetag .= "_" . $self->annot();}

		my $SetFeat = Bio::SeqFeature::Generic->new(
							    -seq_id      => $query_name || "Results_$nresults",
							    -source_tag  => $sourcetag,
							    -primary_tag => $self->getApplInfo('set'),
							    -start       => $SetStart,
							    -end         => $SetStop,
							    -score       => $hit->significance(),
							    -strand      => $SetStrand || '.',
							    -frame       => $SetFrame || '.',
							    -tag         => $SetTag
							   );
		
		next unless($self->applyFilters($SetFeat));

		unless($self->flush()){
		    push(@$features, $SetFeat, @PartFeat);
		}else{
		    print $writer->write_feature($SetFeat, @PartFeat);
		}
		
	    } #end elsif($bioperl)
	    else{
		$self->_exitOnError("Problem during paring, we don't know if we need to use bioperl or not....");
	    }
	    #$nhsps = 0;

	    last if($maxhits && $nhits >= $maxhits);

	} #end next_hits()
	
	#Progress Bar, to complete to 100% the progress bar
	$progress->update($result->num_hits()) if($self->verbose());# if($result->num_hits >= $next_update);

	$nhits = 0;

    } #end next_result()

    $features->{nresults} = $nresults - 1 unless($self->usebioperl());

    return($features);
}




1;
