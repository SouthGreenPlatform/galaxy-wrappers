=head1 NAME

RepeatMasker - Module to parse RepeatMasker results.

=head1 SYNOPSIS

#example here
use GPI::Parser::RepeatMasker;

=head1 DESCRIPTION

The parsing is made line by line and store the results (if requiered) in a Bio::SeqFeature::Generic object.

=head1 VERSIONS

$Id: RepeatMasker.pm,v 1.20 2007/03/01 09:28:04 equevill Exp $
$Log: RepeatMasker.pm,v $
Revision 1.20  2007/03/01 09:28:04  equevill
Added CVS log in the header of the file


Copyright (c) INRA/URGI 2005

=head1 AUTHORS / ACKNOWLEDGEMENTS

Emmanuel Quevillon <emmanuel.quevillon@versailles.inra.fr>

=cut


package GPI::Parser::RepeatMasker;
use lib '/apps/GnpAnnot/lib';
use strict;
use English;
use Data::Dumper;
use vars qw(@ISA);
use GPI::Parser;
@ISA = qw( GPI::Parser );


=head1 new

Descrption: Create a new GPI::Parser::RepeatMasker object.
 Arguments: $params reference to hash table with parameters.
   Returns: Reference to GPI::Parser::Blast object or undef.

=cut

sub new (@){

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

Description: Parse RepeatMasker results using Bio::SearchIO
             object.
             It returns a hash table containing each result,
             hits and hsp in a specific structure.
  Arguments: $in sim4 results file.
    Returns: 0, $features on success
             1, message on error

=cut

sub parse ($$){

    require Bio::SeqFeature::Generic;

    my ($self, $ifile) = @_;

    unless($self->_checkApplication($self->appl())){
	my $list = join(', ', sort keys %${GPI::ApplInfo2GFFFeat});
	$self->_exitOnError("Application name not set or not supported!: \n\tSupported: $list\n");
    }

    #$self->_exitOnError("Fitlers not supported for RepeatMasker.") if($self->haveFilters());

    my $writer;
    my $features;
    my $bioperl = $self->usebioperl();

    if($bioperl){
	#require Bio::SeqFeature::Generic;
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

    unless(open(IN, "<$ifile")){
	$self->_exitOnError("Could not open $ifile for parsing : $!");
    }

    local $OUTPUT_AUTOFLUSH = 1;

    my($nresults, $nhits, $nhsps) = (0, 0, 0);
    my $backName;
    my $backLength;

    while (<IN>) {

	chomp;

	if (/no repetitive sequences detected/) {
	    print STDERR "RepeatMasker didn't find any repetitive sequences\n";
	    return (0, $features);
	}

	if (/\d+/) { ##ignore introductory lines
	    my @element = split;
	    ## ignore features with negatives
	    next if ($element[11-13] =~ /-/);

	    ### TODO recuperer 2,3,4 colonnes.
	    ### DONE but what to do with it.

	    my ($score, $percdiv, $percdel, $percins, $query_name, $query_start, $query_end, $rest, $strand, $repeat_name, $repeat_class) = (split)[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
	    $rest =~ s/\(|\)//g;

	    my $ctg_length = ($rest + $query_end);

	    if($backLength != $ctg_length && $self->verbose()){
		$backLength = $ctg_length;
	    }

	    my ($hit_start,$hit_end, $hit_length);
	    if ($strand eq '+') {
		($hit_start, $hit_end, $hit_length) = (split)[11, 12, 13];
		$strand = 1;
	    }
	    elsif ($strand eq 'C') {
		($hit_length, $hit_end, $hit_start) = (split)[11, 12, 13];
		$strand = -1;
	    }else{
		return (1, "Problem while parsing strand.");
	    }

	    $hit_length =~ s/\(|\)//g;
	    $hit_length += $hit_end;

	    #We now parse class and subclass from the repeat, so we don't next when simple_repeat or low_complexity
	    #next if($repeat_class eq 'Simple_repeat' || $repeat_class eq 'Low_complexity');
	    my ($rclass, $rsubclass) = split(/\//, $repeat_class);

	    if($backName ne $query_name){
		print STDERR "\tReading results $nresults [", $query_name, "] " if($self->verbose());
		print STDERR "." if($self->verbose());
		$nhits = 0;
		$nresults++;
		$backName = $query_name;
	    }

	    $nhits++;

	    if($bioperl){

		#To avoid repetitive definition of the query
		my $rf = undef;

		if($backName ne $query_name){

		    $rf = Bio::SeqFeature::Generic->new(
							-seq_id      => $query_name,
							-source_tag  => $self->appl()  || 'Generic_source',
							-primary_tag => $self->qtype() || $self->getApplInfo('region'),
							-start       => 1,
							-end         => $ctg_length,
							-score       => '.',
							-strand      => $strand,
							-frame       => '.',
						       );
		}

		my $rf2 = Bio::SeqFeature::Generic->new(
							-seq_id      => $query_name,
							-source_tag  => $self->appl()  || 'Generic_source',
							-primary_tag => $self->qtype() || $self->getApplInfo('set'),
							-start       => $query_start,
							-end         => $query_end,
							-score       => $score,
							-strand      => $strand,
							-frame       => '.',
						       );
							

		#$self->printDebug('result', $result) if($self->debug());

		my $SetTag = { };
		$SetTag->{ID} = $query_name;

		my $PartTag = { };
		$PartTag->{ID} = join("_",
				      $self->appl(),
				      $query_name,
				      $repeat_name,
				      $self->getApplInfo('set') . sprintf('%0'.$NAMELEN.'d', $nhits)
				     );

		$PartTag->{Name}          = $rf2->seq_id();
		$PartTag->{target_length} = $hit_length if($self->addnote() && $hit_length);
		#$PartTag->{Note} = "HitLen:$hit_length";
		$PartTag->{Target}        = join("+",
						 $repeat_name,
						 $hit_start,
						 $hit_end
						);
		$PartTag->{repeat_class}    = $rclass;
		$PartTag->{repeat_subclass} = $rsubclass if($rsubclass);

		if($backName ne $query_name){
		    $rf->set_attributes(
					-tag         => $SetTag,
					-primary_tag => $self->getApplInfo('region')
				       );
		}
		$rf2->set_attributes(
				     -tag         => $PartTag,
				     -primary_tag => $self->getApplInfo('set')
				    );

		unless($self->flush()){
		    defined($rf) ? push(@$features, $rf, $rf2) : push(@$features, $rf2);
		}else{
		    defined($rf) ? print $writer->write_feature($rf, $rf2) : print $writer->write_feature($rf2);
		}

	    } #end if($bioperl)

	    else{

		my $hits = { };

		$hits->{name}    = $repeat_name;
		#$hits->{evalue}  = $score;
		$hits->{strandq} = '+';
		$hits->{strands} = $strand =~ /-/ ? '-' : '+';
		$hits->{strand}  = $strand;
		$hits->{starts}  = $hit_start;
		$hits->{startq}  = $query_start;
		$hits->{ends}    = $hit_end;
		$hits->{endq}    = $query_end;
		$hits->{length}  = $hit_length;
		$hits->{frame}   = '.';
		$hits->{id}      = $nhits;
		push(@{$hits->{_private_tags}}, "repeat_class=$rclass");
		push(@{$hits->{_private_tags}}, "repeat_subclass=$rsubclass") if($rsubclass);

		$features->{$nresults} = { } unless($features->{$nresults});;
		$features->{$nresults}->{topdesc}->{title} = "##sequence-region";
		$features->{$nresults}->{topdesc}->{name}  = $query_name ? $query_name : "Region_$nresults";
		$features->{$nresults}->{topdesc}->{name}  =~ s/.+\/([^\/]+)$/$1/;
		$features->{$nresults}->{topdesc}->{start} = 1;
		$features->{$nresults}->{topdesc}->{end}   = $ctg_length;
		$features->{$nresults}->{hits}->{$nhits}   = $hits;

	    }

	    if($self->getParam('backname') ne $query_name){ $self->setParam('backname', $query_name);}

	    if($self->flush() && !$bioperl){
		my ($res, $msg) = $self->writeGFF3($features, $self->fh());
		$self->_topSeen(1);
		$self->_exitOnError($msg) if($res);
		$features = { };
	    }	

	    print STDERR " DONE\n" if($self->verbose());

	} #end if(/(\d+)/)

    } #end while	

    $features->{nresults} = $nresults - 1 unless($self->usebioperl());
	
    return($features);;

}




1;
