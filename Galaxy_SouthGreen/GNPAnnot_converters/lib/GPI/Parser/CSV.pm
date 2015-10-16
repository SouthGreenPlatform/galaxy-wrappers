=head1 NAME

CSV - Module to read csv file and create GFF3 output.

=head1 SYNOPSIS

#example here
use GPI::Parser::CSV;

=head1 DESCRIPTION

This module was firstly created to read csv format file from Broad Institute and create GFF3 output.

=head1 VERSIONS

$Id: CSV.pm,v 1.4 2007/03/01 09:13:36 equevill Exp $
$Log: CSV.pm,v $
Revision 1.4  2007/03/01 09:13:36  equevill
CVS log added into the header of the file


Copyright (c) INRA/URGI 2005

=head1 AUTHORS / ACKNOWLEDGEMENTS

Emmanuel Quevillon <emmanuel.quevillon@versailles.inra.fr>

=cut

package GPI::Parser::CSV;
use lib '/apps/GnpAnnot/lib';
use strict;
use Data::Dumper;
use vars qw(@ISA);
use GPI::Parser;
@ISA = qw( GPI::Parser );

our $CSVLIST = {
		'refname' => 'The reference name the exon/gene is mapped to.',
		'estart ' => 'The start position of the exon.',
		'estop  ' => 'The stop position of the exon.',
		'estrand' => 'The strand of the exon.',
		'gstart ' => 'The start position of the gene = The start of the first exon somehow.',
		'gstop  ' => 'The stop position of the gene = The stop of the last exon somehow.',
		'gstrand' => 'Strand of the gene on the reference.',
		'gname  ' => 'Name of the gene.(Optional)',
		'glength' => 'Length of the gene.(Optional)',
		'gid    ' => 'ID of the gene.',
		'eid    ' => 'ID of the exon.(Optional).',
		'score  ' => 'A score to be associated with',
		'phase  ' => 'Phase of the exon',
		'source ' => 'Source of the annotation',
		'type   ' => 'Type of the feature.(Optional. default=gene).',
	      };

my $ARGS = [qw(refname gstart gstop gid gstrand)];

=head1 new

 Description: Creates a new GPI::Parser::CSV object
   Arguments: csvlist : The order the elements are in the csv file.
     Returns: blessed reference

=cut

sub new {

    my($class, @params) = @_;

    my $self = bless { }, $class;

    my %params =  $self->_init(\@params);

    $self->setParams(\%params);

    if($self->getParam('debug')){
	for my $key (sort keys %{$self->getParams()}){
	    $self->_print_debug("$key => $params{$key}");
	}
    }

    $self->usebioperl(1);

    unless($self->csvlist()){
	$self->_exitOnError("You need to provide a list of field in order mapping you csv file.");
    }else{
	$self->_init_list();
    }

    my($ret, $msg) = $self->_check_list();
    $self->_exitOnError($msg) if($ret);

    return $self;

}

=head1 _init_list

 Description: Initialize and set into memory csv list
   Arguments:
    Returns : 

=cut

sub _init_list {

    my($self) = @_;

    $self->_print_debug("_init_list csv ... ");

    $self->{_csvlist} = { };
    my $i = 0;
    for(split(/,/, $self->csvlist())){
	$self->{_csvlist}->{$_} = $i;
	$self->_print_debug("\tOption: [$_] $i");
	$i++;
    }

    return;
}

=head1 _check_list

 Description: Check that spme required arguments have been passed on the command line
   Arguments:
     Returns: mess on error

=cut

sub _check_list {

    my($self) = @_;

    for(@$ARGS){
	unless(exists($self->{_csvlist}->{$_})){
	    $self->_exitOnError("Arguments $_ is required for cvslist.");
	}
    }

    return;
}

sub getkey {

    my($self, $key) = @_;

    if(defined($key)){
	return $self->{_csvlist}->{$key};
    }

    return undef;
}

=head1 parse

 Description: Parses a csv file
   Arguments: file - csv input file
     Returns: hash reference

=cut

sub parse ($){

    my($self, $ifile) = @_;

    unless($ifile){
	$self->_exitOnError("No input file given to parse");
    }

    unless($self->_checkApplication($self->appl())){
	my $list = join(', ', sort keys %${GPI::ApplInfo2GFFFeat});
	$self->_exitOnError("Application name not set or not supported!: \n\tSupported: $list\n");
    }

    my $writer;
    my $features = [ ];
    my $bioperl = $self->usebioperl();

    require Bio::SeqFeature::Generic;
    use GPI::GFF;

    if($self->flush()){
	$writer = GPI::GFF->new(
					    -gff_version => 3,
					    -fh          => $self->fh() || \*STDOUT
					   );
    }

    my $csvfh;
    unless(open($csvfh, "<$ifile")){
	$self->_exitOnError("Could not open $ifile : $!");
    }

    my $seen_genes = { };
    my $genes = 0;
    my $exons = 0;
    while(defined(my $line = <$csvfh>)){

	next if $line =~ /^#/;
	chomp $line;
	my(@fields) =  split(/,/, $line);

	s/"//g for(@fields);

	#We are in presence of exons, so we expect that we have some genes as well, e.g.(Broad):
	#EXON_START,EXON_STOP,EXON_STRAND,FULL_GENE_NAME,GENE_START,GENE_STOP,GENE_LENGTH,GENE_STRAND,CONTIG_NAME,GENE_LOCUS,GENE_LOCUS_VERSION
	if($self->getkey('estart') && $self->getkey('estop') && $self->getkey('gstart') && $self->getkey('gstop')){

	    #We already seen the gene with some exons?
	    unless(exists($seen_genes->{$fields[$self->getkey('gid')] || $fields[$self->getkey('gname')]})){

		################
		## GENE INFOS ##
		################
		$genes++;
		my $tag = {
			   'ID' => join("_",
					$self->annot ? join("_", $self->appl(), $self->annot()) : $self->appl(),
					$fields[$self->getkey('refname')],
					$fields[$self->getkey('gid')] || $fields[$self->getkey('gname')] || "GENE".sprintf('%0'.$NAMELEN.'d', $genes),
					$self->getApplInfo('set') . sprintf('%0'.$NAMELEN.'d', $genes)
				       ),
			  };
		my $feat = Bio::SeqFeature::Generic->new(
							 -seq_id      => $self->getkey('refname') ? $fields[$self->getkey('refname')] : '.',
							 -start       => $self->getkey('gstart')  ? $fields[$self->getkey('gstart')]  : '.',
							 -end         => $self->getkey('gstop')   ? $fields[$self->getkey('gstop')]   : '.',
							 -source_tag  => $self->getkey('source')  ? $fields[$self->getkey('source')]  : $self->appl(),
							 -primary_tag => $self->qtype() || $self->getApplInfo('set'),
							 -frame       => $self->getkey('frame')   ? $fields[$self->getkey('frame')]   : '.',
							 -score       => $self->getkey('score')   ? $fields[$self->getkey('score')]   : '.',
							 -strand      => $self->getkey('strand')  ? $fields[$self->getkey('gstrand')]  : '.',
							 -tag         => $tag,
							);
		#We keep trace of the gene having multiple exons
		$seen_genes->{$fields[$self->getkey('gid') || $self->getkey('gname')]} = 1;

		if($self->flush()){
		    $writer->write_feature($feat);
		}
		else{
		    push(@$features, $feat);
		}
	    } #end unless

	    #################
	    ## EXONS INFOS ##
	    #################
	    $exons++;
	    my $tag = {
		       'ID' => join("_",
				    $self->annot ? join("_", $self->appl(), $self->annot()) : $self->appl(),
				    $fields[$self->getkey('refname')],
				    $fields[$self->getkey('gid')] || $fields[$self->getkey('gname')] || "GENE".sprintf('%0'.$NAMELEN.'d', $genes),
				    $self->getApplInfo('set') . sprintf('%0'.$NAMELEN.'d', $genes),
				    $self->getApplInfo('part') . sprintf('%0'.$NAMELEN.'d', $exons),
				       ),
			  };
	    my $feat = Bio::SeqFeature::Generic->new(
						     -seq_id      => $self->getkey('refname') ? $fields[$self->getkey('refname')] : '.',
						     -start       => $self->getkey('estart')  ? $fields[$self->getkey('estart')]  : '.',
						     -end         => $self->getkey('estop')   ? $fields[$self->getkey('estop')]   : '.',
						     -source_tag  => $self->getkey('source')  ? $fields[$self->getkey('source')]  : $self->appl(),
						     -primary_tag => $self->qtype() || $self->getApplInfo('part'),
						     -frame       => $self->getkey('frame')   ? $fields[$self->getkey('frame')]   : '.',
						     -score       => $self->getkey('score')   ? $fields[$self->getkey('score')]   : '.',
						     -strand      => $self->getkey('strand')  ? $fields[$self->getkey('estrand')]  : '.',
						     -tag         => $tag,
						    );

	    if($self->flush()){
		$writer->write_feature($feat);
	    }
	    else{
		push(@$features, $feat);
	    }
	}
	else{
	    $genes++;
	    my $tag = {
		       'ID' => join("_",
				    $self->annot ? join("_", $self->appl(), $self->annot()) : $self->appl(),
				    $fields[$self->getkey('refname')],
				    $fields[$self->getkey('gid')] || $fields[$self->getkey('gname')] || "GENE".sprintf('%0'.$NAMELEN.'d', $genes),
				    $self->getApplInfo('set') . sprintf('%0'.$NAMELEN.'d', $genes)
				   ),
		      };
	    #We create the 'gene' feature
	    my $feat = Bio::SeqFeature::Generic->new(
						     -seq_id      => $self->getkey('refname') ? $fields[$self->getkey('refname')] : '.',
						     -start       => $self->getkey('gstart')  ? $fields[$self->getkey('gstart')]  : '.',
						     -end         => $self->getkey('gstop')   ? $fields[$self->getkey('gstop')]   : '.',
						     -source_tag  => $self->getkey('source')  ? $fields[$self->getkey('source')]  : $self->appl(),
						     -primary_tag => $self->qtype() || $self->getApplInfo('set'),
						     -frame       => $self->getkey('frame')   ? $fields[$self->getkey('frame')]   : '.',
						     -score       => $self->getkey('score')   ? $fields[$self->getkey('score')]   : '.',
						     -strand      => $self->getkey('strand')  ? $fields[$self->getkey('gstrand')]  : '.',
						     -tag         => $tag,
						    );

	    if($self->flush()){
		$writer->write_feature($feat);
	    }
	    else{
		push(@$features, $feat);
	    }
	}
    }

    return($features);

}
