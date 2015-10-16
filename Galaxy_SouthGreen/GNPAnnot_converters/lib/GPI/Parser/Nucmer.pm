=head1 NAME

Sim4 - Module to parse Nucmer .coords file results.

=head1 SYNOPSIS

#example here
use GPI::Parser::Nucmer;

=head1 DESCRIPTION

This modules is used to parse .coords file produces by nucmer with option
--coords or file that was produce with command 'show-coords' from a .delta
file from Mummer package.

=head1 VERSIONS

$Id: Nucmer.pm,v 1.3 2007/03/01 09:28:04 equevill Exp $
$Log: Nucmer.pm,v $
Revision 1.3  2007/03/01 09:28:04  equevill
Added CVS log in the header of the file


Copyright (c) INRA/URGI 2007

=head1 AUTHORS / ACKNOWLEDGEMENTS

Emmanuel Quevillon <emmanuel.quevillon@versailles.inra.fr>

=cut


package GPI::Parser::Nucmer;

use strict;
use English;
use Data::Dumper;
use vars qw(@ISA);
use GPI::Parser;
@ISA = qw( GPI::Parser );


=head1 new

Descrption: Create a new GPI::Parser::Nucmer object.
 Arguments: $params reference to hash table with parameters.
   Returns: Reference to GPI::Parser::Nucmer object or undef.

=cut

sub new ($;$$){

    my($class, @params) = @_;

    my $self = bless { }, $class;

    my %params =  $self->_init(\@params);

    $self->setParams(\%params);

    if($self->debug()){
	for my $key (sort keys %{$self->getParams()}){
	    $self->printDebug("$key => $params{$key}");
	}
    }

    $self->usebioperl(1);

    return $self;

}



=head1 parse

Description: Parse Nucmer(.coords) results.
             It returns a hash table containing each result,
  Arguments: $in results file.
    Returns: 0, $features on success
             1, message on error

=cut

sub parse($$){

    my ($self, $ifile) = @_;

    unless($self->_checkApplication($self->appl())){
	my $list = join(', ', sort keys %${GPI::ApplInfo2GFFFeat});
	return(1, "Application name not set or not supported!: \n\tSupported: $list\n");
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
	return(1, "No input file given to parse");
    }

    my $ifh;
    unless(open($ifh, "<$ifile")){
	return(1, "Could not open input file : $!");
    }

    local $OUTPUT_AUTOFLUSH = 1;
    my $hash = { };

#Output example
#    /home/projects/magnaporthe/data/db/fasta/Magnaporthe/Broad/MG_2_BCS.fasta /home/projects/magnaporthe/data/Data_Transfered/Agilent/Oligos/AgV1_Oligos.fasta
#NUCMER

#    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS]
#=====================================================================================
#    1556     1615  |        1       60  |       60       60  |   100.00  | MG_CONTIG_2.100      AZM10392
#    4529     4588  |        1       60  |       60       60  |   100.00  | MG_CONTIG_2.100      AZM09304
#    6356     6415  |       60        1  |       60       60  |   100.00  | MG_CONTIG_2.100      AZM08214
#    1427     1486  |        1       60  |       60       60  |   100.00  | MG_CONTIG_2.1000     AZM11478
#    1578     1637  |        1       60  |       60       60  |   100.00  | MG_CONTIG_2.1001     AZM12573
#    8914     8973  |       60        1  |       60       60  |   100.00  | MG_CONTIG_2.1001     AZM01104


    my $bool = 0;
    my $QUERY = 0;
    my $SET = my $PART = 0;

    while(my $line = <$ifh>){

	chomp;

	if($bool){
	    my($s1, $e1, $s2, $e2, $len1, $len2, $pid, $query, $subj) = ($line =~ /^\s+(\S+)\s+(\S+)\s+\|\s+(\S+)\s+(\S+)\s+\|\s+(\S+)\s+(\S+)\s+\|\s+(\S+)\s+\|\s+(\S+)\s+(\S+)/);

	    if(exists($hash->{$query})){

		####################
		## Match Set Tags ##
		####################
		my $stag = { };
		$stag->{ID} = join("_",
				   $query,
				   $subj,
				   $self->getApplInfo('set') .
				   sprintf('%0'.$GPI::Parser::NAMELEN.'d', ++$SET)
				  );
		$stag->{Name} = $subj;
		$stag->{target_pident} = $pid;
		$stag->{target_length} = $len2;
		$stag->{target_start}  = $s2 > $e2 ? $e2 : $s2;
		$stag->{target_end}    = $s2 > $e2 ? $s2 : $e2;
		$stag->{target_id}     = $subj;
		$stag->{program}       = $self->appl();
		$stag->{Target}        = join("+", $subj, $s2 > $e2 ? $e2 : $s2, $s2 > $e2 ? $s2 : $e2);

		#####################
		## Match Part Tags ##
		#####################
		my $ptag = { };
		$ptag->{ID} = join("_",
				   $query,
				   $subj,
				   $self->getApplInfo('set') .
				   sprintf('%0'.$GPI::Parser::NAMELEN.'d', $SET),
				   $self->getApplInfo('part') .
				   sprintf('%0'.$GPI::Parser::NAMELEN.'d', ++$PART)
				  );
		$ptag->{Name} = $subj;
		$ptag->{Target}        = join("+", $subj, $s2 > $e2 ? $e2 : $s2, $s2 > $e2 ? $s2 : $e2);
		$ptag->{Parent}        = $stag->{ID};

		my $s = Bio::SeqFeature::Generic->new(
						      -seq_id      => $query,
						      -source      => $self->appl() || 'Numcer',
						      -primary_tag => $self->qtype() || $self->getApplInfo('set'),
						      -start       => $s1,
						      -end         => $e1,
						      -score       => '.',
						      -strand      => $s2 > $e2 ? '-' : '+',
						      -frame       => '.',
						      -tag         => $stag,
						     );

		my $p = Bio::SeqFeature::Generic->new(
						      -seq_id      => $query,
						      -source      => $self->appl() || 'Numcer',
						      -primary_tag => $self->qtype() || $self->getApplInfo('part'),
						      -start       => $s1,
						      -end         => $e1,
						      -score       => '.',
						      -strand      => $s2 > $e2 ? '-' : '+',
						      -frame       => '.',
						      -tag         => $ptag,
						     );

		if($self->flush()){
		    print $writer->write_feature($s);
		    print $writer->write_feature($p);
		}else{
		    push @$features, $s, $p;
		}
	
	    }else{
		$hash->{$query} = 1;

		############
		## Region ##
		############
		#my $ref = Bio::SeqFeature::Generic->new(
		#					-seq_id      => $query,
		#					-source      => $self->appl() || 'Numcer',
		#					-primary_tag => $self->qtype() || $self->getApplInfo('region'),
		#					#-start      => 1,
		#					#-end        => $result->query_length(),
		#					-score      => '.',
		#					-strand     => '+',
		#					-frame      => '.',
		#					-tag        => {
		#							ID => $query . sprintf('%0'.$GPI::Parser::NAMELEN.'d', ++$QUERY),
		#						       },
		#				       );

		####################
		## Match Set Tags ##
		####################
		my $stag = { };
		$stag->{ID} = join("_",
				   $query,
				   $subj,
				   $self->getApplInfo('set') .
				   sprintf('%0'.$GPI::Parser::NAMELEN.'d', ++$SET)
				  );
		$stag->{Name} = $subj;
		$stag->{target_pident} = $pid;
		$stag->{target_length} = $len2;
		$stag->{target_start}  = $s2 > $e2 ? $e2 : $s2;
		$stag->{target_end}    = $s2 > $e2 ? $s2 : $e2;
		$stag->{target_id}     = $subj;
		$stag->{program}       = $self->appl();
		$stag->{Target}        = join("+", $subj, $s2 > $e2 ? $e2 : $s2, $s2 > $e2 ? $s2 : $e2);

		#####################
		## Match Part Tags ##
		#####################
		my $ptag = { };
		$ptag->{ID} = join("_",
				   $query,
				   $subj,
				   $self->getApplInfo('set') .
				   sprintf('%0'.$GPI::Parser::NAMELEN.'d', $SET),
				   $self->getApplInfo('part') .
				   sprintf('%0'.$GPI::Parser::NAMELEN.'d', ++$PART)
				  );
		$ptag->{Name} = $subj;
		$ptag->{Target}        = join("+", $subj, $s2 > $e2 ? $e2 : $s2, $s2 > $e2 ? $s2 : $e2);
		$ptag->{Parent}        = $stag->{ID};

		my $s = Bio::SeqFeature::Generic->new(
						      -seq_id      => $query,
						      -source      => $self->appl() || 'Numcer',
						      -primary_tag => $self->qtype() || $self->getApplInfo('set'),
						      -start       => $s1,
						      -end         => $e1,
						      -score       => '.',
						      -strand      => $s2 > $e2 ? '-' : '+',
						      -frame       => '.',
						      -tag         => $stag,
						     );
		my $p = Bio::SeqFeature::Generic->new(
						      -seq_id      => $query,
						      -source      => $self->appl() || 'Numcer',
						      -primary_tag => $self->qtype() || $self->getApplInfo('part'),
						      -start       => $s1,
						      -end         => $e1,
						      -score       => '.',
						      -strand      => $s2 > $e2 ? '-' : '+',
						      -frame       => '.',
						      -tag         => $ptag,
						     );

		if($self->flush()){
		    #print $writer->write_feature($ref);
		    print $writer->write_feature($s);
		    print $writer->write_feature($p);
		}else{
		    #push @$features, $ref, $s, $p;
		    push @$features, $s, $p;
		}
	    }

	    next;
	}

	if($line =~ /^\s*================/){
	    $bool = 1;
	}

    }

    return(0, $features);

} #end parse
