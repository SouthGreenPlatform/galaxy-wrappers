=head1 NAME

InterProScan - Module to parse InterProScan results, from the merged.raw file.

=head1 SYNOPSIS

#example here
use GPI::Parser::InterProScan;

=head1 DESCRIPTION

This module can be used only by using the merged results file called merged.raw from
InterProScan.

=head1 VERSIONS

$Id: InterProScan.pm,v 1.22 2007/03/01 09:28:03 equevill Exp $
$Log: InterProScan.pm,v $
Revision 1.22  2007/03/01 09:28:03  equevill
Added CVS log in the header of the file


Copyright (c) INRA/URGI 2005

=head1 AUTHORS / ACKNOWLEDGEMENTS

Emmanuel Quevillon <emmanuel.quevillon@versailles.inra.fr>

=cut


package GPI::Parser::InterProScan;
use lib '/apps/GnpAnnot/lib';
use strict;
use Data::Dumper;
use vars qw(@ISA);
use GPI::Parser;
@ISA = qw( GPI::Parser );


=head1 new

Descrption: Create a new GPI::Parser::InterProScan object.
 Arguments: $params reference to hash table with parameters.
   Returns: Reference to GPI::Parser::InterProScan object or undef.

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

Description: Parses InterProScan results.
  Arguments: $in InterProScan results file.
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
    my $bioperl = $self->usebioperl();

    require Bio::SeqFeature::Generic;
    use GPI::GFF;

    if($self->flush()){
	$writer = GPI::GFF->new(
					    -gff_version => 3,
					    -fh          => $self->fh() || \*STDOUT
					   );
    }

    my ($results, $seq_info) = $self->_parse_raw($ifile);

    my $databases = {
		     'HMMPfam'     => {
				       'db'     => 'PFAM',
				       'dbtype' => 'hmmmodel',
				       'prog'   => 'hmmpfam',
				      },
		     'HMMTigr'     => {
				       'db'     => 'TIGRFAMs',
				       'dbtype' => 'hmmmodel',
				       'prog'   => 'hmmpfam',
				      },
		     'HMMSmart'    => {
				       'db'     => 'SMART',
				       'dbtype' => 'hmmmodel',
				       'prog'   => 'hmmpfam',
				      },
		     'HMMPIR'      => {
				       'db'     => 'PIR',
				       'dbtype' => 'hmmmodel',
				       'prog'   => 'hmmpfam',
				      },
		     'HMMPanther'  => {
				       'db'     => 'PANTHER',
				       'dbtype' => 'hmmmodel',
				       'prog'   => 'hmmpfam',
				      },
		     'superfamily' => {
				       'db'     => 'SUPERFAMILY',
				       'dbtype' => 'hmmmodel',
				       'prog'   => 'hmmpfam',
				      },
		     'Gene3D'      => {
				       'db'     => 'GENE3D',
				       'dbtype' => 'hmmmodel',
				       'prog'   => 'hmmpfam',
				      },
		     'FPrintScan'  => {
				       'db'     => 'PRINTS',
				       'dbtype' => 'matrix',
				       'prog'   => 'FingerPRINTScan',
				      },
		     'SignalPHMM'  => {
				       'db'     => 'SIGNALP',
				       'dbtype' => 'model',
				       'prog'   => 'signalp',
				      },
		     'SignalP'     => {
				       'db'     => 'SIGNALP',
				       'dbtype' => 'model',
				       'prog'   => 'signalp',
				      },
		     'TMHMM'       => {
				       'db'     => 'TMHMM',
				       'dbtype' => 'model',
				       'prog'   => 'tmhmm',
				      },
		     'Coil'        => {
				       'db'     => 'COIL',
				       'dbtype' => 'matrix',
				       'prog'   => 'ncoils',
				      },
		     'Seg'         => {
				       'db'     => 'SEG',
				       'dbtype' => 'na',
				       'prog'   => 'seg',
				      },
		     'ScanRegExp'  => {
				       'db'     => 'PROSITE',
				       'dbtype' => 'strings',
				       'prog'   => 'scanregexp',
				      },
		     'ProfileScan' => {
				       'db'     => 'PROFILE',
				       'dbtype' => 'strings',
				       'prog'   => 'pfscan',
				      },
		     'PfScan'      => {
				       'db'     => 'PROFILE',
				       'dbtype' => 'strings',
				       'prog'   => 'pfscan',
				      },
		     'BlastProDom' => {
				       'db'     => 'PRODOM',
				       'dbtype' => 'sequences',
				       'prog'   => 'blastp',
				      },
		    };


    ####################################
    ## GFF3 formatting ... #############
    ####################################

    my $delim         = "\t";
    my $LineSeparator = "###\n";
    my $source        = "InterProScan";
    my $RefType       = "polypeptide";
    my $SetType       = "polypeptide";
    my $PartType      = "polypeptide_domain";
    my $RefTypeCnt    = 0;
    my $PartTypeCnt   = 1;
    my $numLen        = 4;

    #print "##gff-version 3\n";

    $self->_print_verbose("\tGoing through results ...");

    foreach my $seq (sort keys %{$results}) {
	
	next unless($seq);

	$RefTypeCnt++;
	my $RefNumFeat     = sprintf('%06d', $RefTypeCnt);
	my $RefNumTypeFeat = join("_", $seq, $RefType, $RefNumFeat);
	
	my $RefTag = { };
	$RefTag->{ID} = $seq || $self->getApplInfo('region') . sprintf('%0'.$NAMELEN.'d', $RefTypeCnt);

	my $RefFeat = Bio::SeqFeature::Generic->new(
						    -seq_id      => $seq,
						    -primary_tag => $self->qtype() || $self->getApplInfo('region'),
						    -source_tag  => $self->appl() || "InterProScan",
						    -start       => 1,
						    -end         => $seq_info->{$seq}{len},
						    -frame       => '.',
						    -score       => '.',
						    -strand      => '+',
						    -tag         => $RefTag
						   );


	unless($self->flush()){
	    push(@$features, $RefFeat);
	}else{
	    print $writer->write_feature($RefFeat);
	}

	$PartTypeCnt = 1;

	foreach my $k (sort keys %{$results->{$seq}}) {

	    next unless($k);

	    if ($k =~ /^(IPR|NULL)/) {
		
		my $FeatureNote;
		
		if ($results->{$seq}{$k}{'name'} && $results->{$seq}{$k}{'name'} ne 'NULL') {
		    my $name = $results->{$seq}{$k}{'name'};
		    $name =~ s/[&\|\%\#\~\,\;]//g;
		    $FeatureNote = join(" ", "InterPro", $k, $name);
		}

		my $FeatureGoClass;

		if ($results->{$seq}{$k}{'class'}) {
		    my $GoClass = [ ];
		    @$GoClass =  grep { s/\((GO:\d+)\)/$1/m } split(" ", $results->{$seq}{$k}{'class'});
		    $FeatureGoClass = join("", @$GoClass);
		}
		

		foreach my $mm (sort keys %{$results->{$seq}{$k}{'meth'}}) {

		    ## We don't use this part anymore as we only have 2 type of features:
		    # - polypeptide
		    # - polypeptide_domain

		    foreach my $hit (sort keys %{$results->{$seq}{$k}{'meth'}{$mm}}) {
			
			foreach my $loc (@{$results->{$seq}{$k}{'meth'}{$mm}{$hit}{'loc'}}){
			    
			    next unless($loc);
			    
			    my ($start, $end, $score) = $loc =~ /\[(\d+)-(\d+)\]\s+(\S+)?/;
			    $score = $score eq 'NA' ? '.' : $score;

			    my $PartFeat;
			    my $PartTag = { };

			    #Autogenerated number...
			    my $PartNumFeat     = sprintf('%0'.$numLen.'d', $PartTypeCnt);
			    my $PartNumTypeFeat = $PartType . $PartNumFeat;
			    my $PartFeatID      = join("_", $RefNumTypeFeat, $PartNumTypeFeat);
			    my $PartFeatName    = join("_", $mm, $results->{$seq}{$k}{'meth'}{$mm}{$hit}{'name'}, $PartNumTypeFeat);
			    $PartFeatName       =~ s/[&\|\%\#\~\,\;]//g;

			    #my $PartTargetName = join("_", $SetTargetName, $PartNumFeat, $start, $end);
			    my $PartTargetName = join("_", $hit, $RefNumFeat, $PartNumFeat, $start, $end);
			    $PartTag->{ID}            = $PartFeatID;
			    #$PartTag->{Target}        = join("+", $PartTargetName, $start, $end);
			    #$PartTag->{Name}          = join("_", $PartFeatName, $RefNumFeat);
			    $PartTag->{target_id}     = $hit;
			    $PartTag->{target_desc}   = $results->{$seq}{$k}{'meth'}{$mm}{$hit}{'name'};
			    $PartTag->{Name}          = join("_", $seq, $hit);
			    $PartTag->{Note}          = $FeatureNote if($FeatureNote);
			    $PartTag->{Ontology_term} = $FeatureGoClass if($FeatureGoClass);
			    $PartTag->{program}       = $databases->{$mm}->{prog};
			    $PartTag->{lib}           = $databases->{$mm}->{dbtype};
			    $PartTag->{Dbxref}        = $databases->{$mm}->{db} . ":$hit";
			    $PartTag->{Dbxref}       .= "," . $self->gbrowseconf() . ":$seq" if($self->gbrowseconf());


			    if(exists($databases->{$mm})){

				my $ptag = $databases->{$mm}->{prog} =~ /signalp/ ? 'signal_peptide' : $self->getApplInfo('part');

				my $PartBSG = Bio::SeqFeature::Generic->new(
									    -seq_id      => $seq,
									    -source_tag  => $self->appl()  || "InterProScan",
									    -primary_tag => $self->qtype() || $ptag,
									    -start       => $start,
									    -end         => $end,
									    -score       => $score,
									    -frame       => '.',
									    -strand      => '+',
									    -tag         => $PartTag
									   );

								
				unless($self->flush()){
				    push(@$features, $PartBSG);
				}else{
				    $writer->write_feature($PartBSG);
				}
				
			    }
			    else{
				$self->_exitOnError("Database [$mm] not supported or does not fit supported databases.");
			    }
			    $PartTypeCnt++;
			} #end foreach
		    }
		}
	    }
	    else {

		foreach my $hit (sort keys %{$results->{$seq}{$k}}){

		    foreach my $loc (@{$results->{$seq}{$k}{$hit}{'loc'}}) {
			
			my ($start, $end, $score) = $loc =~ /\[(\d+)-(\d+)\]\s+(\S+)?/;
			$score = $score eq 'NA' ? '.' : $score;

			#Autogenerated number...
			my $PartFeat;
			my $PartTag = { };
			my $PartNumFeat     = sprintf('%0'.$numLen.'d', $PartTypeCnt);
			my $PartNumTypeFeat = $PartType . $PartNumFeat;
			my $PartFeatID      = join("_", $RefNumTypeFeat, $PartNumTypeFeat);
			my $PartFeatName    = join("_", $k, $results->{$seq}{$k}{$hit}{'name'}, $PartNumTypeFeat);
			$PartFeatName       =~ s/[&\|\%\#\~\,\;]//g;
			my $PartTargetName = join("_", $hit, $RefNumFeat, $PartNumFeat, $start, $end);
			
			$PartFeatName = $PartFeatName =~ /no description/i ? join("_", $k, $hit, $PartNumTypeFeat) :  $PartFeatName;
			$PartFeat     = "ID=$PartFeatID;Target=$PartTargetName+$start+$end;Name=" . join("_", $PartFeatName, $RefNumFeat);

			$PartTag->{ID}            = $PartFeatID;
			#$PartTag->{Target}        = join("+", $PartTargetName, $start, $end);
			#$PartTag->{Name}          = join("_", $PartFeatName, $RefNumFeat);
			$PartTag->{target_id}     = $hit;
			$PartTag->{target_desc}   = $results->{$seq}{$k}{$hit}{'name'};
			$PartTag->{Name}          = join("_", $seq, $hit);
			#$PartTag->{Note}          = $SetFeatNote if($SetFeatNote);
			#PartTag->{Ontology_term} = $SetFeatGoClass if($SetFeatGoClass);
			$PartTag->{program}       = $databases->{$k}->{prog};
			$PartTag->{lib}           = $databases->{$k}->{dbtype};
			$PartTag->{Dbxref}        = $databases->{$k}->{db} . ":$hit";
			$PartTag->{Dbxref}       .= "," . $self->gbrowseconf() . ":$seq" if($self->gbrowseconf());

			my $PartBSG = Bio::SeqFeature::Generic->new(
								    -seq_id      => $seq,
								    -source_tag  => $self->appl() || "InterProScan",
								    -primary_tag => $self->qtype() || $self->getApplInfo('part'),
								    -start       => $start,
								    -end         => $end,
								    -score       => $score,
								    -frame       => '.',
								    -strand      => '+',
								    -tag         => $PartTag
								   );

								
			unless($self->flush()){
			    push(@$features, $PartBSG);
			}else{
			    print $writer->write_feature($PartBSG);
			}
			$PartTypeCnt++;
		    }
		    #print $LineSeparator;
		} #end foreach my $hit
	    } # end else
	} #end foreach my $k
    } #end foreach my $seq

    $self->_print_verbose("\tDONE");

    return($features);

}

=head1 _parse_raw

 Description: Parses InterProScan output results file (raw format) tab delimited
 Arguments:   The raw output results file
 Returns:     Hash table with parsed results.

=cut


sub _parse_raw ($){

    my($self, $ifile) = @_;


    my %results;
    my %seq_info;
    my %appl;
    my %ipr_classification;
    my $go;

    unless(open(IN,"< $ifile")){
	$self->_exitOnError("Could not open $ifile : $!");
    }


    $self->_print_verbose("\tReading raw file $ifile ... ");

    while (<IN>)
    {
	my ($seq_ac, $seq_crc, $seq_len, $_raw, $meth, $beg, $ipr) = /^(\S+)\s+(\S+)\s+(\d+)\s+((\S+)\s+\S+\t[\w\/\[\]\(\)\.\:\"\'\+\=\,\;\-\& ]*\t(\d+)\s+\d+\s+\S*\s*\S+\s+\d+\-\S+\-\d+\s*(\S+)*\s*([^\n]+)?)/;

	my $key;
	$appl{$seq_ac}{$meth} = 1 if(!exists($appl{$seq_ac}{$meth}));########### 07-10-03
	if (defined $ipr) {
	    $key = $ipr;
	} else {
	    $key = $meth;
	}

	if (defined $results{$seq_ac}{$key}) {
	    $results{$seq_ac}{$key} .= "\n". $_raw;
	} else {
	    $results{$seq_ac}{$key} = $_raw ;
	}
 	
	# $seq_ac is supposed to be unique
	$seq_info{$seq_ac}{'crc'} = $seq_crc;
	$seq_info{$seq_ac}{'len'} = $seq_len;
    }

    close(IN);


    #parsing input
    foreach my $ac (sort keys %results) {
	next unless $ac;
        foreach my $key (sort keys %{$results{$ac}}) {

            if ($key =~ /^(IPR|NULL)/) {	#if interpro
		$results{$ac}{$key} =~ /\S+\s+\S+\t[\w\/\[\]\(\)\,\;\.\:\"\'\+\=\-\& ]*\t\d+\s+\d+\s+\S*\s*\S+\s+\d+\-\S+\-\d+\s+\S+\s+([^\n\t]+)(\t+([^\n\t]+))?/;
		my ($ipr_mapping, $go_mapping_long, $go_mapping_short) = ($1, $2, $3);
		my (%ipr) = ('name' => "$ipr_mapping");
		
		if ( defined($go_mapping_short) ) {
		    $ipr{'class'} = "$go_mapping_short";
		    $go = 1;
		}
		
		my (%meth) = ();
		while ( $results{$ac}{$key} =~ /(\S+)\s+(\S+\t[\w\/\[\]\(\)\.\:\"\'\+\=\,\;\-\& ]*\t\d+\s+\d+\s+\S*\s*\S+\s+\d+\-\S+\-\d+)/g ) {
		    my ($method, $match) = ($1, $2);
		    
		    if ( defined($meth{$method}) ) {
			$meth{$method} .= "\n". "$match";
		    } else {
			$meth{$method} = "$match";
		    }
		}
		
		foreach my $m (sort keys %meth) {
		    my %hits;
		    while ($meth{$m} =~ /(\S+)\t([\w\/\[\]\(\)\.\:\"\'\+\=\,\;\-\& ]*)\t(\d+)\s+(\d+)\s+(\S*)\s*(\S+)\s+(\d+\-\S+\-\d+)/g) {
			my ($methAC, $methName, $start, $end, $methScore, $status, $date) = ($1, $2, $3, $4, $5, $6, $7);
			my $hit_desc = "${status}\[$start-$end]";
			$hit_desc .= " $methScore" if ( defined($methScore) && ($methScore !~ /^\s*$/) );
			
			if (defined($hits{$methAC})) {	# hit_ac is uniq  
			    push(@{$hits{$methAC}{'loc'}}, "$hit_desc");
			} else {
			    $hits{$methAC} = {'name' => "$methName", 'loc' => [$hit_desc]};
			}
		    }
		    $ipr{'meth'}{$m} = \%hits;
		}
		$results{$ac}{$key} = \%ipr;
		
	    } else {
		my %hits;
		
		while ($results{$ac}{$key} =~ /(\S+)\s+(\S+)\t([\w\/\[\]\(\)\.\:\"\'\+\=\,\;\-\& ]*)\t(\d+)\s+(\d+)\s+(\S*)\s*(\S+)\s+(\d+\-\S+\-\d+)/g ) {
		    if (defined $hits{$2}){	#hit_ac is uniq 
			push(@{$hits{$2}{'loc'}}, "${7}\[$4-$5\] $6");
		    } else {
			$hits{$2} = {'name' => $3, 'loc' => ["${7}[$4-$5] $6"]};	  
		    }
		}
		$results{$ac}{$key} = \%hits;
	    }
	}
    }

    $self->_print_verbose("\tDONE");

    return (\%results, \%seq_info);
}



1;
