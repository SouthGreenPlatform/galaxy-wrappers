=head1 NAME

Parser.pm - Interface modules to parse Blast(n/x), Sim4 results and create GFF3 output.

=head1 SYNOPSIS

#example here
use GPI::Parser;

=head1 DESCRIPTION

It uses BioPerl modules to parse different output results.

=head1 VERSIONS

$Id: Parser.pm,v 1.73 2007/03/16 09:10:28 equevill Exp $
$Log: Parser.pm,v $
Revision 1.73  2007/03/16 09:10:28  equevill
*** empty log message ***


Copyright (c) INRA/URGI 2005

=head1 AUTHORS / ACKNOWLEDGEMENTS

Emmanuel Quevillon <emmanuel.quevillon@versailles.inra.fr>

=cut


package GPI::Parser;
use lib '/apps/GnpAnnot/lib';
use 5.008000;
use strict;
use warnings;
use Math::BigFloat;
use Term::ProgressBar;
use Exporter;

#Version
our $VERSION = "1.45";

#Length for automtic generated names.
our $NAMELEN         = 4;
our $DELIM           = "\t";
our $MAX_SUBJ_LENGTH = 1_000_000;
our $SET_SEPARATOR   = "###\n";

use English;
use vars qw($AUTOLOAD @ISA);
use Data::Dumper;

our @ISA = qw(Exporter);

#We export by default those variables to be accessible for children modules.
our @EXPORT = qw($NAMELEN $DELIM $MAX_SUBJ_LENGTH);

our $ApplInfo2GFFFeat = {

			 'sim4'         => {
					    'id'     => 'match_part',
					    'source' => 'match_set',
					    'set'    => 'match_set',
					    'part'   => 'match_part',
					    'region' => 'region'
					   },
			 'blastn'       => {
					    'id'     => 'match_part',
					    'source' => 'match_set',
					    'set'    => 'match_set',
					    'part'   => 'match_part',
					    'region' => 'region'
					   },
			 'blastp'       => {
                                            'id'     => 'match_part',
                                            'source' => 'match_set',
                                            'set'    => 'match_set',
                                            'part'   => 'match_part',
                                            'region' => 'region'
                                           },
			 'blastx'       => {
					    'id'     => 'match_part',
					    'source' => 'match_set',
					    'set'    => 'match_set',
					    'part'   => 'match_part',
					    'region' => 'region'
					   },
			 'tblastx'       => {
					    'id'     => 'match_part',
					    'source' => 'match_set',
					    'set'    => 'match_set',
					    'part'   => 'match_part',
					    'region' => 'region'
					   },
			 'garnier'      => {
					    'id'     => 'match_part',
					    'source' => 'match_set',
					    'set'    => 'match_set',
					    'part'   => 'match_part',
					    'region' => 'region'
					   },
			 'RepeatMasker' => {
					    'id'     => 'match_part',
					    'source' => 'match_set',
					    'set'    => 'match',
					    'part'   => 'match',
					    'region' => 'region'
					   },
			 'trf'          => {
					    'id'     => 'match_part',
					    'source' => 'match_set',
					    'set'    => 'match',
					    'part'   => 'match',
					    'region' => 'region'
					   },
			 'iprscan'      => {
					    'id'     => 'match_part',
					    'source' => 'match_set',
					    'set'    => 'polypeptide',
					    'part'   => 'polypeptide_domain',
					    'region' => 'polypeptide'
					   },
			 'InterProScan' => {
					    'id'     => 'match_part',
					    'source' => 'match_set',
					    'set'    => 'polypeptide',
					    'part'   => 'polypeptide_domain',
					    'region' => 'polypeptide'
					   },
			 'csv'          => {
					    'id'     => 'match_part',
					    'source' => 'match_set',
					    'set'    => 'match_set',
					    'part'   => 'match_part',
					    'region' => 'region'
					   },
			 'nucmer'       => {
					    'id'     => 'match_part',
					    'source' => 'match_set',
					    'set'    => 'match_set',
					    'part'   => 'match_part',
					    'region' => 'region'
					   },
			 'hmmtop'       => {
					    'id'     => 'match_part',
					    'source' => 'match_set',
					    'set'    => 'match',
					    'part'   => 'polypptide_domain',
					    'region' => 'polypeptide'
					   },
			 'est2genome'   => {
					    'id'     => 'match_part',
					    'source' => 'match_set',
					    'set'    => 'match_set',
					    'part'   => 'match_part',
					    'region' => 'region',
					   },
		       };

my $ChildToParent = {
		     'exon'       => 'transcript',
		     'mRNA'       => 'transcript',
		     'HSP'        => 'match',
		     'match_part' => 'match_set',
		     'CDS'        => 'mRNA',
		     'intron'     => 'gene',
		     'cDNA'       => 'mRNA',
		     'cDNA_match' => 'match',
		     'EST'        => 'mRNA',
		     'EST_match'  => 'match'
		    };


=head1 new

Description: Creates a new object GPI::Parser. But you cannot
             use this function, you need to use the subclass
             new function, e.g. GPI::Parser::Blast->new.
  Arguments: none
    Returns: Reference to the object or undef if failed.

=cut

sub new {
    my($class) = shift;

    print STDERR "You cannot use this module directly to create an object, use subclass method.";
    return undef;
}


=head1 writeGFF3

Description: Writes a gff3 format file.
  Arguments: $hres Reference to a hash table.
             $ofh  File handle to write results (default STDOUT).
    Returns: 0, '' o success
             1, $msg on error

=cut

sub writeGFF3 ($$){

    my($self, $hres, $ofh) = @_;

    #Parser::Fasta already writes gff3 when parsing to avoid loading too much the memory
    return(0, '') if($self->appl() eq 'fasta');

    unless(ref($hres)){
	return(1, "writeGFF3 need a reference to an array or a hash table");
    }

    #To store all the seqids seen if user wants to add them add the enf of GFF3
    my $seqid = { };

    #We use BioPerl to write gff3 format
    if(ref($hres) eq 'ARRAY'){
	use GPI::GFF;

	my $writer = GPI::GFF->new(
					       -gff_version => 3,
					       -fh          => $ofh || \*STDOUT
					      );

	#Prevent to close $ofh when $writer destroyed
	$writer->noclose(1);

	return(1, "Could not create GPI::GFF object") unless($writer);

	for(@$hres){
	    #User wants FASTA sequences at the end of GFF3 ?
	    if($self->addseq() && !exists($seqid->{$_->seq_id()})){
		$seqid->{$_->seq_id()} = 1;
		$seqid->{nseq}++;
	    }

	    print $ofh "###\n" if($_->primary_tag() =~ /match_set|region/ && !$writer->{_first});
	    print $writer->write_feature($_);
	}

    }
    elsif(ref($hres) eq 'HASH'){
	#Second file handle to write the region definition into a separate file.
	my $ofh2 = $self->fh2();

	unless($self->isFileHandle($ofh)){
	    $ofh  = \*STDOUT;
	    $ofh2 = \*STDOUT;
	}

	#We set again the value in case..
	$self->fh($ofh);
	$self->fh2($ofh2);

	my($GFFwriter, $GFFwriter2, $RefTag, $SetTag, $PartTag);

	$seqid->{nseq} = $hres->{nresults};

	#To avoid complains from 'sort {$a <=> $b} keys %$hres'
	delete $hres->{nresults};

	my $appl  = $self->appl() || 'UNKNOWN';

	my $source = $appl;

	$source .= "_" . $self->annot() if($self->annot());

	print $ofh "##gff-version 3\n" unless($self->_topSeen());
	print $ofh2 "##gff-version 3\n" if ( ! $self->usebioperl() && $self->nogbrowse());

	foreach my $result (sort {$a <=> $b} keys %$hres){

	    my $hits = $hres->{$result}->{hits};
	    next unless(keys (%$hits));

	    my $top = $hres->{$result}->{topdesc};

	    if($self->addseq()){
		
		unless(exists($seqid->{$top->{name}})){

		    $seqid->{$top->{name}} = 1;
		    print STDERR "##### Name ", $top->{name}, " Fasta sequence to add later\n" if($self->verbose());

		}else{

		    print STDERR "\t********* WARNING ***********\n";
		    print STDERR "\t    ", $seqid->{$top->{name}}, "  already exists (duplicated?)\n";
		    print STDERR "\t*****************************\n\n";

		}
		
	    }

	    foreach my $nhit (sort {$a <=> $b} keys %$hits){
		
		my $hit = $hits->{$nhit};
		
		#We redefine the sequence region
		my($RefID, $RefFeat, $RefNum);
		my $hitid = $hit->{id} || 0;

		$RefNum  = sprintf('%0'.$NAMELEN.'d', $hitid);
		$RefID   = $top->{name} || $self->getApplInfo('region') . $RefNum;
		$RefFeat = "ID=$RefID\n";
		
		unless(defined($self->getParam('backname')) && $self->getParam('backname') eq $top->{name}){

		    $self->setParam('backname', $top->{name});
		    $self->printTop($top, $RefFeat);
		}

		#Now we are defining the match_set.
		my($SetStrand, $SetID, $SetNum, $SetQStart, $SetQEnd, $SetSStart, $SetSEnd, $SetEvalue, $SetFrame, $SetName, $SetTarget);
		my($SetFeat, $SetFeatID, $SetFeatTarget, $SetFeatName);
		
		$SetStrand = $hit->{strand};
		$SetStrand =~ s/\d+//;
		$SetStrand = $SetStrand eq '' ? '+' : '-';
		$SetEvalue = $hit->{evalue};

		#We extend the automatic ID by prefixing it with the annotation name (e.g. Broad, Agilent....)
		$SetNum = sprintf('%0'.$NAMELEN.'d', $hitid);
		
		#$SetID .= $self->annot() ."_" if($self->annot());
		#remove $RefID from the join, do we really need this reference id as it is associated with the ID of the region in the genome browser.
		$SetID .= join("_", $source, $RefID);
		$SetID .= "_" . $hit->{name} if($hit->{name});
		$SetID .= "_" . $hit->{extra} if($hit->{extra}); #Extra annotation
		$SetID .= "_" . join("_", $self->getApplInfo('set') . $SetNum);

		if(defined(@{$hit->{hsps}})){

		    #We sort hsps by subject start and end position just in case we need it for later
		    my $hspss = [ ];
		    my $hspse = [ ];
		    @$hspss = sort { $a->{starts} <=> $b->{starts} } @{$hit->{hsps}};
		    @$hspse = sort { $a->{ends} <=> $b->{ends} } @{$hit->{hsps}};

		    #We sort hsps by query start and end position just in case we need it for later
		    my $hspqs = [ ];
		    my $hspqe = [ ];
		    @$hspqs = sort { $a->{startq} <=> $b->{startq} } @{$hit->{hsps}};
		    @$hspqe = sort { $a->{endq} <=> $b->{endq} } @{$hit->{hsps}};
		
		    $SetQStart = $hspqs->[0]->{startq};
		    $SetQEnd   = $hspqe->[-1]->{endq};
		    $SetSStart = $hspss->[0]->{starts};
		    $SetSEnd   = $hspse->[-1]->{ends};
		    $SetEvalue = $hit->{evalue} ? $hit->{evalue} : '.';
		    #my $SetFrame  = $hit->{frame} ? $hit->{frame} : '.';
		    $SetFrame  = '.';
		
		}else{
		    $SetSStart = $hit->{starts};
		    $SetSEnd   = $hit->{ends};
		    $SetQStart = $hit->{startq};
		    $SetQEnd   = $hit->{endq};
		    $SetFrame  = $hit->{frame};
		}

		$SetFeatID  = "ID=$SetID";
		
		if($hit->{name}){

		    $SetTarget     = join('+', $hit->{name}, $SetSStart, $SetSEnd);
		    $SetFeatTarget = ";Target=$SetTarget";

		    $SetName       = $top->{name} . "_" . $hit->{name};
		    $SetFeatName   = ";Name=$SetName";
		}
		
		$SetFeat = $SetFeatID . $SetFeatTarget . $SetFeatName;
		
		if($self->haveFilters()){

		    my @notes = ( );

		    #Target Informations
		    if($hit->{name}){
			$SetFeat .= ";target_id=" . $hit->{name} if(defined($hit->{name}));
		    }
		    if($SetSStart){
			$SetFeat .= ";target_start=" . $SetSStart if(defined($SetSStart));
		    }
		    if($SetSEnd){
			$SetFeat .= ";target_end=" . $SetSEnd if(defined($SetSEnd));
		    }
		    if(defined($self->addnote()) || defined($self->pcoverage())){
			$SetFeat .= ";target_pcover=" . sprintf("%.2f", $hit->{pcoverage}) if(defined($hit->{pcoverage}));
		    }
		    if(defined($self->addnote()) || defined($self->pidentity())){
			$SetFeat .= ";target_pident=" . sprintf("%.2f", $hit->{pidentity}) if(defined($hit->{pidentity}));
		    }
		    if(defined($self->addnote()) || defined($self->nconserved())){
			$SetFeat .= ";target_pcons=" . sprintf("%.2f", $hit->{nconserved}) if(defined($hit->{nconserved}));
		    }
		    if(defined($self->addnote()) || defined($self->evalue())){
			$SetFeat .= ";evalue=" . $hit->{evalue} if(defined($hit->{evalue}));
		    }
		    if(defined($self->desc())){
			$SetFeat .= ";target_description=" . $hit->{desc} if(defined($hit->{desc}));
		    }
		    if($hit->{length}){
			$SetFeat .= ";target_length=" . $hit->{length} if(defined($hit->{length}));
		    }

		    #Query Informations
		    if($self->qacc() && defined($hit->{qacc})){
			$SetFeat .= ";query_acc=".$hit->{qacc};
		    }	
		    if($self->qdesc() && defined($hit->{qdesc})){
			$SetFeat .= ";query_desc=".$hit->{qdesc};
		    }	
		    if($hit->{gapst}){
			push @notes, "Gaps:$hit->{gapst}";
		    }
		    if($self->annot()){
			$SetFeat .= ";lib=" . $self->annot() if(defined($self->annot()));
		    }
		    if($appl){
			$SetFeat .= ";program=" . $appl if(defined($appl));
		    }

		    #Our list of private tags
		    if(exists($hit->{_private_tags}) && scalar(@{$hit->{_private_tags}})){
			$SetFeat .= ";" . join(";", @{$hit->{_private_tags}});
		    }
		    if(scalar(@notes)){
			$SetFeat .= ";Note=" . join(",", @notes);
		    }
		}

		if($self->gbrowseconf()){
		    $SetFeat .= ";Dbxref=" . $self->gbrowseconf() . ":" . $SetName;

		}

		print $ofh join($DELIM,
				$top->{name},
				$source,
				$self->stype() || $self->getApplInfo('set'),
				$SetQStart,
				$SetQEnd,
				$SetEvalue || '.',
				$SetStrand,
				$SetFrame || '.',
				"$SetFeat\n");
		
		next unless(defined(@{$hit->{hsps}}));

		#We sort the hsp to have them order with their respecting start regarding query start.
		foreach my $hsp (sort { $a->{startq} <=> $b->{startq} } @{$hit->{hsps}}){

		    #Here we define the match_part

		    my $PartStrand = $hsp->{strand};
		    $PartStrand =~ s/\d+//;
		    $PartStrand = $PartStrand eq '' ? '+' : '-';

		    my $PartFrame  = $hsp->{frame} || '.';
		    $PartFrame     =~ s/\D+// if($PartFrame);
		    #$PartFrame     = '.';

		    my $PartSStart = $hsp->{starts};
		    my $PartSEnd   = $hsp->{ends};
		    my $PartQStart = $hsp->{startq};
		    my $PartQEnd   = $hsp->{endq};
		    my $PartEvalue = $hsp->{evalue} ? $hsp->{evalue} : '.';

		    my $PartID;
		    my $PartNum = sprintf('%0'.$NAMELEN.'d', $hsp->{id});
		    $PartID .= join("_",
				    $SetID,
				    $self->getApplInfo('part') . $PartNum
				   );

		    #We split the 9th field can be usefull? We don't write a 'Name' part as we use 'Parent' instead.
		    my ($PartFeat, $PartFeatID, $PartFeatParent, $PartFeatTarget);
		    my($PartParent, $PartTarget);

		    $PartFeatID     = "ID=$PartID";

		    if($hsp->{name}){
			$PartTarget     = join('+', $hsp->{name}, $PartSStart, $PartSEnd);
			$PartFeatTarget = ";Target=" . join('+', $hsp->{name}, $PartSStart, $PartSEnd);
		    }

		    $PartParent     = $SetID;
		    $PartFeatParent = ";Parent=$PartParent";

		    $PartFeat = $PartFeatID . $PartFeatTarget . $PartFeatParent;
			
		    print $ofh join($DELIM,
				    $top->{name},
				    $source,
				    $self->stype() ? $self->stype() . "_part" : $self->getApplInfo('part'),
				    $PartQStart || 0,
				    $PartQEnd || 0,
				    $PartEvalue || '.',
				    $PartStrand || '.',
				    $PartFrame || '.',
				    "$PartFeat\n");
		} #end foreach hsp
		print $ofh $SET_SEPARATOR;
	    } #end foreach hit
	}

    }else{
	return(1, "You need to provide either a hash reference or an array reference");
    }

    #We want to add Fasta sequences at the end of the gff file.
    if($self->addseq()){
	my($res, $mess) = $self->write_sequences($seqid);
	return (1, "Problem while writing FASTA sequences..\n$mess") if($res);
    }

    return(0, '');
}


=head1 write_sequences

Description: Write FASTA sequences for GFF
  Arguments: Array reference of Bio::Seq objects
             or
             Hash reference with seqid, but requires the fasta file sequences
    Returns: 1, message on error
             0, '' otherwise

=cut


sub write_sequences ($) {

    require Bio::SeqIO;

    my($self, $seqs) = @_;


    #Seq writer
    my $seqo = Bio::SeqIO->new(-format => 'fasta',
			       -fh     => $self->fh()
			      );

    my $ofh = $self->fh();

    if($seqs && ref($seqs) eq 'ARRAY'){

	#FASTA delimiter for gff3
	print $ofh "##FASTA\n";

	foreach my $seq (@$seqs){
	    $seqo->write_seq($seq);
	}

	$seqo->close();

    }elsif($seqs && ref($seqs) eq 'HASH'){

	my $seqi = Bio::SeqIO->new(
				   -file   => $self->addseq(),
				   -format => 'fasta'
				  );

	unless($seqi){
	    return(1, "Could not initialise Bio::SeqIO object when writing fasta sequence : $!");
	}

	my $seqSeen = 0;
	my $res     = $seqs->{nseq};

	print STDERR "\tAdding [$res] Fasta sequences : " if($self->verbose());

	#FASTA delimiter for gff3
	print $ofh "##FASTA\n";

	#Using Files
	while(my $seq = $seqi->next_seq()){

	    if(exists($seqs->{$seq->id()})){

		$seqSeen = 1;

		#Writes the Fasta formated sequence on the file handle (default STDOUT if no file given).
		$seqo->write_seq($seq);

		$seqs->{nseq}--;

	    }
		
	    last if($seqs->{nseq} == 0);
	}

	print STDERR "Fasta sequences added.\n" if($self->verbose() && $seqSeen);
	$seqo->close();

	return (1, "No seq id seen in the output parsed file have been found in the file ". $self->addseq() . "\n") unless($seqSeen);

	
    }
    else{
	return (1, "write_sequences: Need an array or hash reference.");
    }

    return (0, '');
}

=head1 printTop

 Description: Print the top of the GFF3 (feature definition)
              NOTE: printRegion is called at the end to set the
              value to 0 to avoid duplicated lines.
   Arguments: $top hash reference
              $RefFeat An optional 9th field
     Returns: 1, error
             0, formatted string

=cut

sub printTop ($) {

    my($self, $top, $RefFeat) = @_;

    return unless($self->printRegion());

    die ("Need a  HASH reference to print top.") unless(ref($top) eq 'HASH');

    my $ofh2 = $self->fh2();
    my $source = $self->appl();
    $source .= "_" . $self->annot() if($self->annot());

    print $ofh2 join($DELIM,
		     $top->{name},
		     $source,
		     $self->qtype() || $self->getApplInfo('region'),
		     $top->{start},
		     $top->{end},
		     '.',
		     '+',
		     '.',
		     $RefFeat || '',
		    );

    #If uncommented, print region only once at the begining of the file.
    #$self->printRegion(0);
}

=head1 isFileHandle

Description: Checks whether parameter is a file handle.
  Arguments: $ofh File handle
    Returns: 1 for FileHandle reference or FileHandle glob.
             0 otherwise

=cut

sub isFileHandle ($) {

  my($self, $fh) = @_;

  my $is_handle = (ref($fh) ? (ref($fh) eq 'GLOB'
			       || UNIVERSAL::isa($fh, 'GLOB')
			       || UNIVERSAL::isa($fh, 'IO::Handle'))
		   : (ref(\$fh) eq 'GLOB'));

  return $is_handle;

}

=head1 getApplInfo

Description: Get Informations for application part (set, part, region, ...)
  Arguments: Informations wanted
    Returns: a string

=cut

sub getApplInfo ($) {

    my($self, $name) = @_;

    return('UNKNOWN') unless(defined($name));

    return $ApplInfo2GFFFeat->{$self->appl()}->{$name} || 'UNKNOWN';

}

=head1 getParams

 Description: Get all the parameters.
   Arguments: None
     Returns: A reference to the parameters hash table

=cut

sub getParams (){

  my($self) = @_;

  return unless($self->{params});

  my $backup = $self->{params};

  return $backup;

  }

=head1 getParam

Description: Get a specific parameter.
  Arguments: A key refering to a parameter
    Returns: The value of the specific key
             Nothing is key does not exists.

=cut

sub getParam ($){

  my($self, $key) = @_;

  return unless defined($key);

  my $backup = $self->getParams();

  return $backup->{$key};

  }

=head1 setParams

Description: Set a parameters hash table.
  Arguments: A reference to a hash table.
    Returns: Nothing

=cut

sub setParams ($){

  my($self, $params) = @_;

  return unless($params && ref($params) eq 'HASH');

  $self->{params} = $params;

  return;

  }


=head1 setParam

 Description: Set a specific parameter.
   Arguments: A key and the associated value.
     Returns: Nothing

=cut

sub setParam ($$){

  my($self, $key, $value) = @_;

  return unless (defined($key) && defined($value));

  my $backup = $self->getParams();

  $backup->{$key} = $value;

  $self->setParams($backup);

  return;

  }


=head2 checkParams

 Description: Check parameters passed to the constructor
   Arguments: none
     Returns: 0, '' if no errors
              1, message otherwise

=cut

sub checkParams {

    my($self) = @_;

    $self->_exitOnError("No parameters set") unless($self->getParams());

    my $params = $self->getParams();

    foreach my $param (keys %$params){
	next unless defined($param);

	if($param =~ /fh/){
	    $self->_exitOnError("[$param] is not a real file handle : " . ref($params->{$param})) unless($self->isFileHandle($params->{$param}));
	}elsif($param eq 'in'){
	    $self->_exitOnError("$params->{$param} could not be found : $!") unless(-f $params->{$param} && -e $params->{$param});
	}elsif($param =~ /pidentity|pcoverage|nconserved|maxhit|evalue/){
	    unless(!defined($params->{$param})){
		$self->_exitOnError("$param must be greater than 0 ($params->{$param})") unless($params->{$param} >= 0);
	    }
	}
    }

    return(0, '');

}

=head2 printDebug

Description: Prints debugging messages.
  Arguments: 'key' What you want to print
             $key  Reference to the object you want to print
             infos from
    Retruns: debugging messages

=cut

sub printDebug ($$) {

    my($self, $ref) = @_;

    return unless($ref);

    local $OUTPUT_FIELD_SEPARATOR = ' ';

    if(ref($ref) =~ /Result/){
	my $tab = "";
	$self->_print_debug( "${tab}======================================");
	$self->_print_debug( "${tab}Query Name       : ". $ref->query_name());
	$self->_print_debug( "${tab}Query Accession  : ". $ref->query_accession());
	$self->_print_debug( "${tab}Query Length     : ". $ref->query_length());
	$self->_print_debug( "${tab}Query Description: ". $ref->query_description());
	$self->_print_debug( "${tab}Database Name    : ". $ref->database_name());
	$self->_print_debug( "${tab}Database Entries : ". $ref->database_entries());
	$self->_print_debug( "${tab}Database Letters : ". $ref->database_letters());
	$self->_print_debug( "========================================\n");
    }
    elsif(ref($ref) =~ /Hit/){
	my $tab = "\t";
	$self->_print_debug( "${tab}**************************************");
	$self->_print_debug( "${tab}Hit Name         : ". $ref->name());
	$self->_print_debug( "${tab}Hit Description  : ". $ref->description());
	$self->_print_debug( "${tab}Hit Length total : ". $ref->length());
	$self->_print_debug( "${tab}Hit Algorithm    : ". $ref->algorithm());
	$self->_print_debug( "${tab}Hit Score        : ". $ref->raw_score());
	$self->_print_debug( "${tab}Hit Score(bits)  : ". $ref->bits() || "");
	$self->_print_debug( "${tab}Hit Signifiance  : ". $ref->significance() || "" );
	$self->_print_debug( "${tab}Hit N of Refs    : ". $ref->n());
	$self->_print_debug( "${tab}Hit Frame        : ". $ref->frame());
	$self->_print_debug( "${tab}Hit Matches(id)  : ". $ref->matches('id'));
	$self->_print_debug( "${tab}Hit Matches(cons): ". $ref->matches('cons'));
	#$self->_print_debug( "${tab}Hit Start (subj) : ". $ref->start('subject'));
	#$self->_print_debug( "${tab}Hit End (subj)   : ". $ref->end('subject'));
	#$self->_print_debug( "${tab}Hit Start (query): ". $ref->start('query'));
	#$self->_print_debug( "${tab}Hit End (query)  : ". $ref->end('query'));
	$self->_print_debug( "${tab}Hit Seq IndiIden : ". join(":", $ref->seq_inds('hit', 'identical',1)));
	$self->_print_debug( "${tab}Hit Seq IndiCons : ". join(":", $ref->seq_inds('hit', 'conserved',1)));	
	$self->_print_debug( "${tab}Hit Strand query : ". $ref->strand('query'));
	$self->_print_debug( "${tab}Hit Strand subj  : ". $ref->strand('subject'));
	$self->_print_debug( "${tab}Hit Rank         : ". $ref->rank());
	$self->_print_debug( "${tab}**************************************\n");
    }elsif(ref($ref) =~ /HSP/){
	my $tab = "\t\t";
	$self->_print_debug( "${tab}--------------------------------------");
	$self->_print_debug( "${tab}Hsp Algorithm              : ". $ref->algorithm());
	$self->_print_debug( "${tab}Hsp Evalue                 : ". $ref->evalue() || '');
	$self->_print_debug( "${tab}Hsp Frac_identical total   : ". $ref->frac_identical('total'));
	$self->_print_debug( "${tab}Hsp Frac_identical query   : ". $ref->frac_identical('query'));
	$self->_print_debug( "${tab}Hsp Frac_identical subject : ". $ref->frac_identical('hit'));
	$self->_print_debug( "${tab}Hsp Frac_conserved total   : ". $ref->frac_conserved('total'));
	$self->_print_debug( "${tab}Hsp Frac_conserved query   : ". $ref->frac_conserved('query'));
	$self->_print_debug( "${tab}Hsp Frac_conserved subject : ". $ref->frac_conserved('hit'));
	$self->_print_debug( "${tab}Hsp Num_identical          : ". $ref->num_identical());
	$self->_print_debug( "${tab}Hsp Num_conserved          : ". $ref->num_conserved());
	$self->_print_debug( "${tab}Hsp Gaps total             : ". $ref->gaps('total'));
	$self->_print_debug( "${tab}Hsp Gaps query             : ". $ref->gaps('query'));
	$self->_print_debug( "${tab}Hsp Gaps subject           : ". $ref->gaps('hit'));
	$self->_print_debug( "${tab}Hsp Query String           : ". $ref->query_string());
	$self->_print_debug( "${tab}Hsp Homology String        : ". $ref->homology_string());
	$self->_print_debug( "${tab}Hsp Hit String             : ". $ref->hit_string());
	$self->_print_debug( "${tab}Hsp Length total           : ". $ref->length('total'));
	$self->_print_debug( "${tab}Hsp Length query           : ". $ref->length('query'));
	$self->_print_debug( "${tab}Hsp Length subject         : ". $ref->length('hit'));
	$self->_print_debug( "${tab}Hsp Percent Identity       : ". $ref->percent_identity());
	$self->_print_debug( "${tab}Hsp Get Aln                : ". $ref->get_aln());
	$self->_print_debug( "${tab}Hsp Seq Ind 'query' ident  : ". join(":", $ref->seq_inds('query', 'identical', 1)));
	$self->_print_debug( "${tab}Hsp Seq Ind 'query' cons   : ". join(":", $ref->seq_inds('query', 'conserved', 1)));
	$self->_print_debug( "${tab}Hsp Seq Ind 'subj'  ident  : ". join(":", $ref->seq_inds('subject', 'identical', 1)));
	$self->_print_debug( "${tab}Hsp Seq Ind 'subj'  cons   : ". join(":", $ref->seq_inds('subject', 'conserved', 1)));
	$self->_print_debug( "${tab}Hsp Strand query           : ". $ref->strand('query'));
	$self->_print_debug( "${tab}Hsp Strand subject         : ". $ref->strand('subject'));
	$self->_print_debug( "${tab}Hsp Strand list            : ". $ref->strand('list'));
	$self->_print_debug( "${tab}Hsp Start query            : ". $ref->start('query'));
	$self->_print_debug( "${tab}Hsp End   query            : ". $ref->end('query'));
	$self->_print_debug( "${tab}Hsp Start subject          : ". $ref->start('subject'));
	$self->_print_debug( "${tab}Hsp End   subject          : ". $ref->end('subject'));
	$self->_print_debug( "${tab}Hsp Seq Str                : ". $ref->seq_str());
	$self->_print_debug( "${tab}Hsp Rank                   : ". $ref->rank());
	$self->_print_debug( "${tab}Hsp Matches                : ". $ref->matches());
	$self->_print_debug( "${tab}Hsp Matches query          : ". $ref->matches('query'));
	$self->_print_debug( "${tab}Hsp Matches subject        : ". $ref->matches('subject'));
	$self->_print_debug( "${tab}Hsp N                      : ". $ref->n());
	$self->_print_debug( "${tab}Hsp Range query            : ". $ref->range('query'));
	$self->_print_debug( "${tab}Hsp Range subject          : ". $ref->range('subject'));
	$self->_print_debug( "${tab}--------------------------------------\n");
    }else{
	return;
    }

}


=head1 FILTERING METHODS

=cut

=head2 haveFilters

Description: Check if any filter or addnote is required by the user
Argument:
Returns: 1 if filters required
         0 otherwise

=cut

sub haveFilters {

    my($self) = @_;

    if(defined($self->addnote()) || defined($self->evalue()) || defined($self->pcoverage()) || defined($self->pidentity())
       || defined($self->nconserved()) || defined($self->desc())){
	return 1;
    }

    return 0;

}


=head2 applyFilters

Description: Checks if any filters need to be applied.
  Arguments: 
   Returns : 1 if not filters need to be applied or if the results pass the filter(s)
             0 otherwise

=cut

sub applyFilters ($) {

    my ($self) = shift;
    my ($feat) = shift;

    unless(ref($feat) eq 'HASH'){
	return $self->_BP_applyFilters($feat);
    }


    if(defined($self->pcoverage()) || defined($self->pidentity()) || defined($self->nconserved()) || defined($self->evalue())){

	#We consider that all the filters are ok by default.
	my $ok = 1;

	if(defined($self->pidentity())){
	    if(defined($feat->{pidentity})){
		unless($feat->{pidentity}  >= $self->pidentity()){
		    $ok = 0;
		}
	    }else{
		$ok = 0;
	    }
	}
	if(defined($self->nconserved())){
	    if(defined($feat->{nconserved})){
		unless($feat->{nconserved} >= $self->nconserved()){
		    $ok = 0;
		}
	    }else{
		$ok = 0;
	    }
	}
	if(defined($self->pcoverage())){
	    if(defined($feat->{pcoverage})){
		unless($feat->{pcoverage}  >= $self->pcoverage()){
		    $ok = 0;
		}
	    }else{
		$ok = 0;
	    }
	}
	if(defined($self->evalue())){
	    if(defined($feat->{evalue})){
		
		my $fevalue = $feat->{evalue};
		my $uevalue = $self->evalue();

		#We transform the feature evalue to a float. Evalue Param is already converted
		my $bigfloat = Math::BigFloat->new($fevalue);
		$fevalue = $bigfloat->bstr();

		unless($fevalue <= $uevalue){
		    $ok = 0;
		}
	    }else{
		$ok = 0;
	    }
	}
	
	return $ok;
	
    }else{
	return 1;
    }

}

=head2 _BP_applyFilters

Description: Checks if any filters need to be applied. Using Bio::SeqFeature::Generic objects
  Arguments: 
   Returns : 1 if no filters need to be applied or if the results pass the filter(s)
             0 otherwise

=cut

sub _BP_applyFilters ($){

    my($self, $feat) = @_;

    unless($feat->isa("Bio::SeqFeature::Generic")){
	warn("Feature passed is not a Bio::SeqFeature::Generic : " . ref($feat));
	return 0;
    }


    if(defined($self->pcoverage()) || defined($self->pidentity()) || defined($self->nconserved()) || defined($self->evalue())){

	#We consider that all the filters are ok by default.
	my $ok = 1;
	
	if(defined($self->pidentity())){
	    my $pident = ($feat->get_tag_values('Note'))[0];
	    if(defined($pident)){
		$pident =~ s/.*Pident:([\d\.]+).*/$1/;
		chomp($pident);
		$ok = 0 if($pident <= $self->pidentity());
	    }else{
		$ok = 0;
	    }
	}
	if(defined($self->nconserved())){
	    my $ncons = ($feat->get_tag_values('Note'))[0];
	    if(defined($ncons)){
		$ncons =~ s/.*Ncons:([\d\.]+).*/$1/;
		chomp($ncons);
		$ok = 0 if($ncons <= $self->nconserved());
	    }else{
		$ok = 0;
	    }
	}
	if(defined($self->pcoverage())){
	    my $pcover =  ($feat->get_tag_values('Note'))[0];
	    if(defined($pcover)){
		$pcover =~ s/.*Pcover:([\d\.]+).*/$1/;
		chomp($pcover);
		$ok = 0 if($pcover <= $self->pcoverage());
	    }else{
		$ok = 0;
	    }
	}
	if(defined($self->evalue())){
	    my $fevalue =  ($feat->get_tag_values('Note'))[0];
	    if(defined($fevalue)){
		$fevalue =~ s/.*Evalue:([\+\-e\d\.]+).*/$1/;
		chomp($fevalue);

		my $bigfloat = Math::BigFloat->new($fevalue);
		$fevalue = $bigfloat->bstr();
		my $uevalue = $self->evalue();

		$ok = 0 if($fevalue <= $uevalue);

	    }else{
		$ok = 0;
	    }
	}

	return $ok;

    }else{
	return 1;
    }

}

=head1 initArraySeq

Description: Initailize an array with zeros corresponding to the sequence length
  Arguments: Length of the sequence
    Returns: array reference.

=cut 

sub initArraySeq ($) {

    my($self, $length) = @_;
    return unless($length);

    if($length > $MAX_SUBJ_LENGTH){
	my $msg = "HIT sequence is too large ($length). Maximum allowed for memory and speed performance is $MAX_SUBJ_LENGTH.\n";
	$msg   .= "\tEither switch off the filter(s) or run your blast in the opposite way, switch query <=> subject\n";

	return (undef, $msg);
    }

    my $ref_array = [ ];
    $ref_array->[$_] = 0 for(0..$length-1);

    return $ref_array;

}

=head1 fillGapsForCoverage

Description: This method is used to fill the gaps we can have in 2 regions having a match between sequences.
             This allow to get the total length of the coverage by switching empty array cell (0 -> 1), then
             when we get the statitistics (sub getStats) we know the coverage of the subject.
  Arguments: $hit, hash reference.
    Returns: 1 on error
             0 on success

=cut

sub fillGapsForCoverage ($){

    my($self, $hits) = @_;

    return(1) unless($hits && ref($hits) eq 'HASH');

    #Get the length of the hit (subject)
    my $length = $hits->{length};

    #Get the array representing the sequence
    my $sequence = $hits->{sequence};

    return unless(defined($hits->{hsps}));

    foreach my $hsp (@{$hits->{hsps}}){
	
	my $start = $hsp->{rstarts} - 1;
	my $end   = $hsp->{rends}   - 1;

	for(my $i = $start; $i <= $end; $i++){

	    $sequence->[$i] = 1 unless($sequence->[$i] > 0); #Avoid to overwrite identity or conserved

	}

    }

    return 0;

}

=head1 _BP_fillGapsForCoverage

Description: This method is used to fill the gaps we can have in 2 regions having a match between sequences.
             This allow to get the total length of the coverage by switching empty array cell (0 -> 1), then
             when we get the statitistics (sub getStats) we know the coverage of the subject.
  Arguments: $hits, hash reference.
             $BioSeqF array reference of Bio::SeqFeature::Generic objects
    Returns: 1 on error
             0 on success

=cut

sub _BP_fillGapsForCoverage ($$){

    my($self, $hits, $BioSeqF) = @_;

    return(1) unless($hits && ref($hits) eq 'HASH');

    return(1) unless($BioSeqF->[0]->isa("Bio::SeqFeature::Generic"));

    my $sequence = $hits->{sequence};


    foreach my $hsp (@$BioSeqF){
      if($hsp->get_tag_values('Starts') && $hsp->get_tag_values('Ends')){	
	#Get the start and stop from the tag values.
	my $start = ($hsp->get_tag_values('Starts'))[0] - 1  ;
	my $end   = ($hsp->get_tag_values('Ends'))[0]  - 1  ;


	for(my $i = $start; $i <= $end; $i++){

	    $sequence->[$i] = 1 unless($sequence->[$i] > 0); #Avoid to overwrite identity or conserved

	}

      }
    }
    return 0;

}

=head1 getStats

Description: This method is used to calculate the coverage, identity and conserved percentage
  Arguments: Hit reference hash table
    Returns: Coverage, Identity, Conservation percentages

=cut

sub getStats ($;$$) {

    my($self, $hits, $BioSeqF) = @_;

    #We use the BioPerl objects to store the objects.
    if(ref($BioSeqF) eq 'ARRAY' && $BioSeqF->[0]->isa("Bio::SeqFeature::Generic")){
	return if($self->_BP_fillGapsForCoverage($hits, $BioSeqF));
    }
    else{
	return unless($hits && ref($hits) eq 'HASH');
	return if($self->fillGapsForCoverage($hits));
    }

    my $length = $hits->{length};

    return unless($length);

    my $coverage  = grep { !/0/  } @{$hits->{sequence}};
    my $identity  = grep { /[3]/ } @{$hits->{sequence}};
    my $conserved = grep { /[2]/ } @{$hits->{sequence}};

    $conserved += $identity if($self->appl() =~ /blastx/i);

#    print STDERR "getStats: Length = $length - ", $hits->{length}, " Coverage = $coverage Id = $identity [$id] Cons = $conserved [$cons]\n";

    return(($coverage / $length) * 100, ($identity / $coverage) * 100, ($conserved / $coverage) * 100);

}


=head2 printRegion

Descrption: Get/Set method to know if need to print the Region part of GFF3.
Argument: $printRegion (boolean)
Return: boolean

=cut

sub printRegion ($){

    my($self, $print) = @_;

    if(defined($print)){
	$self->{_printregion} = $print;
    }

    return $self->{_printregion};
}

=head1 PRIVATE METHODS

=cut

=head1 _init

Descrption: When new methods passes arguments, removes first '-'
 Arguments: Reference to an array with parameters
   Returns: Array

=cut

sub _init ($){

  my($self, $array) = @_;

  return unless($array && ref($array) eq 'ARRAY');

  $self->{_topseen} = 0;
  $self->{_printregion} = 1;

  return @$array;
}

=head1 _topSeen

Description: Get/Set method. Checks whether we need to print ##gff-version 3 or not.
   Argument: $seen 
    Returns: value

=cut

sub _topSeen (;$){

  my($self, $seen) = @_;

  if(defined($seen)){
      $self->{_topseen} = $seen;
  }
  return $self->{_topseen};

}


=head1 _checkApplication

Description: check if the application passed as an argument is supported
Arguemtns: Application name
Returns : 1 if supported
          0 otherwise

=cut

sub _checkApplication ($) {

    my($self, $appl) = @_;

    return 0 unless $appl;

    for(keys %$ApplInfo2GFFFeat){
	if(lc($appl) eq lc($_)){
	    #As the name of the program is stored in $ApplInfo2GFFFeat, we automatically set this name to write GFF3.
	    $self->appl($_);
	    return 1;
	}
    }

    return 0;
}

=head1 _print_debug

 Description: prints debugging message
   Arguments: message to print
     Returns:

=cut

sub _print_debug (;$){

    my($self, $mess) = @_;

    if($mess){
	print STDERR "[DEBUG] $mess\n" if($self->debug());
    }

    return;

}

=head1 _exitOnError

 Description: Prints a message on STDERR and die.
   Arguments: Message to print
     Returns:

=cut

sub _exitOnError (;$){

    my ($self, $mess) = @_;

    if($mess){
	print STDERR "[ERROR] $mess\n";
    }

    die;
}

=head1 _print_verbose

 Description: Prints verbose message.
   Arguments: Message to print
     Returns:

=cut

sub _print_verbose (;$){

    my ($self, $mess) = @_;

    if($mess){
	print STDERR "[VERBOSE] $mess\n" if($self->verbose());
    }

    return;
}


#
#Empty subroutine
#
sub DESTROY {
    #do nothing
    #it avoids warning from AUTOLOAD
}


#
# Allow the user to by pass use of getParam, use name of the parameter directlry as a function name.
# The user can also set a value by giving an argument
#

sub AUTOLOAD {
    my ($self, $arg) = @_;

    $AUTOLOAD =~ /.*::(\w+)/ || return(0 , "No such method implemented ($AUTOLOAD).");
    my $attr = $1;

    my $hash = $self->getParams();

    if(exists($hash->{$attr})){
	if(defined($arg)){
	    $self->setParam($attr, $arg);
	}
	return $self->getParam($attr);
    }else{
	my $msg = "\n\t- " . join("\n\t- ", keys %{$self->getParams()}) . "\n";
	die "\nCould not map and create function $AUTOLOAD\n\tAUTOLOAD methods allowed : $msg";
    }
}

1;
