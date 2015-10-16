#
# Bioperl module for EuGene
# Cared for by Fabrice Legeai <flegeai@infobiogen.fr> and Joelle Amselem <jamselem@infobiogen.fr>
#
# Copyright Fabrice Legeai, Joelle Amselem
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

EuGene - Results of one EuGene run

=head1 DESCRIPTION

The EuGene module provides a parser for EuGene gene structure prediction
output. It parses gene predictions into Bio::Tools::Prediction::Gene object 
(a Bio::SeqFeature::Gene::Transcript).

This module also implements the Bio::SeqAnalysisParserI interface, and thus
can be used wherever such an object fits. See L<Bio::SeqAnalysisParserI>.

=head1 SYNOPSIS

   use EuGene;

   #filename for  Eugene result and Bio::Seq object for  Eugene query sequence :
   $eugene = EuGene->new(-file => 'result.eugene', -query_seq => 'seq');

   # filehandle and Bio::Seq object for  Eugene query sequence :
   $eugene = EuGene->new( -fh  => \*INPUT , -query_seq => 'seq');

   # parse the results
   # note: this class is-a Bio::Tools::AnalysisResult which implements
   # Bio::SeqAnalysisParserI, i.e., $eugene->next_feature() is the same

   # while($pred = $eugene->next_prediction()) {
       # $pred is an instance of Bio::Tools::Prediction::Gene
       # which inherits off Bio::SeqFeature::Gene::Transcript.
       #  You can access to the following methods to get the predictions features

       # $pred->exons() returns an array of Bio::SeqFeature::Gene::Exon
       # all features:
       @all_features = $pred->exons();

       # initial cds only
       @init_cds = $pred->exons('Initial');
       # internal cds only
       @intrl_cds = $pred->exons('Internal');
       # terminal cds  only
       @term_cds = $pred->exons('Terminal');
       # singleton  only
       @single_cds = $pred->exons('Single');

       # utrs only
       @utrs = $pred->utrs();

       # utr3 or utr5 only according to the SOFA (Sequence Ontology)
       @utrs = $pred->utrs('three_prime_UTR');
       @utrs = $pred->utrs('five_prime_UTR');

       # exons (hybrids of coding and non-coding provided by utr and cds merging)
       @exons = $pred->exons('exon');
   }

   # essential if you gave a filename at initialization (otherwise the file
   # will stay open)
   $eugene->close();


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bio.perl.org/MailList.html     - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via email
or the web:

  bioperl-bugs@bio.perl.org
  http://bugzilla.bioperl.org/

=head1 AUTHORS - Fabrice Legeai, Joelle Amselem

Email flegeai@infobiogen.fr, jamselem@infobiogen.fr

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package EuGene;
use vars qw(@ISA);
use strict;
use Bio::Tools::AnalysisResult;
use Bio::Tools::Prediction::Gene;
use Bio::SeqFeature::Gene::UTR;
use Data::Dumper;

@ISA = qw(Bio::Tools::AnalysisResult);


my @authorized_types = qw(Utr3 Utr5 Init Term Intr Sngl);
my %ExonTags = (
		'Init' => 'Initial',
		'Intr' => 'Internal',
		'Term' => 'Terminal',
		'Sngl' => 'Single'
	       );

=head2 _initialize_state

 Title   : _initialize_state
 Usage   : n/a; usually called by _initialize() itself called by new()
 Function: This method is supposed to reset the state such that any 'history'
           is lost. State information that does not change during object
           lifetime is not considered as history, e.g. parent, name, etc shall
           not be reset. An inheriting object should only be concerned with
           state information it introduces itself, and for everything else
           call SUPER::_initialize_state(@args).

           The argument syntax is the same as for new() and _initialize(),
           i.e., named parameters following the -name=>$value convention.
           The following parameters are dealt with by the implementation
           provided here:
              -INPUT, -FH, -FILE
           (tags are case-insensitive).
 Example :
 Returns :
 Args    :

=cut

#--------------------
sub _initialize_state  {
#--------------------

  my ($self,@args,) = @_;
  # first call the inherited method!
  $self->SUPER::_initialize_state(@args);
  # our private state variables
  $self->{'_parsed'} = 0;
  # array of pre-parsed predictions
  $self->{'_preds'} = [];
  # a Bio::Seq object (the Eugene query sequence)
  my ($query_seq) = $self->_rearrange([qw(QUERY_SEQ)], @args);
  $self->analysis_query($query_seq);
  $self->analysis_method('Eugene');
  return 1;
}

=head2 analysis_query

 Usage     : $eugene->analysis_query();
 Purpose   : returns the name query sequence
 Returns   : String
 Argument  : n/a

=cut

#------------------
sub analysis_query  {
#------------------
    my ($self, $obj) = @_;
    if($obj) {
	$self->{'_analysis_query'} = $obj;
    }
    return $self->{'_analysis_query'};

}
=head2 analysis_method

 Usage     : $eugene->analysis_method();
 Purpose   : Inherited method. Overridden to ensure that the name matches
             /eugene/i.
 Returns   : String
 Argument  : n/a

=cut

#-------------
sub analysis_method {
#-------------
    my ($self, $method) = @_;
    if($method && ($method !~ /eugene/i)) {
	$self->throw("method $method not supported in " . ref($self));
    }
    return $self->SUPER::analysis_method($method);
}


=head2 next_feature

 Title   : next_feature
 Usage   : while($gene = $eugene->next_feature()) {
                  # do something
           }
 Function: Returns the next gene prediction of the EuGene result
           file. Call this method repeatedly until FALSE is returned.

 Example :
 Returns : Bio::Tools::Prediction::Gene
 Args    :

=cut

#-------------
sub next_feature {
#-------------

  my ($self,@args) = @_;
  # even though next_prediction doesn't expect any args (and this method
  # does neither), we pass on args in order to be prepared if this changes
  # ever
  return $self->next_structure(@args);
}


=head2 next_prediction

 Title   : next_prediction
 Usage   : while($gene = $eugene->next_prediction()) {
                  # do something
           }
 Function: Returns the next gene prediction of the EuGene result
           file. Call this method repeatedly until FALSE is returned.

 Example :
 Returns : Bio::Tools::Prediction::Gene
 Args    :

=cut

#-------------
sub next_prediction {
#-------------
  my ($self) = @_;
  unless ($self->_parsed) {
    #print " not yet parsed\n";           #debug
    $self->_parse;
    #print " all is parsed\n";        #debug
  }

  my $pred = shift(@{$self->{'_preds'}});
  if($pred){
    $pred=$self->_merge_in_exon($pred);

#    Execute in the script calling Eugene.pm
#    $pred->add_tag_value('primary_transcript',$pred->seq()->seq());
#    my @cds = $pred->exons();
#    my @cds_sorted = sort{$a->start <=> $b->start}  @cds;
#    my $cds_seq="";
#    foreach my $exon (@cds_sorted) {
#      if($exon->is_coding()) {
#	$cds_seq .= $exon->seq()->seq;
#      }
#    }
#    $pred->add_tag_value('cds_sequence',$cds_seq);

  }
  #print " \n next prediction called \n";      #debug
  return  $pred;
}


=head2 _parse

 Title   : _parse
 Usage   : private method
           # do something

 Function: this method is called by next_prediction; it completely parses Eugene results file to generate a list of Bio::SeqFeature::Gene::GeneStructure, available through the structure method

 Example :
 Returns : 1 if success 0 otherwise
 Args    : none

=cut


#-------------
sub _parse {
#-------------
  my $self=shift;

  while (my $line = $self->_readline()) {
    # if parse the feature
    $self->_parse_feature($line);

  }

  # stores the last prediction
 # $self->_add_pred($self->_current_pred);
  $self->_parsed('TRUE');
}


=head2 _parse_feature

 Title   : _parse_feature
 Usage   : private method
           # do something

 Function: this method parse a eugene result line
 Example :
 Returns : and returns a Bio::SeqFeature::Generic
 Args    : a string

=cut

#------------
sub _parse_feature {
#------------
  my $self=shift;
  my $line = shift;

  # parse the line 
  #seq.1.1       Init    +        927    1014       88     +1      +3     926    1015     0.0
  my ($feat_tag, $type, $s, $lend, $rend, $length, $phase, $frame, $ac, $do, $pr) = split /\s+/, $line;

  $frame =~ s/NA/\./;
  $frame =~ s/(\+|\-)//;
  my $feat=undef;
  unless (grep {$_ eq $type} @authorized_types) { return undef}

  # split $feat_tag in transcript_name and feat_id

  my ($egn_transcript_name,$feat_id) = ($feat_tag =~ /(^.*)\.(\d+)$/);

  # transcript_id will be used as mRNA tag ID and other feats tag Parent

  my ($transcript_id) = ($egn_transcript_name =~ /^.*\.(\d+)$/);
  $transcript_id=sprintf "%2.2d", $transcript_id;
  $transcript_id = $transcript_id ;
  
  #  my $transcript_name = $self->analysis_query->display_name.'_'.$self->analysis_method.'_'.$transcript_id;
  
  my $transcript_name = $self->analysis_query->display_name.'.'.$transcript_id;
  if (defined $self->_current_pred) {
    if ($transcript_name ne  $self->_current_pred->_tag_value('Name')){
      
#      my $ID=$self->analysis_query->display_name.'_'.$self->analysis_method.'_'.$transcript_id;
      my $ID=$self->analysis_query->display_name.'.'.$transcript_id;
      my $pred = Bio::Tools::Prediction::Gene->new(
						   -seq_id => $self->analysis_query->display_name,
						   -primary => 'mRNA',
						   -source => $self->analysis_method,
						   -tag => {
							    ID => $ID,
							    Name => $transcript_name
							   }
						  );
      $pred->attach_seq($self->analysis_query());
      $self->_add_pred($pred);
    }
  }
  else {
 #   my $ID=$self->analysis_query->display_name.'_'.$self->analysis_method.'_'.$transcript_id;
    my $ID=$self->analysis_query->display_name.'.'.$transcript_id;
    my $pred = Bio::Tools::Prediction::Gene->new(
						 -seq_id => $self->analysis_query->display_name,
						 -primary => 'mRNA',
						 -source =>$self->analysis_method,
						 -tag => {
							  ID => $ID,
							  Name => $transcript_name
							 }
						);
    $pred->attach_seq($self->analysis_query());
    $self->_add_pred($pred);
  }
  my $pred = $self->_current_pred;
  $pred->display_name($transcript_name);
  {
    #print "   parse next feature ",$type,"\n";         #debug
    #-------- utr
    if ($type =~ /Utr/) {
      if($type =~ /Utr5/){$type = "5'-UTR";}
      else{$type = "3'-UTR";}
      my @parentIDs=$pred->get_tag_values('ID');
      my $parentID=$parentIDs[0];
      my $strNumber=sprintf "%2.2d",scalar($pred->utrs($type)) + 1 ;
      my $ID=$parentID;
      $feat = Bio::SeqFeature::Gene::UTR->new(
					      -seq_id => $self->analysis_query->display_name,
					      -verbose => -1,
					      -primary_tag => $type,
					      -is_coding => 0,
					      -start => $lend,
					      -end => $rend,
					      -strand => $s,
					      -frame => '.',
					      -score => $pr,
					      -source =>$self->analysis_method,
					      -tag => {
						       ID=> $ID,
						       Parent => $pred->get_tag_values('ID')
						      }
					     );
      $pred->add_utr($feat);
      last;
    } # endif type Utr

    # otherwise it is a CDS (Initial, Internal, Terminal or Single)---------
    my @parentIDs=$pred->get_tag_values('ID');
    my $parentID=$parentIDs[0];
    my $strNumber=sprintf "%2.2d", $feat_id;
    my $ID=$parentID;#.'_'.$type.'_'.'cds_'.$strNumber;
    #my $ID=$self->analysis_query->display_name.'_'.$self->analysis_method.'_'.'cds_'.$feat_id;
    $feat = Bio::SeqFeature::Gene::Exon->new(
					     -seq_id => $self->analysis_query->display_name,
					     -primary => $ExonTags{$type},
					     -is_coding => 1,
					     -start => $lend,
					     -end => $rend,
					     -strand => $s,
					     -frame => '.',
					     -score => $pr,
					     -source => $self->analysis_method,
					     -tag => {
						      ID => $ID,
						      Parent => $pred->get_tag_values('ID'),
						     }
					    );
    #	Name => 'protein_'.$transcript_name
    $pred->add_exon($feat, $ExonTags{$type});

  } #endof cases
  #print Dumper $feat,"\n";
  #return($feat);
}

=head2 _merge_in_exon

 Title   : _merge_in_exon
 Usage   : private method

 Function: this method merge utr and cds if necessary to generate exons and add them to the prediction

 Example : _merge_in_exon($pred)
 Returns : a Bio::Tools::Prediction::Gene
 Args    : a Bio::Tools::Prediction::Gene

=cut


#----------
sub _merge_in_exon {
  my $self=shift;
  my $pred=shift;

  my @features = $pred->exons('');
  my @features_sorted = sort{$a->start <=> $b->start}  @features;

  my $exon_feat;
  for(my $i=0;$i<(scalar(@features_sorted)-1);$i++){
    my $j=$i;
    while($features_sorted[$j+1] && $features_sorted[$j]->end+1 eq $features_sorted[$j+1]->start){
      $j++;
    }
    my @parentIDs=$pred->get_tag_values('ID');
    my $parentID=$parentIDs[0];
    my $strNumber=sprintf "%2.2d", scalar($pred->exons('exon')) + 1;
    my $ID=$parentID;#.'_'."exon_".$strNumber;
    #my $ID=$self->analysis_query->display_name.'_'.$self->analysis_method.'_'."exon_".scalar($pred->exons('exon')+1);
    $exon_feat = Bio::SeqFeature::Gene::Exon->new(
						  -seq_id => $self->analysis_query->display_name,
						  -primary => 'exon',
						  -is_coding => 0,
						  -start => $features_sorted[$i]->start,
						  -end => $features_sorted[$j]->end,
						  -strand => $features_sorted[$i]->strand,
						  -frame => '.',
						  -score => '.',
						  -source => $self->analysis_method,
						  -tag => {
							   ID => $ID,
							   Parent => $features_sorted[$i]->get_tag_values('Parent')
							  }
						 );
    $pred->add_exon($exon_feat, 'exon');
    $i=$j;
  }
  return($pred);
}

=head2  _current_pred
 Title   : _current_pred 
 Usage   : private method
           # do something

 Function: this method returns the Bio::SeqFeature::TranscriptI currently parsed. If arguments the current_transcript replaced.

 Example : 
 Returns : a Bio::Tools::Prediction::Gene
 Args    : a Bio::Tools::Prediction::Gene  or nothing

=cut


#----------
sub _current_pred {
#----------
  my $self=shift;
  if (@_) {
    $self->{'current_pred'}= shift;
  }
  return $self->{'current_pred'};
}


=head2 _add_pred

 Title   : _add_pred
 Usage   : private method
           # do something

 Function: adds a prediction to the predictions list and put it as the current one

 Example : 
 Returns : 
 Args    : a Bio::Tools::Prediction::Gene

=cut

#------------
sub _add_pred {
#------------
  my $self=shift;
  my $pred = shift;

  if (defined $self->{'_preds'}) {
    push  @{$self->{'_preds'}}, $pred;
  }
  else {
    $self->{'_preds'}=[$pred];
  }
  $self->_current_pred($pred);
} 


=head2  _parsed

 Title   : _parsed
 Usage   : private method
           # do something

 Function: if the  file has been parsed this  method returns 1 0 otherwise. If  argument is the string  'TRUE', this method returns 1 0 otherwise

 Example : 
 Returns : 1 if results have been already parsed, 0 otherwise
 Args    : TRUE or nothing

=cut


#----------
sub _parsed {
#----------
  my $self=shift;

  if (@_) {
    my $arg = shift;
    if ($arg eq 'TRUE') {
      $self->{'parsed'}=1;
      return(1);
    }
    else {
      $self->throw("Error the function _parsed should receive the string 'TRUE' or nothing");
    }
  }
  if (defined $self->{'parsed'}) {
    return       $self->{'parsed'};
  }
}
