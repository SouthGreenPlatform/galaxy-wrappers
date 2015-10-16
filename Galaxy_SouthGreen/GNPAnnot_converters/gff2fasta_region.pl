#!/usr/bin/perl 
use Carp qw (cluck confess croak); 
use Pod::Usage;
use Error;
use Fatal qw (open close); 
use Getopt::Long;
use Bio::Tools::GFF;
use Bio::SeqIO;

 

my $file_fasta = shift;
my $file_gff = shift;
my $file_fasta_region = shift;
my $file_cds = shift;
my $file_prot_seq = shift;
  
my $seqin = new Bio::SeqIO(
	-file   => $file_fasta,
	-format => 'fasta'
);
my $gffin = new Bio::Tools::GFF(
	-file   	 => $file_gff,
	-gff_version => 3
);
 
my $seqout = Bio::SeqIO->new(
	-file   => ">$file_fasta_region",
	-format => 'fasta'
); 
my $out_cds_seq = Bio::SeqIO->new(
	-file   => ">$file_cds",
    -format => 'fasta'
);
my $out_pep_seq = Bio::SeqIO->new(
    -file   => ">$file_prot_seq",
    -format => 'fasta'
); 
my %sequence; 
while( my $seqobj = $seqin->next_seq) { 
	$sequence{$seqobj->display_id} = $seqobj;
}
$seqin->close; 

my %count;
my %gff;
while (my $feature = $gffin->next_feature){
	if ($feature->primary_tag() eq "gene") {	
		push @{$gff{$feature->seq_id}},$feature;
		$count{$feature->seq_id}++;
	}
	if ($feature->primary_tag() eq "CDS") {
		my ($gene_id) = $feature->get_tag_values("Parent");
		push @{$cds{$feature->seq_id}{$gene_id}}, $feature;
	}
}
$gffin->close; 

foreach my $seq_id (sort {$a cmp $b} keys %gff){ 
	if ($sequence{$seq_id}){
		my $seqobj = $sequence{$seq_id};
		my $end_chr = $seqobj->length(); 
		my $overlap   = 5000;
		my $num_gene = $count{$seq_id}; 
		foreach my $gene ( keys %{$cds{$seq_id}}){ 
			my $seq_cds;
			my $strand = $cds{$seq_id}{$gene}[0]->strand;  
			foreach my $feat (sort {$a->start<=>$b->start} @{$cds{$seq_id}{$gene}} ){
            	$seq_cds .= $seqobj->subseq( $feat->start,$feat->end );  
			} 
			my $seq_obj_cds = Bio::PrimarySeq->new(
				-display_id => $gene,
				-seq        => $seq_cds
			);
       		if ( $strand =~ /-/ ) {
            	$seq_obj_cds = $seq_obj_cds->revcom();
       		}
        	$out_cds_seq->write_seq($seq_obj_cds);
        	my $seq_obj_prot = $seq_obj_cds->translate(); 
        	$out_pep_seq->write_seq($seq_obj_prot);  
		}
		my @features = sort {$a->start<=>$b->start} @{$gff{$seq_id}};
		my $cpt = 0;
		for ( my $i = 0 ; $i <= $#features ; $i++ ) {
			$cpt++;
			my ( $start_of_region, $end_of_region );
			my $start     = $features[$i]->start;
			my $end       = $features[$i]->end;
			my $next      = $i + 1;
			my $locus_tag = ( $features[$i]->get_tag_values('Name') )[0];
			my $previous_start;
			my $previous_end;
			my $next_start;
			my $previous;
			if ( $i == 0 ) {
				$next_start = $features[$next]->start;
				$start_of_region = $start - $overlap < 0 ? "1" : $start - $overlap + 1;
				$end_of_region = $next_start - $end - 1 < $overlap ? $end + $next_start - $end - 1 : $end + $overlap;
			}
			elsif ( $next == $#features + 1 ) {
				$previous     = $i - 1;
				$previous_end = $features[$previous]->end;
				$start_of_region = $start - $previous_end < $overlap  ? $previous_end + 1 : $start - $overlap + 1;
				$end_of_region   = $end_chr - $end - 1 < $overlap ? $end + $end_chr - $end : $end + $overlap;
			}
			else {
				$previous     = $i - 1;
				$next_start   = $features[$next]->start;
				$previous_end = $features[$previous]->end;
				$start_of_region = $start - $previous_end < $overlap ? $previous_end + 1 : $start - $overlap;
				$end_of_region = $next_start - $end - 1 < $overlap ? $end + $next_start - $end - 1 : $end + $overlap;
			}
			my $query = join( "-", $locus_tag, $features[0]->seq_id,$start_of_region, $end_of_region );
			my $seq_region = $seqobj->subseq( $start_of_region, $end_of_region );
			my $seq_region_obj    = Bio::PrimarySeq->new(
				-display_id => $query,
				-seq        => $seq_region
			);
			$seqout->write_seq($seq_region_obj); 	
		}
	}	
}
$seqout->close;	

 
