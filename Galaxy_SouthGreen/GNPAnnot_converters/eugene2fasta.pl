#!/usr/bin/perl

use FindBin;                 
use lib "$FindBin::Bin/lib";
use Bio::SearchIO;
use EuGene;
use Bio::SeqIO;
use Bio::Location::Split;
use Bio::Location::Simple;
use Bio::SeqFeature::Generic;
use XML::Simple qw(:strict);
use EuGene;
use File::Basename;
use Data::Dumper;
use Bio::Tools::GFF;


my $file_fasta  = shift;
my $file_eugene = shift;
my $output_gff 	= shift;
my $output_embl = shift;
my $output_embl_utr = shift;
my $output_gene     = shift;
my $output_cds      = shift;
my $output_pep      = shift;
my $output_region   = shift;
my $version         = shift;

&eugene2fasta($file_fasta,$file_eugene,$output_gff,$output_embl,$output_embl_utr,$output_gene,$output_cds,$output_pep,$output_region);

sub eugene2fasta {
    my ($file_fasta,$file_eugene, $output_gff,$output_embl,$output_embl_utr,$output_gene,$output_cds,$output_pep,$output_region) = @_;
    my $seq    = Bio::SeqIO->new(
		-file   => $file_fasta,
		-format => 'fasta' 
	)->next_seq();
    my $eugene = EuGene->new(
		-file => $file_eugene,
		-query_seq => $seq
	);
    my $clone_name    = $seq->display_id();	
    my $gffio         = new Bio::Tools::GFF(
		-file        => ">$output_gff",
		-gff_version => 3
	);
    my $out_genomic_seq = Bio::SeqIO->new(
		-file   => ">$output_gene",
		-format => 'fasta' 
	);
    my $out_cds_seq     = Bio::SeqIO->new(
		-file   => ">$output_cds",
		-format => 'fasta' 
	);
    my $out_prot_seq    = Bio::SeqIO->new(
		-file   => ">$output_pep",
		-format => 'fasta' 
	);
    my $emblio = Bio::SeqIO->new(
		-file   => ">$output_embl",
		-format => "EMBL"
	);
    my $source_tag    = "eugene"; 
    my $species       = "rice";
    my $inference     = "{Eugene rice 3.2}";
    my $CDSnb         = 1;
    my $date          = `date`;
    open(UTR,">$output_embl_utr");
    chomp($date);
    my $feature = new Bio::SeqFeature::Generic(
		-primary_tag => 'source',
		-start       => 1,
		-end         => $seq->length()
	);
    $seq->add_SeqFeature($feature);
    my %cds;
    my %gene;
    my %mrna;
    my %utr;
    my %exon;
    my %utr3;       
    my %utr5;
    my ($gene_id,$mrna_id,$polypeptide_id);
	while ( my $prediction = $eugene->next_prediction() ) {
		my @cds = $prediction->exons();
		my @cds_sorted = sort { $a->start <=> $b->start } @cds;
		my @cds_coding;
		foreach my $exon (@cds_sorted) {
			if ( $exon->is_coding() ) {
				push @cds_coding, $exon;
			}
		}
		my $nb_exon = scalar(@cds_coding);
		my ( $query, $source, $type, $start, $end, $score, $strand, $frame,$id ) = split( "\t", $prediction->gff_string() );
		if ( $nb_exon > 0 ) {
			my $gene_id = sprintf( "%s_g%04d", $clone_name, $CDSnb * 10 );
			my $mrna_id = sprintf( "%s_t%04d", $clone_name, $CDSnb * 10 );
			my $polypeptide_id = sprintf( "%s_p%04d", $clone_name, $CDSnb * 10 );
			my $nb_cds   = 0;
			my $location = new Bio::Location::Split();
			my $seq_cds;
			my ( @start, @end );
			foreach my $exon (@cds) { 
				if ( $exon->is_coding() ) {
					$exon->primary_tag('CDS');
				}
				$exon->remove_tag('ID');
				my ( $query, $source, $type,  $start, $end, $score, $strand, $frame, $id) = split( "\t", $exon->gff_string() );
			 
				my $strand_embl = $strand eq "+" ? "+1" : "-1";
				if ( $type eq "CDS" ) {
					push @{ $cds{$gene_id}{start} },  $start;
					push @{ $cds{$gene_id}{end} },    $end;
					push @{ $cds{$gene_id}{strand} }, $strand; 
					push @{ $cds{$gene_id}{frame} }, $frame; 
				}
				if ( $type eq "exon" ) {
					push @{ $exon{$gene_id}{start} }, $start;
					push @{ $exon{$gene_id}{end} },   $end;
				}
				if ( $type eq "utr3prime" ) {
					push @{ $utr3{$gene_id}{start} }, $start;
					push @{ $utr3{$gene_id}{end} },   $end;
					if ( $strand eq "+" ) {
						print UTR "FT   3'UTR	    $start..$end\nFT		       /note=\"$gene_id\"\n";
						print UTR "FT		    /color=103\n";
					}
					else {
						print UTR "FT   3'UTR	    complement($start..$end)\nFT		   /note=\"$gene_id\"\n";
						print UTR "FT		    /color=103\n";
					}
				}
				if ( $type eq "utr5prime" ) {
					push @{ $utr5{$gene_id}{start} }, $start;
					push @{ $utr5{$gene_id}{end} },   $end;
					if ( $strand eq "+" ) {
						print UTR "FT   5'UTR	    $start..$end\nFT		       /note=\"$gene_id\"\n";
						print UTR "FT		    /color=103\n";
					}
					else {
						print UTR "FT   5'UTR	    complement($start..$end)\nFT		   /note=\"$gene_id\"\n";
						print UTR "FT		    /color=103\n";
					}
				}
			}
		}
		$CDSnb++;
	}
	my $cpt_gene = 1;
    foreach my $gene_id ( sort keys %cds ) {
		$mrna_id        = sprintf( "%s_t%04d", $clone_name, $cpt_gene * 10 );
		$polypeptide_id = sprintf( "%s_p%04d", $clone_name, $cpt_gene * 10 );
		$cpt_gene++;
		my @cds_start   = @{ $cds{$gene_id}{start} };
		my @cds_end     = @{ $cds{$gene_id}{end} };
		my @cds_strand  = @{ $cds{$gene_id}{strand} };
		my @cds_frame  = @{ $cds{$gene_id}{frame} };
		my @exon_start  = @{ $exon{$gene_id}{start} };
		my @exon_end    = @{ $exon{$gene_id}{end} };
		my @utr5_start  = @{ $utr5{$gene_id}{start} };
		my @utr5_end    = @{ $utr5{$gene_id}{end} };
		my @utr3_start  = @{ $utr3{$gene_id}{start} };
		my @utr3_end    = @{ $utr3{$gene_id}{end} };
		my $strand_gene = $cds_strand[0];
		my $strand_mrna = $strand_gene; 
		my $start_gene;
		my $end_gene;
		my $start_mrna = $cds_start[0];
		my $end_mrna   = $cds_end[-1];

		if ( $strand_gene eq "-" ) {
			if (@utr3_start) {
				$start_gene = _min(@utr3_start);
			}
			else {
				$start_gene = $exon_start[0];
			}
			if (@utr5_start) {
				$end_gene = _max(@utr5_end);
			}
			else {
				$end_gene = $exon_end[-1];
			}
		}
		else {
			if (@utr3_start) {
				$start_gene = _min(@utr5_start);
			}
			else {
				$start_gene = $exon_start[0];
			}
			if (@utr5_start) {
				$end_gene = _max(@utr3_end);
			}
			else {
				$end_gene = $exon_end[-1];
			}
		} 
        my $seq_cds;
		my $location = new Bio::Location::Split();
		my $strand1;
		my $min_start;
		my $max_end; 
		for ( my $i = 0 ; $i <= $#cds_start ; $i++ ) { 
			$location->add_sub_Location(
				new Bio::Location::Simple(
					-start  => $cds_start[$i],
					-end    => $cds_end[$i],
					-strand => $strand_gene
				)
			); 
			
			$seq_cds .= $seq->subseq( $cds_start[$i], $cds_end[$i] ); 
			$seq_cds=~s/X/N/g;
		}
		my $location1 = new Bio::Location::Split();
		for ( my $i = 0 ; $i <= $#exon_start ; $i++ ) {
			$location1->add_sub_Location(
				new Bio::Location::Simple(
					-start  => $exon_start[$i],
					-end    => $exon_end[$i],
					-strand => $strand_gene
				)
			);
		}
		my $strand_embl = $strand_gene eq "+" ? "+1" : "-1";
		$location1->guide_strand($strand_embl);
		$location->guide_strand($strand_embl);
		my $seq_genomic = $seq->subseq( $location1->start, $location1->end );
		$seq_genomic=~s/X/N/g;
        my $seq_obj_genomic = Bio::PrimarySeq->new (
			-display_id  => $gene_id, 
			-seq         => $seq_genomic
		);
        my $seq_obj_cds = Bio::PrimarySeq->new(
			-display_id => $mrna_id,
			-desc       => $polypeptide_id,
			-seq        => $seq_cds
		);
		if ( $strand_gene eq "-" ) {
			$seq_obj_cds     = $seq_obj_cds->revcom();
			$seq_obj_genomic = $seq_obj_genomic->revcom();
		} 
		my $seq_obj_prot = $seq_obj_cds->translate();  
 		$out_genomic_seq->write_seq($seq_obj_genomic);
		$out_cds_seq->write_seq($seq_obj_cds);
		$seq_obj_prot->display_id($polypeptide_id); 
		$out_prot_seq->write_seq($seq_obj_prot);
		my $length_aa = $seq_obj_prot->length();
		my $locstr    = &encod( $location->to_FTstring() );
		my %label     = (
			locus_tag            => $gene_id,
			inference            => $inference,
			owner                => $self->{species},
			date                 => $date,
			original_location    => $locstr,
			annotator_comment    => "missing_annotator_comment",
			alternative_location => "no_alternative_location",
			Ontology_term        => "CC_status:in_progress",
			Ontology_term        => "CC_evidence:automatic",
			Ontology_term        => "CC_evidence_code:ISS"
		); 
        my $feature_gene_manual =  new Bio::SeqFeature::Generic(
			-seq_id      => $clone_name,
			-source_tag  => "manual_curation",
			-primary_tag => 'gene',
			-start       => $start_gene,
			-end         => $end_gene,
			-strand      => $strand_gene,
			-tag         =>  {
				ID        => $gene_id,
				Name      => $gene_id,
				locus_tag => $gene_id
			}
		);
        my $feature_mrna_manual =  new Bio::SeqFeature::Generic(
			-seq_id      => $clone_name,
			-source_tag  => "manual_curation",
			-primary_tag => 'mRNA',
			-start       => $start_gene,
			-end         => $end_gene,
			-strand      => $strand_mrna,
			-tag         =>  {
				ID        => $mrna_id,
				Name      => $gene_id,
				Parent    => $gene_id,
				locus_tag => $mrna_id
			}
		);
        my $feature_polypeptide_manual =  new Bio::SeqFeature::Generic(
			-seq_id      => $clone_name,
			-source_tag  => "manual_curation",
			-primary_tag => 'polypeptide',
			-start       => $cds_start[0],
			-end         => $cds_end[-1],
			-strand      => $strand_gene,
			-tag         =>  {	
				ID           => $polypeptide_id,
				Name         => $polypeptide_id,
				Derives_from => $mrna_id,
				%label
			}
		);	
		my $feature_exon_manual = new Bio::SeqFeature::Generic(
			-seq_id      => $clone_name,
			-source_tag  => "manual_curation",
			-primary_tag => 'exon',
			-location    => $location1, 
			-tag         =>  {
				Parent => $mrna_id
			}
		);	
		my $feature_cds_manual = new Bio::SeqFeature::Generic(
			-seq_id      => $clone_name,
			-source_tag  => "manual_curation",
			-primary_tag => 'CDS',
			-location    => $location, 
			-tag         =>  {
				Parent => $mrna_id
			}
		);	
        my $feature_cds_embl = new Bio::SeqFeature::Generic(
			-primary_tag => 'CDS',
			-location    => $location, 
			-tag         =>  {%label }
		);
        my $feature_mrna_embl =  new Bio::SeqFeature::Generic(
			-primary_tag => 'mRNA',
			-location    => $location,
			-tag         =>  {
				ID        => $mrna_id,
				Name      => $gene_id,
				Parent    => $gene_id,
				locus_tag => $mrna_id
			}
		);
        $seq->add_SeqFeature($feature_gene_manual);
        #$seq->add_SeqFeature($feature_mrna_embl);
        $seq->add_SeqFeature($feature_cds_embl);
		$gffio->write_feature($feature_gene_manual);
		$gffio->write_feature($feature_mrna_manual);
        $gffio->write_feature($feature_polypeptide_manual);
		$gffio->write_feature($feature_exon_manual);	
		$gffio->write_feature($feature_cds_manual);		
		 
    }
    $emblio->write_seq($seq);
    $gffio->close();
    $eugene->close();
    $emblio->close();
    close UTR;
    $output_region = &embl2fna($output_embl,$seq->length,$output_region,$clone_name);  
}


sub _max (@) {
    my $i = shift;
    foreach (@_) {
	$i = $_ if $_ > $i;
    }
    return $i;
}


sub _min (@) {
    my $i = shift;
    foreach (@_) {
	$i = $_ if $_ < $i;
    }
    return $i;
}


sub encod {
    my $encod = shift; 
    $encod =~ s/([^a-zA-Z0-9_. :?^*\(\)\[\]@!-])/uc sprintf("%%%02x",ord($1))/eg;
#    $encod =~ s/([,;=%&()'"])/sprintf "%%%02X", ord($1)/gei;
    return $encod;
}


sub embl2fna {
    my ($file,$end_bac,$file_fna,$clone_name) = @_;
    my $source_tag = "eugene";
    my $seq = Bio::SeqIO->new(
		-format => 'EMBL',
		-file   => $file
	)->next_seq();
    my $overlap = 5000;
    my $start_bac = 1;
    my @features;
    my $CDSnb = 0;  
    my $seqobj   = Bio::SeqIO->new(
		-file   => ">$file_fna",
		-format => 'fasta' 
	);
    foreach my $feat ($seq->all_SeqFeatures())  {
        my $primary_tag = $feat->primary_tag();
        if ($primary_tag eq "CDS") {
            push @features,$feat;
            $CDSnb ++;
        }
    }
    my $nb_gene = scalar(@features);
    if ($nb_gene > 1) {
        for (my $i=0;$i<=$#features;$i++) {
            my ($start_of_region,$end_of_region);
            my $start = $features[$i]->start;
            my $end = $features[$i]->end ;
            my $next = $i + 1;
            my $locus_tag = ($features[$i]->get_tag_values('locus_tag'))[0];
            if ($i == 0) {
                my $next_start = $features[$next]->start;
                $start_of_region = $start - $overlap < 0 ? "1" : $start - $overlap + 1;
                $end_of_region = $next_start - $end - 1 < $overlap ? $end + $next_start - $end - 1 : $end + $overlap;
            } 
	    	elsif ($next == $#features + 1) {
                my $previous = $i - 1;
                my $previous_end = $features[$previous]->end;
                $start_of_region = $start - $previous_end < $overlap ? $previous_end + 1 : $start - $overlap + 1;
                $end_of_region = $end_bac - $end - 1 < $overlap ? $end + $end_bac - $end  : $end + $overlap;
            } 
	    	else {
                my $previous = $i - 1;
                my $next_start = $features[$next]->start;
                my $previous_end = $features[$previous]->end;
                $start_of_region = $start - $previous_end < $overlap ? $previous_end + 1 : $start - $overlap;
                $end_of_region = $next_start - $end - 1 < $overlap ? $end + $next_start - $end - 1 : $end + $overlap;
            }
            my $query = join("-",$locus_tag , $start_of_region ,$end_of_region);
			my $query = join( "-", $locus_tag,$clone_name,$start_of_region, $end_of_region );
            my $seq_region     = $seq->subseq($start_of_region,$end_of_region);
            my $seq_region_obj = Bio::PrimarySeq->new (
				-display_id  => $query,                                                     
				-seq         => $seq_region
			);
            $seqobj->write_seq($seq_region_obj);
        } 
    }
    else {
        my ($start_of_region,$end_of_region);
        my $start = $features[0]->start;
        my $end   = $features[0]->end ;
        my $locus_tag = ($features[0]->get_tag_values('locus_tag'))[0];
        $start_of_region = $start - $overlap < 0 ? "1" : $start - $overlap + 1;
        $end_of_region   = $end + $overlap > $end_bac ? $end_bac : $end+$overlap;	
		my $query = join( "-", $locus_tag, $clone_name,$start_of_region, $end_of_region );
        my $seq_region     = $seq->subseq($start_of_region,$end_of_region); 
        my $seq_region_obj = Bio::PrimarySeq->new (
			-display_id  => $query,
			-seq         => $seq_region
		);
        $seqobj->write_seq($seq_region_obj);
    }
    return $file_fna;
}


1;
