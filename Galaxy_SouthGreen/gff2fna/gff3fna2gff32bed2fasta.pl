#!/usr/bin/perl

=pod

=head1 NAME

[gff3fna2gff32bed2fasta.pl - short description]

=head1 SYNOPSIS

    perl gff3fna2gff32bed2fasta.pl <gff3_input_file> <genome_sequence_file> <fna_output_file> <bed_output_file> -sort=on|off -bed=all|gene|mRNA|polypeptide|... -begin=0|-x|+x -end=0|-x|+x -flank=prom|term

=head1 REQUIRES

[Perl5.004, POSIX, Win32] #+++

=head1 DESCRIPTION

#--- What this script can do, how to use it, ...
[Script description]. #+++

=cut

use strict;
use Carp qw (cluck confess croak);
use warnings;
use Getopt::Long;
use Pod::Usage;
#--- when using try {...} catch XXX::YYY with {...} otherwise {...};
#+++ use Error qw(:try);
#--- when working with files
use Fatal qw/:void open close/;

use Bio::Tools::GFF;
use Bio::FeatureIO::bed;
use List::MoreUtils qw(uniq);
use Bio::SeqIO;
use Bio::Seq;
use Bio::DB::Fasta;
use JSON;
use Statistics::Descriptive;
use Bio::DB::SeqFeature;
use Env qw(HOME);
use Data::Dumper;

# Script global constants
##########################

=pod

=head1 CONSTANTS

B<$CONSTANT_NAME>: ([constant nature]) #+++

[constant description and use]. #+++

B<CONSTANT_NAME2>: (sub returning a [constant return nature]) #+++

[constant description and use]. #+++

#--- Example:
#--- B<PI>: (sub returning a real)
#---
#--- used by trigonometrical functions;
#--- can also be used to define an angle value.
#---
#--- ...
#---
#--- sub PI() { return (4 * atan2(1, 1)); }
#---

=cut

#+++ my [$CONSTANT_NAME] = ["value"];
#+++ sub [CONSTANT_NAME2]() { return ["value2"]; }


# Script global variables
##########################

=pod

=head1 VARIABLES

B<variable_name>: ([variable nature]) #+++

[variable description, use and default value]. #+++

#--- Example:
#--- B<$output_method>: (integer)
#---
#--- used to store current output method;
#--- 0=text on screen (default), 1=graphic on screen, 2=text in a file, 3=graphic to the printer.
#---
#---     ...
#---
#--- my $output_method = 0;
#---

=cut

my $gff3_file=$ARGV[0];  #gff3 input file
my $genome=$ARGV[1]; #genome sequence directory
my $fna_output_directory= $1 if $ARGV[2] =~ /(.*\/)[^\/]*\.?.*/i;
my $fna_output_file = $ARGV[2];
my $bed_output_directory= $1 if $ARGV[3] =~ /(.*\/)[^\/]*\.?.*/i;
my $bed_output_file = $ARGV[3];

my $title;
my $title2;
my $title3;

open my $FH_IN, "<", $gff3_file or die "Could not open file '$gff3_file' $!";
while (my $ligne = <$FH_IN>) {     # lit le fichier en entrée ligne par ligne
     if (($ligne !~ /^##/i) && ($ligne =~ /^[A-Z]{5}/i)){     
         my @mots = split('\t', $ligne);
         if (($mots[0] =~ /^([A-Z]{5})/i) && ($mots[1] !~ /manual_curation/i)){
        	$title = "$1-"."$mots[1]-sequence_feature" ;
         	$title2 = "$1-"."$mots[1]";
        	$title3 = "$1";
         }
     }     
}
close $FH_IN;

# Script global functions
##########################

=pod

=head1 FUNCTIONS

=head2 [subName] #+++

B<Description>: [function description]. #+++

B<ArgsCount>: [count of arguments] #+++

=over 4

=item variable_name: ([variable nature]) ([requirements]) #+++ see below

#--- requirement can be:
#--- (R)=required,
#--- (O)=optional
#--- (U)=optional and must be undef if omitted
[variable description]. #+++

=item variable_name2: ([variable nature]) ([requirements]) #+++

[variable description]. #+++

=back

B<Return>: ([return type]) #+++

[return description]. #+++

B<Exception>:

=over 4

[=item * exception_type:

description, case when it occurs.] #+++ see below

=back

#--- Example:
#---
#---=over 4
#---
#---=item * Range error:
#---
#---thrown when "argument x" is outside the supported range, for instance "-1".
#---
#---=back

B<Example>:

    [example code] #+++

=cut

#+++sub [subName] #+++
#+++{
#+++    my ([...]) = @_; #+++ add missing arguments
#--- if needed, parameters check:
#+++    # parameters check
#+++    if (0 != @_)
#+++    {
#+++        confess "usage: subName();"; #+++
#+++    }
#+++}
my ($flank, $begin, $end, $bed);
sub sortGff3{

	my $gff = new Bio::Tools::GFF(
	-file => $gff3_file,
	-gff_version => 3
	);
	
	my $cpt = 0;
	my (%feature, $id, $mrna_id, %gene_refseq, %mrna, %polypeptide, %other, %other_RNA, %cds, %transposable_gene, $transp_id);
	
	while(my $feature = $gff->next_feature) {
		if ($feature->primary_tag() eq "gene") {
			$cpt++;
			#last if $cpt ==20;
			($id) = $feature->get_tag_values("ID");	
			$feature{$id} = $feature;			
			push @{$gene_refseq{$feature->seq_id}} , $feature;
		}
		elsif ($feature->primary_tag() eq "pseudogene"){	
			push @{$gene_refseq{$feature->seq_id}} , $feature;		
		}
		elsif ($feature->primary_tag =~ /transposable_element|repeat_region/){	
			push @{$gene_refseq{$feature->seq_id}} , $feature;		
		}
		elsif ($feature->primary_tag eq "transposable_element_gene") {		 
			push @{$gene_refseq{$feature->seq_id}} , $feature;			
		}
		elsif ($feature->primary_tag eq "pseudogenic_transcript") {
			($mrna_id) = $feature->get_tag_values("ID");
			($id)  = $feature->get_tag_values("Parent"); 
			$mrna{$id}{$mrna_id} = $feature;
		}
		elsif ($feature->primary_tag eq "mRNA") {
			($mrna_id) = $feature->get_tag_values("ID");
			($id)  = $feature->get_tag_values("Parent"); 
			$mrna{$id}{$mrna_id} = $feature;
		}
		elsif ($feature->primary_tag =~ /polypeptide|protein/ ) { 
			($mrna_id) = $feature->get_tag_values("Derives_from"); 
			push @{$polypeptide{$mrna_id}}, $feature;
		}
		elsif ($feature->primary_tag eq "exon") { 	 
			($mrna_id) = $feature->get_tag_values("Parent");  
			push @{$other{$mrna_id}}, $feature;
		}
		elsif ($feature->primary_tag eq "pseudogenic_exon") {
			($mrna_id) = $feature->get_tag_values("Parent");  
			push @{$other{$mrna_id}}, $feature;
		}
		elsif ($feature->primary_tag eq "CDS") { 	 
			($mrna_id) = $feature->get_tag_values("Parent"); 
			push @{$other{$mrna_id}}, $feature;
			push @{$cds{$mrna_id}}, $feature;
		}
		elsif ($feature->primary_tag =~/.+_prime_UTR/i) { 	
			($mrna_id) = $feature->get_tag_values("Parent"); 
			push @{$other{$mrna_id}}, $feature;
		}
		elsif (($feature->primary_tag =~/.+RNA/i) && ($feature->primary_tag !~/^mRNA/i)) { 	
			($mrna_id) = $feature->get_tag_values("Parent"); 
			push @{$other_RNA{$mrna_id}}, $feature;
		}		
		elsif ($feature->primary_tag =~ /transposon_fragment|transp_element|helitron|DNA_transposon|LTR_retrotransposon|non_LTR_retrotransposon|terminal_inverted_repeat_element/) { 	
			($mrna_id) = $feature->get_tag_values("Parent"); 
			push @{$other{$mrna_id}}, $feature;
		}		
	}
	$gff->close;
	
	my $in = new Bio::SeqIO(
		-file => $genome,
		-format => "fasta"
	);
	
	my (%seq,%length); 
	while(my $seqobj = $in->next_seq){
		$seq{$seqobj->display_id} = $seqobj;
	}
	$in->close;
	
	my $file_cds = $fna_output_directory.$title2."-CDS-genfam.fna";
	my $file_prot = $fna_output_directory.$title2."-polypeptide-genfam.faa";
	my $file_exon = $fna_output_directory.$title2."-exon-genfam.fna";
	my $seq_cds;
	my $seq_pep;
	my $seq_exon;
	$seq_cds = new Bio::SeqIO(
		-file => ">$file_cds",
		-format => "fasta"
	);
	$seq_pep = new Bio::SeqIO(
		-file => ">$file_prot",
		-format => "fasta"
	);
	$seq_exon = new Bio::SeqIO(
		-file => ">$file_exon",
		-format => "fasta"
	);
	
	my $outfile = $fna_output_directory.$title."-locus_tag-genfam.gff3";
	print "Print data to $outfile...\n";
	
	my $out = new Bio::Tools::GFF(
		-file => ">$outfile",
		-gff_version => 3
	);
	
	my %h_id = ( "ARATH" => "At", "BRADI" => "Bd", "GLYMA" => "Gm", "GOSRA"   => "Gr", "LOTJA"   => "Lj", "MEDTR"   => "Mt", "MUSAC"   => "Ma",
	"ORYSI"   => "Os", "ORYSJ"   => "Os", "POPTR"   => "Pt", "SOLLC"   => "Sl", "SORBI"   => "Sb", "THECC"   => "Tc", "VITVI"   => "Vv", 
	"MAIZE"   => "Zm", "MALDO"   => "Md", "MANES"   => "Me", "RICCO"   => "Rc", "SETIT"   => "Si", "SOLTU"   => "St", "COFCA"   => "Cc");	

	my %locus;
	my %function;
	my $source_tag = "manual_curation";
	my $is_complete = 0; 
	my $is_incomplete=0;
	my $codon_start=0;
	my $codon_stop=0;
	
	foreach my $seq_id (sort {$a cmp $b} keys%gene_refseq){
		my $count = 0;
		my $seqobj;
		
		if (defined $seq{$seq_id}){
			$seqobj = $seq{$seq_id};
		}
		else{
			$seqobj = 0;
		}
		foreach my $gene (sort {$a->start <=> $b->start} @{$gene_refseq{$seq_id}}){
			if ($gene->primary_tag() eq "gene"){
				my ($name) = $gene->get_tag_values("ID");	
				$count++;		
				my $gene_id = sprintf( "_g%05d", $count * 10 );
				my $first_poly_id = sprintf( "_p%05d", $count * 10 );
				my $first_mrna_id = sprintf( "_t%05d", $count * 10 );
				my $start_gene = $gene->start;
				my $end_gene = $gene->end;
				my ($function) = $function{$name};
	
				my ($L5, $chr);
				if ($gene->seq_id =~ /^([A-Z]{5})(\d+)/i) {
					$L5 = $1;
					$chr = $2;
				} elsif($gene->seq_id =~ /^([A-Z]{5})_scaffold(\d+)/i){
					$L5 = $1;
					$chr = $2;
				} elsif($gene->seq_id =~ /^([A-Z]{5})_contig(\d+)/i){
					$L5 = $1;
					$chr = $2;
				}
	
				my $L2 = $h_id{$L5};
	
	
				my ($Gene_id) =  $gene->get_tag_values("ID");
				my ($Gene_name) =  $gene->get_tag_values("Name");
				my $locus_tag = $L2.$chr.$gene_id."_$L5";

				$locus{$Gene_id} = {   					
					locus_tag => $locus_tag,
					count => $count,
				};

				$gene->add_tag_value("locus_tag", $L2.$chr.$gene_id."_$L5");
				if ($gene->start > $gene->end){
					$gene->start($gene->end);
					$gene->end($gene->start);
				}
				$out->write_feature($gene);
		
				my $count_mrna = 0;
				if (exists $mrna{$name}){
					foreach my $mrna (keys%{$mrna{$name}}){
						$count_mrna++;
						my $mrna_id = $first_mrna_id ."." . $count_mrna;
						my $poly_id = $first_poly_id ."." . $count_mrna;
						my $feat_mrna = $mrna{$name}{$mrna};
						my ($old_mrna_id) =  $feat_mrna->get_tag_values("ID");
						my ($old_mrna_name) =  $feat_mrna->get_tag_values("Name");					
						my $start_poly = 100000000000000000;
						my $end_poly   = -1; 
						my $start_poly_exon = 100000000000000000;
						my $end_poly_exon   = -1; 
						my @cds;
						my @exon;
						my $protein_id;
						foreach my $other (sort {$a->start <=> $b->start} @{$other{$mrna}}){		
							if ($other->primary_tag() eq "CDS") {
								($protein_id) = $other->get_tag_values("protein_id") if $other->has_tag("protein_id");
								$start_poly = $other->start if $other->start < $start_poly;
								$end_poly   = $other->end if $other->end > $end_poly;
								push @cds,$other;
							}			
							if ($other->primary_tag() eq "exon") { 
								push @exon,$other;
								$start_poly_exon = $other->start if $other->start < $start_poly_exon;
								$end_poly_exon   = $other->end if $other->end > $end_poly_exon;
							}	
						}
						my $strand = $feat_mrna->strand;
	
						if ($feat_mrna->seq_id =~ /^([A-Z]{5})(\d+)/i){
							$L5 = $1;
							$chr = $2;
						}elsif($feat_mrna->seq_id =~ /^([A-Z]{5})_scaffold(\d+)/i){
							$L5 = $1;
							$chr = $2;
						}elsif($feat_mrna->seq_id =~ /^([A-Z]{5})_contig(\d+)/i){
							$L5 = $1;
							$chr = $2;
						}
			
						$feat_mrna->add_tag_value("locus_tag", $L2.$chr.$mrna_id."_$L5");
						if ($feat_mrna->start > $feat_mrna->end){
							$feat_mrna->start($feat_mrna->end);
							$feat_mrna->end($feat_mrna->start);
						}
						$out->write_feature($feat_mrna);
				
						unless (@cds ) {
							@cds = @exon;
							$start_poly = $start_poly_exon;
							$end_poly   = $end_poly_exon;
						}
						if (exists $polypeptide{$mrna}){
							foreach my $poly (@{$polypeptide{$mrna}}){
	
								if ($poly->seq_id =~ /^([A-Z]{5})(\d+)/i){
									$L5 = $1;
									$chr = $2;
								}elsif($poly->seq_id =~ /^([A-Z]{5})_scaffold(\d+)/i){
									$L5 = $1;
									$chr = $2;
								}elsif($poly->seq_id =~ /^([A-Z]{5})_contig(\d+)/i){
									$L5 = $1;
									$chr = $2;
								}
								$poly->add_tag_value("locus_tag", $L2.$chr.$poly_id."_$L5");
								if ($poly->start > $poly->end){
									$poly->start($poly->end);
									$poly->end($poly->start);
								}
								$out->write_feature($poly);
							}
						}
						else {
							my $poly = new Bio::SeqFeature::Generic(
								-seq_id 	=> $seq_id,
								-source_tag => $source_tag,
								-primary_tag => 'polypeptide',
								-start       => $start_poly,
								-end         => $end_poly,
								-strand      => $strand,
								-tag 		 => {
									ID	=> $feat_mrna->get_tag_values("ID"),
									Name	=> $feat_mrna->get_tag_values("Name"),
									Derives_from => $feat_mrna->get_tag_values("ID")
								}
							); 
							if ($feat_mrna->has_tag("Dbxref")) {
								my @dbxref = $feat_mrna->get_tag_values("Dbxref");
								foreach my $dbxref (@dbxref) {
									$poly->add_tag_value("Dbxref",$dbxref);
								}
							}
							if ($poly->seq_id =~ /^([A-Z]{5})(\d+)/i){
								$L5 = $1;
								$chr = $2;
							}elsif($poly->seq_id =~ /^([A-Z]{5})_scaffold(\d+)/i){
								$L5 = $1;
								$chr = $2;
							}elsif($poly->seq_id =~ /^([A-Z]{5})_contig(\d+)/i){
								$L5 = $1;
								$chr = $2;
							}
							$poly->add_tag_value("locus_tag", $L2.$chr.$poly_id."_$L5");
							if ($poly->start > $poly->end){
								$poly->start($poly->end);
								$poly->end($poly->start);
							}
							$out->write_feature($poly);
						}
			
						foreach my $other (sort {$a->start <=> $b->start} @{$other{$old_mrna_id}}){
							$out->write_feature($other);						
						}
						
						if ($seqobj != 0){
							foreach my $feature (@cds){
								if ($feature->start > $feature->end){
									$feature->start($feature->end);
									$feature->end($feature->start);
								}
							}
							if ($cds[0]->strand == 1) {
								@cds = sort{$a->start <=> $b->start} @cds;
							}
							else {
								@cds = sort{$b->start <=> $a->start} @cds;
							}
							my $cds;
							foreach my $feature (@cds) {
								if ($feature->end < $seqobj->length()){
									my $seqtrunc_cds = $seqobj->trunc($feature->start,$feature->end);	
									if ($feature->strand == -1) {
										$cds .= $seqtrunc_cds->revcom()->seq;
									}
									else {	
										$cds .= $seqtrunc_cds->seq;
									}
								}
							}
							
							my $seqobj_cds = Bio::PrimarySeq->new (
								-display_id  => $old_mrna_name."_$title3",
								-seq         => $cds,
								-desc		=> $function
							);  
							
							foreach my $feature (@exon){
								if ($feature->start > $feature->end){
									$feature->start($feature->end);
									$feature->end($feature->start);
								}
							}
							if (defined $exon[0]){
								if ($exon[0]->strand == 1) {
									@exon = sort{$a->start <=> $b->start} @exon;
								}
								else {
									@exon = sort{$b->start <=> $a->start} @exon;
								}
								my $exon;
								foreach my $feature (@exon) {
									if ($feature->end < $seqobj->length()){									
										my $seqtrunc_exon = $seqobj->trunc($feature->start,$feature->end);	
										if ($feature->strand == -1) {
											$exon .= $seqtrunc_exon->revcom()->seq;
										}
										else {	
											$exon .= $seqtrunc_exon->seq;
										}
									}
								}
								my $seqobj_exon = Bio::PrimarySeq->new (
									-display_id  => $old_mrna_name."_$title3",
									-seq         => $exon,
									-desc		=> $function
								); 
								$seq_exon->write_seq($seqobj_exon);
							}  
							
							if ($cds =~ /^ATG.*/ &&  $cds =~ /.*(TAG|TAA|TGA)$/){
								$is_complete++;
							}
							else {
								if ($cds =~ /^ATG.*/){
									$codon_stop++;
								}
								elsif ($cds =~ /.*(TAG|TAA|TGA)$/){
									$codon_start++;
								}
								else {
									$is_incomplete++;
									#print $name ,"\t", $old_mrna_id,"\n";
								}
							}
							
							my $seqobj_pep = $seqobj_cds->translate();
							$seqobj_pep->display_id($old_mrna_name."_$title3");
							$seq_cds->write_seq($seqobj_cds);
							$seq_pep->write_seq($seqobj_pep);
						}
					}
				}
				elsif(exists $other_RNA{$name}){
					foreach my $other (sort {$a->start <=> $b->start} @{$other_RNA{$name}}){
						my ($id) = $other->get_tag_values("ID");
						$out->write_feature($other);
						if (defined $other{$id}){
							foreach my $other (sort {$a->start <=> $b->start} @{$other{$id}}){
								$out->write_feature($other);
							}		
						}
					}
				}
				elsif(not exists $mrna{$name}) {
					$count_mrna++;
					my $strand = $gene->strand;
					my $mrna = new Bio::SeqFeature::Generic(
						-seq_id 	=> $seq_id,
						-source_tag => $source_tag,
						-primary_tag => 'mRNA',
						-start       => $start_gene,
						-end         => $end_gene,
						-strand      => $strand,
						-tag 		 => {
							ID	    => "$Gene_id\.$count_mrna",
							Name	=> "$Gene_name\.$count_mrna",
							Parent  => $Gene_id
						}
					);							
			
					my $mrna_id = $first_mrna_id ."." . $count_mrna;
					my $poly_id = $first_poly_id ."." . $count_mrna;
					my $feat_mrna = $mrna;
					my ($old_mrna_id) =  $feat_mrna->get_tag_values("Parent");
					my ($old_mrna_name) =  $feat_mrna->get_tag_values("Name");
					my $start_poly = 100000000000000000;
					my $end_poly   = -1; 
					my $start_poly_exon = 100000000000000000;
					my $end_poly_exon   = -1; 
					my @cds;
					my @exon;
					my $protein_id;
					foreach my $other (sort {$a->start <=> $b->start} @{$other{$old_mrna_id}}){		
						if ($other->primary_tag() eq "CDS") {
							($protein_id) = $other->get_tag_values("protein_id") if $other->has_tag("protein_id");
							$start_poly = $other->start if $other->start < $start_poly;
							$end_poly   = $other->end if $other->end > $end_poly;
							push @cds,$other;
						}			
						if ($other->primary_tag() eq "exon") { 
							push @exon,$other;
							$start_poly_exon = $other->start if $other->start < $start_poly_exon;
							$end_poly_exon   = $other->end if $other->end > $end_poly_exon;
						}	
					}
		
					if ($feat_mrna->seq_id =~ /^([A-Z]{5})(\d+)/i){
						$L5 = $1;
						$chr = $2;
					}elsif($feat_mrna->seq_id =~ /^([A-Z]{5})_scaffold(\d+)/i){
						$L5 = $1;
						$chr = $2;
					}elsif($feat_mrna->seq_id =~ /^([A-Z]{5})_contig(\d+)/i){
						$L5 = $1;
						$chr = $2;
					}

					$feat_mrna->add_tag_value("locus_tag", $L2.$chr.$mrna_id."_$L5");
					if ($feat_mrna->start > $feat_mrna->end){
						$feat_mrna->start($feat_mrna->end);
						$feat_mrna->end($feat_mrna->start);
					}
					$out->write_feature($feat_mrna);
					unless (@cds ) {
						@cds = @exon;
						$start_poly = $start_poly_exon;
						$end_poly   = $end_poly_exon;
					}
					if (exists $polypeptide{$old_mrna_id}){
						foreach my $poly (@{$polypeptide{$old_mrna_id}}){
											
							if ($poly->seq_id =~ /^([A-Z]{5})(\d+)/i){
								$L5 = $1;
								$chr = $2;
							}elsif($poly->seq_id =~ /^([A-Z]{5})_scaffold(\d+)/i){
								$L5 = $1;
								$chr = $2;
							}elsif($poly->seq_id =~ /^([A-Z]{5})_contig(\d+)/i){
								$L5 = $1;
								$chr = $2;
							}
							$poly->add_tag_value("locus_tag", $L2.$chr.$poly_id."_$L5");
							if ($poly->start > $poly->end){
								$poly->start($poly->end);
								$poly->end($poly->start);
							}
							$out->write_feature($poly);
						}
					}
					else {
						my $poly = new Bio::SeqFeature::Generic(
							-seq_id 	=> $seq_id,
							-source_tag => $source_tag,
							-primary_tag => 'polypeptide',
							-start       => $start_poly,
							-end         => $end_poly,
							-strand      => $strand,
							-tag 		 => {
								ID	=> "$Gene_id\.$count_mrna",
								Name	=> "$Gene_name\.$count_mrna",
								Derives_from => "$Gene_id\.$count_mrna",
							}
						); 
						if ($feat_mrna->has_tag("Dbxref")) {
							my @dbxref = $feat_mrna->get_tag_values("Dbxref");
							foreach my $dbxref (@dbxref) {
								$poly->add_tag_value("Dbxref",$dbxref);
							}
						}
						if ($poly->seq_id =~ /^([A-Z]{5})(\d+)/i){
							$L5 = $1;
							$chr = $2;
						}elsif($poly->seq_id =~ /^([A-Z]{5})_scaffold(\d+)/i){
							$L5 = $1;
							$chr = $2;
						}elsif($poly->seq_id =~ /^([A-Z]{5})_contig(\d+)/i){
							$L5 = $1;
							$chr = $2;
						}
						$poly->add_tag_value("locus_tag", $L2.$chr.$poly_id."_$L5");
						if ($poly->start > $poly->end){
							$poly->start($poly->end);
							$poly->end($poly->start);
						}
						$out->write_feature($poly);
					}
					foreach my $other (sort {$a->start <=> $b->start} @{$other{$old_mrna_id}}){
						$out->write_feature($other);						
					}
					
					if ($seqobj != 0){
						foreach my $feature (@cds){
							if ($feature->start > $feature->end){
								$feature->start($feature->end);
								$feature->end($feature->start);
							}
						}
						if ($cds[0]->strand == 1) {
							@cds = sort{$a->start <=> $b->start} @cds;
						}
						else {
							@cds = sort{$b->start <=> $a->start} @cds;
						}
						my $cds;
						foreach my $feature (@cds) {
							my $seqtrunc_cds = $seqobj->trunc($feature->start,$feature->end);	
							if ($feature->strand == -1) {
								$cds .= $seqtrunc_cds->revcom()->seq;
							}
							else {	
								$cds .= $seqtrunc_cds->seq;
							}
						}
						my $seqobj_cds = Bio::PrimarySeq->new (
							-display_id  => $old_mrna_name."_$title3",
							-seq         => $cds,
							-desc		=> $function
						);  
						
						foreach my $feature (@exon){
							if ($feature->start > $feature->end){
								$feature->start($feature->end);
								$feature->end($feature->start);
							}
						}
						if ($exon[0]->strand == 1) {
							@exon = sort{$a->start <=> $b->start} @exon;
						}
						else {
							@exon = sort{$b->start <=> $a->start} @exon;
						}
						my $exon;
						foreach my $feature (@exon) {	
							my $seqtrunc_exon = $seqobj->trunc($feature->start,$feature->end);	
							if ($feature->strand == -1) {
								$exon .= $seqtrunc_exon->revcom()->seq;
							}
							else {	
								$exon .= $seqtrunc_exon->seq;
							}
						}
						my $seqobj_exon = Bio::PrimarySeq->new (
							-display_id  => $old_mrna_name."_$title3",
							-seq         => $exon,
							-desc		=> $function
						);   
						if ($cds =~ /^ATG.*/ &&  $cds =~ /.*(TAG|TAA|TGA)$/){
							$is_complete++;
						}
						else {
							if ($cds =~ /^ATG.*/){
								$codon_stop++;
							}
							elsif ($cds =~ /.*(TAG|TAA|TGA)$/){
								$codon_start++;
							}
							else {
								$is_incomplete++;
								#print $name ,"\t", $old_mrna_id,"\n";
							}
						}
						my $seqobj_pep = $seqobj_cds->translate();
						$seqobj_pep->display_id($old_mrna_name."_$title3");
						$seq_cds->write_seq($seqobj_cds);
						$seq_pep->write_seq($seqobj_pep);
						$seq_exon->write_seq($seqobj_exon);	
					}		
				}											
			}
			elsif  ($gene->primary_tag() eq "pseudogene"){
				my ($name) = $gene->get_tag_values("ID");
				$out->write_feature($gene);
				foreach my $mrna (keys%{$mrna{$name}}){
					my $feat_mrna = $mrna{$name}{$mrna};
					my ($old_mrna_id) =  $feat_mrna->get_tag_values("ID");	
					$out->write_feature($feat_mrna);						
					foreach my $other (@{$other{$old_mrna_id}}){
						$out->write_feature($other);						
					}				
				}
			}
			elsif ($gene->primary_tag() eq "transposable_element_gene"){
				my ($name) = $gene->get_tag_values("ID");
				$out->write_feature($gene);
				foreach my $mrna (keys%{$mrna{$name}}){
					my $feat_mrna = $mrna{$name}{$mrna};
					my ($old_mrna_id) =  $feat_mrna->get_tag_values("ID");	
					$out->write_feature($feat_mrna);						
					foreach my $other (@{$other{$old_mrna_id}}){
						$out->write_feature($other);						
					}				
				}
			}
			elsif (($gene->primary_tag =~ /transposable_element|repeat_region/)){
				my ($name) = $gene->get_tag_values("ID");
				$out->write_feature($gene);						
				foreach my $other (@{$other{$name}}){
					$out->write_feature($other);						
				}				
			}
		}			
	}
	my $json = encode_json \%locus;
	my $file = $fna_output_directory.$title2."-locus_tag.json";
	unless(open FILE, '>'.$file) {
		die "Unable to create $file";
	}
	print FILE $json;
	close FILE;
	$out->close();
}

sub gff3tobed{

	my $new_file = $fna_output_directory.$title."-locus_tag-genfam.gff3";
	my $gffio = Bio::Tools::GFF -> new(-file =>$new_file , -gff_version => 3);

	my $strand=0;
	my $strand1='-';
	
	# ##### Liste les fichiers ayant l'extension .bed 
# 	my @list = glob("*.bed"); 
#   
# 	### recupere le nombre de fichier 
# 	my $numberoffile = scalar(@list);
# 	
# 	for ( my $v = 0; $v < $numberoffile; $v++ ) {
# 		unlink $list[$v]; #supprime les fichiers bed
# 	}
	
	my $feat_start;
	my $feat_end;
	my ($FILE1, $FILE2, $FILE3, $FILE4, $FILE5, $FILE6, $FILE7, $FILE8);
	
	if($bed =~/gene|all/i){
		my $file1 = $fna_output_directory.$title."-gene-genfam.bed";
		open($FILE1,">".$file1);
	}
	
	if($bed =~/mRNA|all/i){
		my $file2 = $fna_output_directory.$title."-mRNA-genfam.bed";
		open($FILE2,">".$file2);
	}
	
	if($bed =~/five_prime_UTR|all/i){
		my $file3 = $fna_output_directory.$title."-five_prime_UTR-genfam.bed";
		open($FILE3,">".$file3);
	}
	
	if($bed =~/three_prime_UTR|all/i){
		my $file4 = $fna_output_directory.$title."-three_prime_UTR-genfam.bed";
		open($FILE4,">".$file4);
	}
	
	if($bed =~/CDS|all/i){
		my $file5 = $fna_output_directory.$title."-CDS-genfam.bed";
		open($FILE5,">".$file5);
	}
	
	if($bed =~/exon|all/i){
		my $file6 = $fna_output_directory.$title."-exon-genfam.bed";
		open($FILE6,">".$file6);
	}
	
	if($bed =~/intron|all/i){
		my $file7 = $fna_output_directory.$title."-intron-genfam.bed";
		open($FILE7,">".$file7);
	}
	
	if($bed =~/polypeptide|all/i){
		my $file8 = $fna_output_directory.$title."-polypeptide-genfam.bed";
		open($FILE8,">".$file8);
	}
	
	while(my $feature = $gffio->next_feature()) {
		if(($feature->primary_tag =~/gene/i) && ($bed =~/gene|all/i)){
			
			# my $file = $fna_output_directory.$title."-gene-genfam.bed";
# 			unless(open FILE, '>>'.$file) {
# 				die "Unable to create $file";
# 			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
			}else{
				$strand1='+';
			}
			if ($feature->start > $feature->end){
				$feat_start = $feature->end;
				$feat_end = $feature->start;
			}
			else{
				$feat_start = $feature->start;
				$feat_end = $feature->end;
			}
			
			#gff2bed:
			print $FILE1 join("\t",$feature->seq_id,$feat_start-1,$feat_end,$feature->{_gsf_tag_hash}->{Name}->[0]."_".$title3,".",$strand1)."\n";
			#close FILE;
		}
	
		if(($feature->primary_tag =~/mRNA/i) && ($bed =~/mRNA|all/i)){

			# my $file = $fna_output_directory.$title."-mRNA-genfam.bed";
# 			unless(open FILE, '>>'.$file) {
# 				die "Unable to create $file";
# 			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
			}else{
				$strand1='+';
			}
			if ($feature->start > $feature->end){
				$feat_start = $feature->end;
				$feat_end = $feature->start;
			}
			else{
				$feat_start = $feature->start;
				$feat_end = $feature->end;
			}
			
			#gff2bed:
			print $FILE2 join("\t",$feature->seq_id,$feat_start-1,$feat_end,$feature->{_gsf_tag_hash}->{Name}->[0]."_".$title3,".",$strand1)."\n";
			#close FILE;
		}
	
		if(($feature->primary_tag =~/five_prime_UTR/i) && ($bed =~/five_prime_UTR|all/i)){

			# my $file = $fna_output_directory.$title."-five_prime_UTR-genfam.bed";
# 			unless(open FILE, '>>'.$file) {
# 				die "Unable to create $file";
# 			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
			}else{
				$strand1='+';
			}
			if ($feature->start > $feature->end){
				$feat_start = $feature->end;
				$feat_end = $feature->start;
			}
			else{
				$feat_start = $feature->start;
				$feat_end = $feature->end;
			}
			
			#gff2bed:
			print $FILE3 join("\t",$feature->seq_id,$feat_start-1,$feat_end,$feature->{_gsf_tag_hash}->{Parent}->[0]."_".$title3,".",$strand1)."\n";
			#close FILE;
		}
	
		if(($feature->primary_tag =~/three_prime_UTR/i) && ($bed =~/three_prime_UTR|all/i)){

			# my $file = $fna_output_directory.$title."-three_prime_UTR-genfam.bed";
# 			unless(open FILE, '>>'.$file) {
# 				die "Unable to create $file";
# 			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
			}else{
				$strand1='+';
			}
			if ($feature->start > $feature->end){
				$feat_start = $feature->end;
				$feat_end = $feature->start;
			}
			else{
				$feat_start = $feature->start;
				$feat_end = $feature->end;
			}
			
			#gff2bed:
			print $FILE4 join("\t",$feature->seq_id,$feat_start-1,$feat_end,$feature->{_gsf_tag_hash}->{Parent}->[0]."_".$title3,".",$strand1)."\n";
			#close FILE;
		}
		
		if(($feature->primary_tag =~/CDS/i) && ($bed =~/CDS|all/i)){

			# my $file = $fna_output_directory.$title."-CDS-genfam.bed";
# 			unless(open FILE, '>>'.$file) {
# 				die "Unable to create $file";
# 			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
			}else{
				$strand1='+';
			}
			if ($feature->start > $feature->end){
				$feat_start = $feature->end;
				$feat_end = $feature->start;
			}
			else{
				$feat_start = $feature->start;
				$feat_end = $feature->end;
			}
			
			#gff2bed:
			print $FILE5 join("\t",$feature->seq_id,$feat_start-1,$feat_end,$feature->{_gsf_tag_hash}->{Parent}->[0]."_".$title3,".",$strand1)."\n";
			#close FILE;
		}
	
		if(($feature->primary_tag =~/exon/i) && ($bed =~/exon|all/i)){

			# my $file = $fna_output_directory.$title."-exon-genfam.bed";
# 			unless(open FILE, '>>'.$file) {
# 				die "Unable to create $file";
# 			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
			}else{
				$strand1='+';
			}
			if ($feature->start > $feature->end){
				$feat_start = $feature->end;
				$feat_end = $feature->start;
			}
			else{
				$feat_start = $feature->start;
				$feat_end = $feature->end;
			}
			
			#gff2bed:
			print $FILE6 join("\t",$feature->seq_id,$feat_start-1,$feat_end,$feature->{_gsf_tag_hash}->{Parent}->[0]."_".$title3,".",$strand1)."\n";
			#close FILE;
		}
	
		if(($feature->primary_tag =~/intron/i) && ($bed =~/intron|all/i)){
# 
# 			my $file = $fna_output_directory.$title."-intron-genfam.bed";
# 			unless(open FILE, '>>'.$file) {
# 				die "Unable to create $file";
# 			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
			}else{
				$strand1='+';
			}
			if ($feature->start > $feature->end){
				$feat_start = $feature->end;
				$feat_end = $feature->start;
			}
			else{
				$feat_start = $feature->start;
				$feat_end = $feature->end;
			}
			
			#gff2bed:
			print $FILE7 join("\t",$feature->seq_id,$feat_start-1,$feat_end,$feature->{_gsf_tag_hash}->{Parent}->[0]."_".$title3,".",$strand1)."\n";
			#close FILE;
		}
	
		if(($feature->primary_tag =~/polypeptide|protein/i) && ($bed =~/polypeptide|protein|all/i)){

			# my $file = $fna_output_directory.$title."-polypeptide-genfam.bed";
# 			unless(open FILE, '>>'.$file) {
# 				die "Unable to create $file";
# 			}
			#change -1 to - and 1 to +
			if($feature->strand=~/-/){
				$strand1='-';
			}else{
				$strand1='+';
			}
			if ($feature->start > $feature->end){
				$feat_start = $feature->end;
				$feat_end = $feature->start;
			}
			else{
				$feat_start = $feature->start;
				$feat_end = $feature->end;
			}
			
			#gff2bed:
			print $FILE8 join("\t",$feature->seq_id,$feat_start-1,$feat_end,$feature->{_gsf_tag_hash}->{Name}->[0]."_".$title3,".",$strand1)."\n";
			#close FILE;
		}
	}
	$gffio->close();
	close $FILE1;
	close $FILE2;
	close $FILE3;
	close $FILE4;
	close $FILE5;
	close $FILE6;
	close $FILE7;
	close $FILE8;
  
}

sub bedtools {

	my @tab_type = ("gene", "mRNA", "intron", "five_prime_UTR", "three_prime_UTR");
	
	#bedtools
	#usage system: system("bedtools getfasta [OPTIONS] -fi <input FASTA> -bed <BED/GFF/VCF> -fo <output FASTA>");
	foreach my $VAR (@tab_type){
		my $fich = $bed_output_directory.$title."-".$VAR."-genfam.bed";
		if (-e $fich){
			#print "bedtools getfasta -s -fi $genome -bed ".$fich." -fo ".$fna_output_directory.$title2."-".$VAR."-genfam.fna -name";
			system("bedtools getfasta -s -fi $genome -bed ".$fich." -fo ".$fna_output_directory.$title2."-".$VAR."-genfam.fna -name");
		}
	}	
}

sub gff3_stat {
	my $file = $fna_output_directory.$title."-locus_tag-genfam.gff3";
	
	my $gff1 = new Bio::Tools::GFF(
		-file => $file,
		-gff_version => 3
	);
	
	my $gff = new Bio::Tools::GFF(
		-file => $file,
		-gff_version => 3
	);
	my (%stats, %gene, %te_gene, %te, %pseudogene, %mrna, %cds, %three, %five, %exon, %chr_length, $seq_obj);
	
	while(my $feature = $gff1->next_feature ){
		$chr_length{$feature->seq_id} = $feature->end;
	}
	$gff1->close; 
	
	my @ftypes = qw(gene transposable_element_gene mrna te pseudogene cds exon five_prime_utr three_prime_utr intron);
	for my $t ( @ftypes ) {
		$stats{$t} = Statistics::Descriptive::Full->new;
	}
	my $id;
	my $mrna_id;	
	while(my $feature = $gff->next_feature) {
		if ($feature->primary_tag() eq "gene") {
			($id) = $feature->get_tag_values("ID");
			push @{$gene{$feature->seq_id}{$id}},$feature;
		}
		if ($feature->primary_tag() eq "mRNA") {
			($id) = $feature->get_tag_values("Parent");
			($mrna_id) = $feature->get_tag_values("ID");
			push @{$mrna{$id}{$mrna_id}},$feature;
		}		
		if ($feature->primary_tag() eq "transposable_element_gene") {
			($id) = $feature->get_tag_values("ID");
			push @{$te_gene{$feature->seq_id}{$id}},$feature;
		}
		if ($feature->primary_tag =~ /transposable_element|repeat_region/) {
			($id) = $feature->get_tag_values("ID");	
			push @{$te{$feature->seq_id}{$id}},$feature;
		}	
		if ($feature->primary_tag() eq "CDS") {
			my ($parent) = $feature->get_tag_values("Parent");
			push @{$cds{$parent}},$feature;
		}	
		if ($feature->primary_tag() eq "three_prime_UTR") {
			my ($parent) = $feature->get_tag_values("Parent");
			push @{$three{$parent}},$feature;
		}	
		if ($feature->primary_tag() eq "five_prime_UTR") {
			my ($parent) = $feature->get_tag_values("Parent");
			push @{$five{$parent}},$feature;
		}	
		if ($feature->primary_tag() eq "exon") {
			my ($parent) = $feature->get_tag_values("Parent");
			push @{$exon{$parent}},$feature;
		}	
		if ($feature->primary_tag eq "pseudogene") {
			($id) = $feature->get_tag_values("ID");				
			push @{$pseudogene{$feature->seq_id}{$id}},$feature;
		}
	}
	$gff->close;
	open(OUT,">".$fna_output_directory."stat.txt");
	my $genome = 0;
	my $gene_sum = 0;
	my $te_sum  = 0;
	my $te_gene_sum = 0;
	my $pseudo_sum = 0;

	print OUT "TABLE:\n";
	print OUT join("\t","Chr","Length (bp)","Num Gene", "Num TE_gene", "Num pseudogene", "Num TE"),"\n";
	foreach my $seq_id (sort {$a cmp $b} keys %chr_length){
		my $length 	= $chr_length{$seq_id};
		my @gene;
		if (defined $gene{$seq_id}){
			@gene = keys%{$gene{$seq_id}};
		}
		
		my @gene_te;
		if (defined $te_gene{$seq_id}){
			@gene_te = keys%{$te_gene{$seq_id}};
		}
		
		my @te;
		if (defined $te{$seq_id}){
			@te = keys%{$te{$seq_id}};
		}
		
		my @pseudo;
		if (defined $pseudogene{$seq_id}){
			@pseudo = keys%{$pseudogene{$seq_id}};
		}
		
		$genome	  += $length;
		$gene_sum += scalar(@gene);
		$te_sum   += scalar(@te);
		$te_gene_sum += scalar(@gene_te);
		$pseudo_sum += scalar(@pseudo);
		
		push @gene,@gene_te,@te,@pseudo; 
		print OUT join("\t",$seq_id,$chr_length{$seq_id},scalar(@gene),scalar(@gene_te),scalar(@pseudo),scalar(@te)),"\n";
		
		my $total = 0; 
		#my $lastfeature;
		#my $cpt = 0;
		#my $genic;
		#my $intergenic;
		foreach my $gene_id (keys%{$gene{$seq_id}}){			
			if (defined $gene{$seq_id}{$gene_id}){
				foreach my $feature_gene (sort{$a->start <=> $b->start} @{$gene{$seq_id}{$gene_id}}) {
					my $gene_length = abs($feature_gene->start - $feature_gene->end);
					$stats{'gene'}->add_data($gene_length); 
					$total++;
				}
			}	
					
			foreach my $mrna_id (keys%{$mrna{$gene_id}}){
				#$cpt++;				
				if (defined $cds{$mrna_id}){
					foreach my $feature_cds (sort{$a->start <=> $b->start} @{$cds{$mrna_id}}) {
						my $cds_length = abs($feature_cds->start - $feature_cds->end);
						$stats{'cds'}->add_data($cds_length); 
						$total++;
					}
				}
				
				if (defined $five{$mrna_id}){
					foreach my $feature_five (sort{$a->start <=> $b->start} @{$five{$mrna_id}}) {
						my $length = abs($feature_five->start - $feature_five->end);
						$stats{'five_prime_utr'}->add_data($length); 
						$total++;
					}
				}
				
				if (defined $three{$mrna_id}){
					foreach my $feature_three (sort{$a->start <=> $b->start} @{$three{$mrna_id}}) {
						my $length = abs($feature_three->start - $feature_three->end);
						$stats{'three_prime_utr'}->add_data($length); 
						$total++;
					}
				}
				my $lastexon;
		
				if (defined $exon{$mrna_id}){
					for my $exon ( sort { $a->start  <=>   $b->start  }  @{$exon{$mrna_id}} ) {
						my $exonlen = abs($exon->start - $exon->end);
						$stats{'exon'}->add_data($exonlen);
						if( $lastexon ) {
							my $intronlen = abs($exon->start - ($lastexon->end)); 
							$stats{'intron'}->add_data($intronlen);
						}
						$lastexon = $exon;
					}
				}
			}
			# if ($lastfeature) {
# 				$intergenic = $feature->start - $lastfeature->end;
# 				if ($intergenic > 0) {
# 					$stats{'intergenic'}->add_data($intergenic);  
# 					#print OUT join("\t",$cpt,"Inter",$seq_id,$lastfeature->end,$feature->start,$intergenic),"\n";  
# 				} 
# 				$genic = abs($feature->start - $feature->end);
# 				if ($feature->primary_tag() eq "mRNA") {
# 					$stats{'gene'}->add_data($genic); 
# 				}
# 				else {
# 					$stats{'transposable_element_gene'}->add_data($genic);
# 					$stats{'te'}->add_data($genic); 
# 					$stats{'pseudogene'}->add_data($genic);
# 				}
# 				$stats{'genic'}->add_data($genic);  
# 				#print OUT join("\t",$cpt,$feature->primary_tag,$seq_id,$feature->start,$feature->end,$genic),"\n"; 
# 				$lastfeature = $feature;
# 			}
# 			else {
# 				$intergenic = $feature->start;
# 				$stats{'intergenic'}->add_data($intergenic);   
# 				$genic = abs($feature->start - $feature->end);
# 				print OUT join("\t",$cpt,$feature->primary_tag,$seq_id,$feature->start,$feature->end,$genic),"\n";
# 				$stats{'genic'}->add_data($genic);  
# 				if ($feature->primary_tag() eq "mRNA") {
# 					$stats{'gene'}->add_data($genic); 
# 				}
# 				else {
# 					$stats{'transposable_element_gene'}->add_data($genic);
# 					$stats{'te'}->add_data($genic); 
# 					$stats{'pseudogene'}->add_data($genic);
# 				}
# 				$lastfeature = $feature;
# 			}
		}
		foreach my $te_gene_id (keys%{$te_gene{$seq_id}}){
			if (defined $te_gene{$seq_id}{$te_gene_id}){
				foreach my $feature_te_gene (sort{$a->start <=> $b->start} @{$te_gene{$seq_id}{$te_gene_id}}) {
					my $te_gene_length = abs($feature_te_gene->start - $feature_te_gene->end);
					$stats{'transposable_element_gene'}->add_data($te_gene_length); 
					$total++;
				}
			}		
		}
		foreach my $te_id (keys%{$te{$seq_id}}){
			if (defined $te{$seq_id}{$te_id}){
				foreach my $feature_te (sort{$a->start <=> $b->start} @{$te{$seq_id}{$te_id}}) {
					my $te_length = abs($feature_te->start - $feature_te->end);
					$stats{'te'}->add_data($te_length); 
					$total++;
				}
			}		
		}
		foreach my $pseudo_id (keys%{$pseudogene{$seq_id}}){
			if (defined $pseudogene{$seq_id}{$pseudo_id}){
				foreach my $feature_pseudo (sort{$a->start <=> $b->start} @{$pseudogene{$seq_id}{$pseudo_id}}) {
					my $pseudo_length = abs($feature_pseudo->start - $feature_pseudo->end);
					$stats{'pseudogene'}->add_data($pseudo_length); 
					$total++;
				}
			}		
		}
		#$intergenic = abs($length - $gene[-1]->end ) ;
		#$stats{'intergenic'}->add_data($intergenic);  
		#print OUT join("\t",$cpt,"Intergenic",$seq_id,$gene[-1]->end +1,$length, $mrna_id),"\n\n";
	}
	print OUT join("\t","All",$genome,$gene_sum,$te_gene_sum,$pseudo_sum,$te_sum),"\n";
	for my $t ( qw (gene transposable_element_gene te pseudogene cds exon five_prime_utr three_prime_utr intron) ) {
		my $percent = 100 * $stats{$t}->sum / $genome;
		print OUT join("\t",$t,$stats{$t}->count,$stats{$t}->sum,$stats{$t}->mean,$percent),"\n";
   
	}
}

sub flanking_region(){

	my $new_file = $gff3_file;
	my $gff = new Bio::Tools::GFF(	
		-file => $new_file,
		-gff_version => 3
	);
	
	my $in = new Bio::SeqIO(
		-file => $genome,
		-format => "fasta"
	);
	
	my (%length, $seq_obj);
	while ($seq_obj = $in->next_seq){  
  		 $length{$seq_obj->primary_id} = $seq_obj->length();
	}
	$in->close;
	
	my $cpt = 0;
	
	my (%gene, %mRNA, %CDS, %exon, %five_prime_UTR, $mrna_id, $id, %transposable_region);
	
	while(my $feature = $gff->next_feature) {
		if ($feature->primary_tag() eq "gene") {
			$cpt++;			
			($id) = $feature->get_tag_values("ID");
			$gene{$feature->seq_id}{$cpt} = $feature;			
		}
		elsif ($feature->primary_tag eq "mRNA") {
			($mrna_id) = $feature->get_tag_values("ID");
			($id)  = $feature->get_tag_values("Parent");
			$mRNA{$id}{$mrna_id} = $feature;
		}
		elsif ($feature->primary_tag eq "exon") { 	 
			($mrna_id) = $feature->get_tag_values("Parent");  
			push @{$exon{$mrna_id}}, $feature;
		}
		elsif ($feature->primary_tag eq "CDS") {
			($mrna_id) = $feature->get_tag_values("Parent"); 
			push @{$CDS{$mrna_id}}, $feature;
		}
		elsif ($feature->primary_tag =~/five_prime_UTR/i) { 	
			($mrna_id) = $feature->get_tag_values("Parent"); 
			push @{$five_prime_UTR{$mrna_id}}, $feature;
		}
		#transposon_fragment|DNA_transposon|LTR_retrotransposon|tRNA|rRNA
		elsif ($feature->primary_tag =~/transposable_element|repeat_region/i) {
			$cpt++; 	
			($mrna_id) = $feature->get_tag_values("ID"); 
			$gene{$feature->seq_id}{$cpt} = $feature;
		}
	}
	$gff->close;

	my $outfile = $bed_output_file;
	print "Print data to $outfile...\n";
	
	open(FILE, ">$outfile") or die "Could not open file '$outfile' $!";
	
	foreach my $seq_id (sort {$a cmp $b} keys%gene){
		foreach my $gene (sort {$a <=> $b} keys%{$gene{$seq_id}}){
			
			my ($feature) = $gene{$seq_id}{$gene};
			my $score = ".";
			my $strand;
			
			if($feature->strand=~/-/){
				$strand='-';
			}else{
				$strand='+';
			}
			
			if (defined($length{$seq_id})){
				#Si region 5’ et brin + ou région 3’ et brin - et non transposable_element_gene
				if (($feature->strand=~/^1/) && ($flank =~ /prom/) && ($feature->primary_tag !~/transposable_element_gene/) 
				|| ($feature->strand=~/^-1/) && ($flank =~ /term/) && ($feature->primary_tag !~/transposable_element_gene/)){
					my ($feature_cds);
					my ($size);	
					my ($mrna_id);
					my ($id) = $gene{$seq_id}{$gene}->get_tag_values("ID");
			
					my @keys = keys%{$mRNA{$id}};	
					$mrna_id = $keys[0];
					if (defined $mrna_id){
						if (defined $CDS{$mrna_id}[0]){
							$size = scalar(@{$CDS{$mrna_id}});
							if ($feature->strand =~/^1/){
								$feature_cds = $CDS{$mrna_id}[0];
							}
							else{
								$feature_cds = $CDS{$mrna_id}[$size-1]
							}						
						}					
					}
					if (0+(keys%{$gene{$seq_id}}) > 1){
						#Cas du premier gène du chromosome
						if (not exists $gene{$seq_id}{$gene-1}){
							if (defined $feature_cds){							
								if($feature_cds->start+$begin > 0){						
									print FILE join("\t",$feature_cds->seq_id,$feature_cds->start-1+$begin,$feature_cds->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";								
								}
								elsif($feature_cds->start+$begin <= 0){
									print FILE join("\t",$feature_cds->seq_id,$feature_cds->start-$feature_cds->start,$feature_cds->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
								}
							}
						}				
						else{
							if (defined $feature_cds){					
								#Cas du gène dans un autre gène (intron)
								if ($gene{$seq_id}{$gene-1}->start < $feature_cds->start){							
									if(($gene{$seq_id}{$gene-1}->end >= $feature_cds->start)){
										my ($id) = $gene{$seq_id}{$gene-1}->get_tag_values("ID");														
										my (@keys) = keys%{$mRNA{$id}};	
										my ($mrna_id) = $keys[0];								
										if ((defined $mrna_id)){
											if (defined $CDS{$mrna_id}[0]){
												my ($cpt)=0;
												my $size = scalar(@{$CDS{$mrna_id}});								
												while (($feature_cds->start	> $CDS{$mrna_id}[$cpt]->end) && ($cpt < $size-1)){
													$cpt += 1;
												}
												if ($CDS{$mrna_id}[0]->end > $feature_cds->start){
													if ($gene{$seq_id}{$gene-1}->start >= $feature_cds->start+$begin){
														print FILE join("\t",$feature_cds->seq_id,$gene{$seq_id}{$gene-1}->start,$feature_cds->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
													}		
													elsif($gene{$seq_id}{$gene-1}->start < $feature_cds->start+$begin){
														print FILE join("\t",$feature_cds->seq_id,$feature_cds->start-1+$begin,$feature_cds->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
													}
												}
												else{
													if ($CDS{$mrna_id}[$cpt-1]->end >= $feature_cds->start+$begin){
														print FILE join("\t",$feature_cds->seq_id,$CDS{$mrna_id}[$cpt-1]->end,$feature_cds->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
													}		
													elsif($CDS{$mrna_id}[$cpt-1]->end < $feature_cds->start+$begin){
														print FILE join("\t",$feature_cds->seq_id,$feature_cds->start-1+$begin,$feature_cds->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
													}
												}	
											}
											else{
												print FILE join("\t",$feature_cds->seq_id,1,1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
											}
										}
										else{
											print FILE join("\t",$feature_cds->seq_id,1,1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
										}				
									}					
									#Cas de 2 gènes qui se suivent
									elsif($gene{$seq_id}{$gene-1}->end < $feature_cds->start){
										my ($id) = $gene{$seq_id}{$gene-1}->get_tag_values("ID");					
										my (@keys) = keys%{$mRNA{$id}};	
										my ($mrna_id) = $keys[0];
										if ((defined $mrna_id)){
											if (defined $CDS{$mrna_id}[0]){
												my $size = scalar(@{$CDS{$mrna_id}});
												if ($CDS{$mrna_id}[$size-1]->end >= $feature_cds->start+$begin){		
													print FILE join("\t",$feature_cds->seq_id,$CDS{$mrna_id}[$size-1]->end,$feature_cds->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
												}	
												elsif($CDS{$mrna_id}[$size-1]->end < $feature_cds->start+$begin){
													print FILE join("\t",$feature_cds->seq_id,$feature_cds->start-1+$begin,$feature_cds->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
												}
											}
											else{
												if($gene{$seq_id}{$gene-1}->end >= $feature_cds->start+$begin){			
													print FILE join("\t",$feature_cds->seq_id,$gene{$seq_id}{$gene-1}->end,$feature_cds->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
												}				
												elsif($gene{$seq_id}{$gene-1}->end < $feature_cds->start+$begin){
													print FILE join("\t",$feature_cds->seq_id,$feature_cds->start-1+$begin,$feature_cds->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
												}									
											}
										}
										else{
											if($gene{$seq_id}{$gene-1}->end >= $feature_cds->start+$begin){			
												print FILE join("\t",$feature_cds->seq_id,$gene{$seq_id}{$gene-1}->end,$feature_cds->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
											}				
											elsif($gene{$seq_id}{$gene-1}->end < $feature_cds->start+$begin){
												print FILE join("\t",$feature_cds->seq_id,$feature_cds->start-1+$begin,$feature_cds->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
											}									
										}
									}
								}
								elsif ($gene{$seq_id}{$gene-1}->start == $feature_cds->start){
									my $cpt1 = 1;
									while ((exists $gene{$seq_id}{$gene-$cpt1}) && ($gene{$seq_id}{$gene-$cpt1}->start == $feature_cds->start)){
										$cpt1 += 1;
									}
									#Cas du début de chromosome -x
									if(not exists $gene{$seq_id}{$gene-$cpt1-1}){
										#Cas de 2 gènes identique en start et placé en début de chromosome
										if(($feature_cds->start+$begin > 0)){

											#gff2bed:
											print FILE join("\t",$feature_cds->seq_id,$feature_cds->start-1+$begin,$feature_cds->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
										}
										elsif (($feature_cds->start+$begin <= 0)){

											#gff2bed:
											print FILE join("\t",$feature_cds->seq_id,$feature_cds->start-$feature_cds->start,$feature_cds->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
										}
									}
									else{							
										if ($gene{$seq_id}{$gene-$cpt1}->end < $feature_cds->start){
											my ($id) = $gene{$seq_id}{$gene-$cpt1}->get_tag_values("ID");					
											my (@keys) = keys%{$mRNA{$id}};	
											my ($mrna_id) = $keys[0];									
											if (defined $mrna_id){
												if (defined $CDS{$mrna_id}[0]){
													my $size = scalar(@{$CDS{$mrna_id}});	
							
													#Cas de x gènes identique en start et région non codante entre le gene -x et notre gène
													if ($CDS{$mrna_id}[$size-1]->end < $feature_cds->start+$begin){									
														print FILE join("\t",$feature_cds->seq_id,$feature_cds->start-1+$begin,$feature_cds->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
													}
													elsif ($CDS{$mrna_id}[$size-1]->end >= $feature_cds->start+$begin){					
														print FILE join("\t",$feature_cds->seq_id,$CDS{$mrna_id}[$size-1]->end,$feature_cds->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
													}
												}
												else{
													print FILE join("\t",$feature_cds->seq_id,$gene{$seq_id}{$gene-$cpt1}->end,$feature_cds->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";									
												}
											}
											else{
												print FILE join("\t",$feature_cds->seq_id,$gene{$seq_id}{$gene-$cpt1}->end,$feature_cds->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";									
											}
										}
										#Cas de x gènes identique en start dans un autre gène
										elsif ($gene{$seq_id}{$gene-$cpt1}->end >= $feature_cds->start){
											my ($id) = $gene{$seq_id}{$gene-$cpt1}->get_tag_values("ID");					
											my (@keys) = keys%{$mRNA{$id}};	
											my ($mrna_id) = $keys[0];									
											if (defined $mrna_id){
												if (defined $CDS{$mrna_id}[0]){
													my ($cpt)=0;
													my $size = scalar(@{$CDS{$mrna_id}});								
													while (($feature_cds->start	> $CDS{$mrna_id}[$cpt]->start) && ($cpt < $size-1)){
														$cpt += 1;
													}
													if ($CDS{$mrna_id}[0]->end > $feature_cds->start){
														if ($gene{$seq_id}{$gene-$cpt1}->start >= $feature_cds->start+$begin){
															print FILE join("\t",$feature_cds->seq_id,$gene{$seq_id}{$gene-$cpt1}->start,$feature_cds->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
														}		
														elsif($gene{$seq_id}{$gene-$cpt1}->start < $feature_cds->start+$begin){
															print FILE join("\t",$feature_cds->seq_id,$feature_cds->start-1+$begin,$feature_cds->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
														}
													}
													else{
														if ($CDS{$mrna_id}[$cpt-1]->end < $feature_cds->start+$begin){
															print FILE join("\t",$feature_cds->seq_id,$feature_cds->start-1+$begin,$feature_cds->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
														}
														elsif ($CDS{$mrna_id}[$cpt-1]->end >= $feature_cds->start+$begin){									
															print FILE join("\t",$feature_cds->seq_id,$CDS{$mrna_id}[$cpt-1]->start-1,$feature_cds->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
														}
													}
												}
												else{
													print FILE join("\t",$feature_cds->seq_id,1,1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
												}
											}
											else{
												print FILE join("\t",$feature_cds->seq_id,1,1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
											}
										}
									}	
								}
							}				
						}
					}
					else{
						if (defined $feature_cds){							
							if($feature_cds->start+$begin > 0){						
								print FILE join("\t",$feature_cds->seq_id,$feature_cds->start-1+$begin,$feature_cds->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";								
							}
							elsif($feature_cds->start+$begin <= 0){
								print FILE join("\t",$feature_cds->seq_id,$feature_cds->start-$feature_cds->start,$feature_cds->start+$end,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
							}							
						}					
					}
				}
				#Si région 5’ et brin - ou région 3’ et brin + et non transposable_element_gene
				elsif (($feature->strand=~/^-1/) && ($flank =~/prom/) && ($feature->primary_tag !~/transposable_element_gene/)||
				($feature->strand=~/^1/) && ($flank =~/term/) && ($feature->primary_tag !~/transposable_element_gene/)){				
					my ($id) = $gene{$seq_id}{$gene}->get_tag_values("ID");					
					my (@keys) = keys%{$mRNA{$id}};	
					my ($mrna_id) = $keys[0];
					my ($feature_cds);
					if (defined $mrna_id){
						if (defined $CDS{$mrna_id}){
							if (defined $CDS{$mrna_id}[0]){
								my ($size) = scalar(@{$CDS{$mrna_id}});								
								if ($feature->strand =~/^1/){
									$feature_cds = $CDS{$mrna_id}[0];
								}
								else{
									$feature_cds = $CDS{$mrna_id}[$size-1]
								}
							}
						}
					}
					if (0+(keys%{$gene{$seq_id}}) > 1){
						#Cas fin chromosome
						if(not exists $gene{$seq_id}{$gene+1}){					
							if (defined $feature_cds){				
								if ($gene{$seq_id}{$gene-1}->end <= $feature_cds->end) {
									if ($length{$seq_id} > $feature_cds->end-$begin){
										print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
									}
									elsif ($length{$seq_id} <= $feature_cds->end-$begin){
										print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$length{$seq_id}-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
									}
								}
								#cas d'un gène à l'intérieur d'un autre
								elsif ($gene{$seq_id}{$gene-1}->end > $feature_cds->end){ 
									my ($id) = $gene{$seq_id}{$gene-1}->get_tag_values("ID");					
									my (@keys) = keys%{$mRNA{$id}};	
									my ($mrna_id) = $keys[0];														
									if (defined $mrna_id){
										if (defined $CDS{$mrna_id}[0]){										
											if ($feature_cds->end < $CDS{$mrna_id}[0]->start){
												if ($CDS{$mrna_id}[0]->start > $feature_cds->end-$begin){
													print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
												}
												elsif($CDS{$mrna_id}[0]->start <= $feature->end-$begin){
													print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$CDS{$mrna_id}[0]->start-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
												}
											}
											elsif ($feature_cds->start > $CDS{$mrna_id}[0]->end){
												my ($size) = scalar(@{$CDS{$mrna_id}});
												my ($cpt)=$size-1;
												if (($CDS{$mrna_id}[$cpt]->end < $feature_cds->start)||($CDS{$mrna_id}[$cpt]->end > $feature_cds->start)&&($CDS{$mrna_id}[$cpt]->end < $feature_cds->end)){
													if ($gene{$seq_id}{$gene-1}->end > $feature_cds->end-$begin){
														print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
													}
													elsif($gene{$seq_id}{$gene-1}->end <= $feature_cds->end-$begin){
														print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$gene{$seq_id}{$gene-1}->end-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
													}
												}
												elsif (($CDS{$mrna_id}[$cpt]->start > $feature_cds->end)){							
													while (($feature_cds->end < $CDS{$mrna_id}[$cpt]->start) && ($cpt > 0 )){
														$cpt -= 1;
													}
								
													if ($CDS{$mrna_id}[$cpt+1]->start > $feature_cds->end-$begin){
														print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
													}
													elsif($CDS{$mrna_id}[$cpt+1]->start <= $feature->end-$begin){
														print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$CDS{$mrna_id}[$cpt+1]->start-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
													}
												}
											}
										}
										else{
											if($gene{$seq_id}{$gene-1}->end > $feature_cds->end-$begin){
												print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
											}
											elsif($gene{$seq_id}{$gene-1}->end <= $feature_cds->end-$begin){
												print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$gene{$seq_id}{$gene-1}->end-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
											}
										}
									}
									else{
										print FILE join("\t",$feature_cds->seq_id,1,1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
									}
								}
							}
						}
						else{
							#Cas début chromosome					
							if (not exists $gene{$seq_id}{$gene-1}){						
								if (defined $feature_cds){
									#
									if(($gene{$seq_id}{$gene+1}->start <= $feature_cds->end) && ($gene{$seq_id}{$gene+1}->end > $feature_cds->end)){
										my ($id) = $gene{$seq_id}{$gene+1}->get_tag_values("ID");					
										my (@keys) = keys%{$mRNA{$id}};	
										my ($mrna_id) = $keys[0];								
										if (defined $mrna_id){	
											if (defined $CDS{$mrna_id}[0]){
												if ($feature_cds->end < $CDS{$mrna_id}[0]->start){
													if ($CDS{$mrna_id}[0]->start > $feature_cds->end-$begin){
														print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
													}
													elsif($CDS{$mrna_id}[0]->start <= $feature->end-$begin){
														print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$CDS{$mrna_id}[0]->start-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
													}
												}
												elsif ($feature_cds->start > $CDS{$mrna_id}[0]->end){												
													my ($size) = scalar(@{$CDS{$mrna_id}});
													my ($cpt)=$size-1;							
													if (($CDS{$mrna_id}[$cpt]->end < $feature_cds->start)||($CDS{$mrna_id}[$cpt]->end > $feature_cds->start)&&($CDS{$mrna_id}[$cpt]->end < $feature_cds->end)){
														if ($gene{$seq_id}{$gene+1}->end > $feature_cds->end-$begin){
															print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
														}
														elsif($gene{$seq_id}{$gene+1}->end <= $feature_cds->end-$begin){
															print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$gene{$seq_id}{$gene+1}->end-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
														}
													}
													elsif (($CDS{$mrna_id}[$cpt]->start > $feature_cds->end)){							
														while (($feature_cds->end < $CDS{$mrna_id}[$cpt]->start) && ($cpt > 0 )){
															$cpt -= 1;
														}								
														if($CDS{$mrna_id}[$cpt+1]->start > $feature_cds->end-$begin){
															print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
														}
														elsif($CDS{$mrna_id}[$cpt+1]->start <= $feature_cds->end-$begin){
															print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$CDS{$mrna_id}[$cpt+1]->start-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
														}
													}
													elsif (($CDS{$mrna_id}[$cpt]->start < $feature_cds->end) && ($CDS{$mrna_id}[$cpt]->end > $feature_cds->end)){							
														print FILE join("\t",$feature_cds->seq_id,1,1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
													}
												}
											}
											else{
												print FILE join("\t",$feature_cds->seq_id,1,1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
											}
										}
										else{
											print FILE join("\t",$feature_cds->seq_id,1,1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
										}
									}
									#
									elsif($gene{$seq_id}{$gene+1}->start > $feature_cds->end){
										my ($id) = $gene{$seq_id}{$gene+1}->get_tag_values("ID");					
										my (@keys) = keys%{$mRNA{$id}};	
										my ($mrna_id) = $keys[0];				
										if (defined $mrna_id){
											if (defined $CDS{$mrna_id}[0]){
												if ($CDS{$mrna_id}[0]->start <= $feature_cds->end-$begin){
													print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$CDS{$mrna_id}[0]->start-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
												}
												elsif($CDS{$mrna_id}[0]->start > $feature_cds->end-$begin){
													print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
												}
											}
											else{
												if ($gene{$seq_id}{$gene+1}->start <= $feature_cds->end-$begin){
													print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$gene{$seq_id}{$gene+1}->start-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
												}
												elsif($gene{$seq_id}{$gene+1}->start > $feature_cds->end-$begin){
													print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
												}
											}
										}
										else{
											if ($gene{$seq_id}{$gene+1}->start <= $feature_cds->end-$begin){
												print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$gene{$seq_id}{$gene+1}->start-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
											}
											elsif($gene{$seq_id}{$gene+1}->start > $feature->end-$begin){
												print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
											}
										}
									}
								}
							}
							else{						
								if (defined $feature_cds){
									#cas d'un gène à l'intérieur d'un autre -1
									if($gene{$seq_id}{$gene-1}->end > $feature_cds->end){
										my ($id) = $gene{$seq_id}{$gene-1}->get_tag_values("ID");					
										my (@keys) = keys%{$mRNA{$id}};	
										my ($mrna_id) = $keys[0];								
										if (defined $mrna_id){	
											if (defined $CDS{$mrna_id}[0]){
												if ($feature_cds->end < $CDS{$mrna_id}[0]->start){
													if ($CDS{$mrna_id}[0]->start > $feature_cds->end-$begin){
														print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
													}
													elsif($CDS{$mrna_id}[0]->start <= $feature->end-$begin){
														print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$CDS{$mrna_id}[0]->start-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
													}
												}
												else{
													my ($size) = scalar(@{$CDS{$mrna_id}});
													my ($cpt)=$size-1;
													if (($CDS{$mrna_id}[$cpt]->end < $feature_cds->start)||($CDS{$mrna_id}[$cpt]->end > $feature_cds->start)&&($CDS{$mrna_id}[$cpt]->end < $feature_cds->end)){
														if ($gene{$seq_id}{$gene-1}->end > $feature_cds->end-$begin){
															print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
														}
														elsif($gene{$seq_id}{$gene-1}->end <= $feature_cds->end-$begin){
															print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$gene{$seq_id}{$gene-1}->end-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
														}
													}
													elsif (($CDS{$mrna_id}[$cpt]->start > $feature_cds->end)){							
														while (($feature_cds->end < $CDS{$mrna_id}[$cpt]->start) && ($cpt > 0 )){
															$cpt -= 1;
														}
														if ($CDS{$mrna_id}[$cpt+1]->start > $feature_cds->end-$begin){
															print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
														}
														elsif($CDS{$mrna_id}[$cpt+1]->start <= $feature_cds->end-$begin){
															print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$CDS{$mrna_id}[$cpt+1]->start-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
														}
													}
												}
											}
											else{
												print FILE join("\t",$feature_cds->seq_id,1,1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
											}
										}
										else{
											print FILE join("\t",$feature_cds->seq_id,1,1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
										}
									}												
									#Cas d'un gène chevauchant 2 gènes
									elsif($gene{$seq_id}{$gene-1}->end <= $feature_cds->end){
										if($gene{$seq_id}{$gene+1}->start <= $feature_cds->end){ 
											if($gene{$seq_id}{$gene+1}->end > $feature_cds->end){
												my ($id) = $gene{$seq_id}{$gene+1}->get_tag_values("ID");					
												my (@keys) = keys%{$mRNA{$id}};	
												my ($mrna_id) = $keys[0];						
												if (defined $mrna_id){	
													if (defined $CDS{$mrna_id}[0]){
														if ($feature_cds->end < $CDS{$mrna_id}[0]->start){
															if ($CDS{$mrna_id}[0]->start > $feature_cds->end-$begin){
																print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
															}
															elsif($CDS{$mrna_id}[0]->start <= $feature->end-$begin){
																print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$CDS{$mrna_id}[0]->start-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
															}
														}
														elsif ($feature_cds->start > $CDS{$mrna_id}[0]->end){
															my ($size) = scalar(@{$CDS{$mrna_id}});
															my ($cpt)=$size-1;
															if (($CDS{$mrna_id}[$cpt]->end < $feature_cds->start)||($CDS{$mrna_id}[$cpt]->end > $feature_cds->start)&&($CDS{$mrna_id}[$cpt]->end < $feature_cds->end)){
																if ($gene{$seq_id}{$gene+1}->end > $feature_cds->end-$begin){
																	print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
																}
																elsif($gene{$seq_id}{$gene+1}->end <= $feature_cds->end-$begin){
																	print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$gene{$seq_id}{$gene+1}->end-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
																}
															}
															elsif (($CDS{$mrna_id}[$cpt]->start > $feature_cds->end)){							
																while (($feature_cds->end < $CDS{$mrna_id}[$cpt]->start) && ($cpt > 0 )){
																	$cpt -= 1;
																}
																if ($CDS{$mrna_id}[$cpt+1]->start > $feature_cds->end-$begin){
																	print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
																}
																elsif($CDS{$mrna_id}[$cpt+1]->start <= $feature->end-$begin){
																	print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$CDS{$mrna_id}[$cpt+1]->start-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
																}
															}
														}
													}
													else{
														print FILE join("\t",$feature_cds->seq_id,1,1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";											
													}											
												}
												else{
													print FILE join("\t",$feature_cds->seq_id,1,1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";											
												}
											}
											elsif($gene{$seq_id}{$gene+1}->end <= $feature_cds->end){
												my $cpt1 = 1;
												while ((exists $gene{$seq_id}{$gene+$cpt1}) && ($gene{$seq_id}{$gene+$cpt1}->end <= $feature_cds->end)){ #attention
													$cpt1 += 1;
												}
												if(not exists $gene{$seq_id}{$gene+$cpt1+1}){													
													if ($length{$seq_id} > $feature_cds->end-$begin){
														print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
													}
													elsif ($length{$seq_id} <= $feature_cds->end-$begin){
														print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$length{$seq_id}-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
													}
												}
												elsif(($gene{$seq_id}{$gene+$cpt1+1}->start <= $feature_cds->end) && ($gene{$seq_id}{$gene+$cpt1+1}->end > $feature_cds->end)){												
													my ($id) = $gene{$seq_id}{$gene+$cpt1+1}->get_tag_values("ID");					
													my (@keys) = keys%{$mRNA{$id}};	
													my ($mrna_id) = $keys[0];								
													if (defined $mrna_id){	
														if (defined $CDS{$mrna_id}[0]){
															if ($feature_cds->end < $CDS{$mrna_id}[0]->start){
																if ($CDS{$mrna_id}[0]->start > $feature_cds->end-$begin){
																	print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
																}
																elsif($CDS{$mrna_id}[0]->start <= $feature->end-$begin){
																	print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$CDS{$mrna_id}[0]->start-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
																}
															}
															elsif ($feature_cds->start > $CDS{$mrna_id}[0]->end){											
																my ($size) = scalar(@{$CDS{$mrna_id}});
																my ($cpt)=$size-1;							
																if (($CDS{$mrna_id}[$cpt]->end < $feature_cds->start)||($CDS{$mrna_id}[$cpt]->end > $feature_cds->start)&&($CDS{$mrna_id}[$cpt]->end < $feature_cds->end)){
																	if ($gene{$seq_id}{$gene+$cpt1+1}->end > $feature_cds->end-$begin){
																		print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
																	}
																	elsif($gene{$seq_id}{$gene+$cpt1+1}->end <= $feature_cds->end-$begin){
																		print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$gene{$seq_id}{$gene+$cpt1+1}->end-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
																	}
																}
																elsif (($CDS{$mrna_id}[$cpt]->start > $feature_cds->end)){							
																	while (($feature_cds->end < $CDS{$mrna_id}[$cpt]->start) && ($cpt > 0 )){
																		$cpt -= 1;
																	}								
																	if($CDS{$mrna_id}[$cpt+1]->start > $feature_cds->end-$begin){
																		print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
																	}
																	elsif($CDS{$mrna_id}[$cpt+1]->start <= $feature_cds->end-$begin){
																		print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$CDS{$mrna_id}[$cpt+1]->start-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
																	}
																}
																elsif (($CDS{$mrna_id}[$cpt]->start < $feature_cds->end) && ($CDS{$mrna_id}[$cpt]->end > $feature_cds->end)){							
																	print FILE join("\t",$feature_cds->seq_id,1,1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
																}
															}
														}
														else{
															print FILE join("\t",$feature_cds->seq_id,1,1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
														}
													}
													else{
														print FILE join("\t",$feature_cds->seq_id,1,1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
													}												
												}
												elsif($gene{$seq_id}{$gene+$cpt1+1}->start > $feature_cds->end){												
													my ($id) = $gene{$seq_id}{$gene+$cpt1+1}->get_tag_values("ID");					
													my (@keys) = keys%{$mRNA{$id}};	
													my ($mrna_id) = $keys[0];								
													if (defined $mrna_id){	
														if (defined $CDS{$mrna_id}[0]){
															if ($CDS{$mrna_id}[0]->start > $feature_cds->end-$begin){
																print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
															}
															elsif($CDS{$mrna_id}[0]->start <= $feature_cds->end-$begin){
																print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$CDS{$mrna_id}[0]->start-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
															}
														}
														else{
															if ($gene{$seq_id}{$gene+$cpt1+1}->start <= $feature_cds->end-$begin){
																print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$gene{$seq_id}{$gene+$cpt1+1}->start-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
															}
															elsif($gene{$seq_id}{$gene+$cpt1+1}->start > $feature_cds->end-$begin){
																print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
															}
														}
													}
													else{
														if ($gene{$seq_id}{$gene+$cpt1+1}->start <= $feature_cds->end-$begin){
															print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$gene{$seq_id}{$gene+$cpt1+1}->start-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
														}
														elsif($gene{$seq_id}{$gene+$cpt1+1}->start > $feature_cds->end-$begin){
															print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
														}
													}
												}									
											}
										}
										#cas d'un gène à l'intérieur d'un autre
										elsif($gene{$seq_id}{$gene+1}->start > $feature_cds->end){
											my ($id) = $gene{$seq_id}{$gene+1}->get_tag_values("ID");					
											my (@keys) = keys%{$mRNA{$id}};	
											my ($mrna_id) = $keys[0];
											if (defined $mrna_id){
												if (defined $CDS{$mrna_id}[0]){
													if($CDS{$mrna_id}[0]->start <= $feature_cds->end-$begin){
														print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$CDS{$mrna_id}[0]->start-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
													}
													elsif($CDS{$mrna_id}[0]->start > $feature_cds->end-$begin){
														print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
													}
												}
												else{
													if($gene{$seq_id}{$gene+1}->start <= $feature_cds->end-$begin){
														print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$gene{$seq_id}{$gene+1}->start-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
													}
													elsif($gene{$seq_id}{$gene+1}->start > $feature_cds->end-$begin){
														print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
													}
												}
											}
											else{
												if($gene{$seq_id}{$gene+1}->start <= $feature_cds->end-$begin){
													print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$gene{$seq_id}{$gene+1}->start-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
												}
												elsif($gene{$seq_id}{$gene+1}->start > $feature->end-$begin){
													print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
												}
											}
										}
									}
								}
							}
						}
					}
					else{
						if (defined $feature_cds){				
							if ($length{$seq_id} > $feature_cds->end-$begin){
								print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$feature_cds->end-$begin,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
							}
							elsif ($length{$seq_id} <= $feature_cds->end-$begin){
								print FILE join("\t",$feature_cds->seq_id,$feature_cds->end-$end-1,$length{$seq_id}-1,$feature->{_gsf_tag_hash}->{ID}->[0]."_".$title3,$score,$strand)."\n";
							}
						}
					}
				}				
			}
		}
	}
	close FILE;
	
	system("bedtools getfasta -s -fi $genome -bed $outfile -fo $fna_output_file -name");

}


# Script options
#################

=pod

=head1 OPTIONS

#--- describes parameters given to the script
#+++ command line syntax
#--- requirement of the option and its parameter can be:
#--- required: name or nature inside <>
#--- optional: name or nature inside []
#--- alternative between 2 elements: elements separated by a |

=head2 Parameters

=over 4

=item B<[option_name]> ([option nature]): #+++

[option description]. #+++
Default: [option default value if one] #+++

=back
#--- Example:
#---
#--- Template.pl [-help | -man]
#---
#--- Template.pl [-debug [debug_level]] [-size <width> [height]]
#---
#--- =over 4
#---
#--- =item B<-help>:
#---
#--- Prints a brief help message and exits.
#---
#--- =item B<-man>:
#---
#--- Prints the manual page and exits.
#---
#--- =item B<-debug> (integer):
#---
#--- Executes the script in debug mode. If an integer value is specified, it will
#--- be the debug level. If "-debug" option was used without specifying a debug
#--- level, level 1 is assumed.
#--- Default: 0 (not in debug mode).
#---
#---=item B<-size> (positive_real) (positive_real):
#---
#--- Set the dimensions of the object that will be drawn. The first value is
#--- the width; the height is the second value if specified, otherwise it will
#--- assume height and width are equal to the first value.
#--- Default: width and height are set to 1.
#---
#---=back

=cut

# CODE START
#############

# options processing
my ($man, $help, $debug, $sort);

# parse options and print usage if there is a syntax error.
GetOptions("help|?"   => \$help,
           "man"      => \$man,
           "debug"    => \$debug, # a flag
           "sort=s"	  => \$sort,
           "bed=s"      => \$bed,
           "flank=s"   => \$flank,
           "begin=i" => \$begin,
           "end=i" => \$end)
#           "length=i" => \$length, # numeric
#           "file=s"   => \$data) # a string
    or pod2usage(2);
if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}


print "Looking for files in $gff3_file\n";
my @files = $gff3_file;
print "Found " . @files . " files to process...\n";

#Fonction pour trier et faire des stats:
if (defined($sort) && $sort eq "on"){
	sortGff3();
	gff3_stat();
}

if (defined($bed)){
	gff3tobed();
}

if (defined($begin) && defined($end) && ($begin == 0) && ($end == 0)){
	bedtools();
}
elsif (defined($begin) && defined($end) && defined($flank)){
	flanking_region();
}
else{
	bedtools();
}

# CODE END
###########


=pod

=head1 DIAGNOSTICS

#--- Give and explain here every error message the the script may generate.
#--- Use: (W)=warning (optional), (D)=deprecation (optional),
#---      (S)=severe warning (mandatory), (F)=fatal error (trappable),
#---      (X)=very fatal error (non-trappable).
#--- example:
#---
#--- =over 4
#---
#--- =item *
#---
#--- Negative radius;
#--- (F) Can not draw a circle with a negative radius.
#---
#---=back #+++

=head1 AUTHORS

Jonathan LORENZO (CIRAD), jonathan.lorenzo@cirad.fr

[author1_name (company), email]#+++

[author2_name (company), email]#+++

=head1 VERSION

Version [version.subversion.build] #+++

Date [DD/MM/YYY] #+++

=head1 SEE ALSO

#--- other documentation or objects related to this package
[perl(1), L<Geometry::Square>] #+++

=cut
