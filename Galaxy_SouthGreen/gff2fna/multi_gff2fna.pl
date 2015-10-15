#!/usr/bin/perl

=pod

=head1 NAME

[multi_gff2fna.pl - short description]

=head1 SYNOPSIS

    perl multi_gff2fna.pl <gff3_input_file> <genome_sequence_file> <fna_output_file> <bed_output_file>  <input_multifasta> prom|term -begin=0|-x|+x -end=0|-x|+x

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
use List::MoreUtils qw(any);

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

my $multi_fna=$ARGV[4];

my $data = read_file($multi_fna);
$data =~ s/\r/\n/g;
write_file($multi_fna, $data);
 
sub read_file {
    my ($filename) = @_;
 
    open my $in, '<:encoding(UTF-8)', $filename or die "Could not open '$filename' for reading $!";
    local $/ = undef;
    my $all = <$in>;
    close $in;
 
    return $all;
}
 
sub write_file {
    my ($filename, $content) = @_;
 
    open my $out, '>:encoding(UTF-8)', $filename or die "Could not open '$filename' for writing $!";;
    print $out $content;
    close $out;
 
    return;
}


my $reader=new Bio::SeqIO(-format=>'fasta',-file=>$multi_fna);

my @list_gene;
while (my $seqRec=$reader->next_seq()){
	my $id = $seqRec->id;
	if ($id =~ /(.+\d+)[A-Z](\d+.*)_[A-Z]{5}/i){
		my $gene = $1."g".$2;
		#print $gene."\n";
		push (@list_gene, $gene);		
	}
	elsif ($id =~ /(.+)_[A-Z]{5}/i){
		my $gene = $1;
		push (@list_gene, $gene);
	}
}

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
my $flank = $ARGV[5];
my $begin;
my $end;

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
	
	my (%gene, %mRNA, %polypeptide, %CDS, %exon, %five_prime_UTR, $mrna_id, $id, %transposable_region);
	
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
		elsif ($feature->primary_tag =~ /polypeptide|protein/ ) { 
			($mrna_id) = $feature->get_tag_values("Derives_from"); 
			push @{$polypeptide{$mrna_id}}, $feature;
		}
		elsif ($feature->primary_tag eq "exon") { 	
			if ($feature->has_tag("Parent")){ 
			($mrna_id) = $feature->get_tag_values("Parent");  
			push @{$exon{$mrna_id}}, $feature;
			}
		}
		elsif ($feature->primary_tag eq "CDS") {
			if ($feature->has_tag("Parent")){
			($mrna_id) = $feature->get_tag_values("Parent"); 
			push @{$CDS{$mrna_id}}, $feature;
			}
		}
		elsif ($feature->primary_tag =~/five_prime_UTR/i) { 
			if ($feature->has_tag("Parent")){	
			($mrna_id) = $feature->get_tag_values("Parent"); 
			push @{$five_prime_UTR{$mrna_id}}, $feature;
			}
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
	
	chmod 0777, $outfile;
	
	foreach my $seq_id (sort {$a cmp $b} keys%gene){
		foreach my $gene (sort {$a <=> $b} keys%{$gene{$seq_id}}){
			
			my ($mrna_id);
			my ($id) = $gene{$seq_id}{$gene}->get_tag_values("ID");
			
			my @keys = keys%{$mRNA{$id}};	
			$mrna_id = $keys[0];
			
			my ($feature) = $gene{$seq_id}{$gene};
			my ($feature_pol) = $polypeptide{$mrna_id}[0];
			my $score = ".";
			my $strand;
			my $element;
			
			if($feature->strand=~/-/){
				$strand='-';
			}else{
				$strand='+';
			}
			
			# if ($list_gene[0] =~/[A-Z]{2}\d+G\d+/i){
# 				if ($feature->has_tag("Name")){
# 					$element = $feature->{_gsf_tag_hash}->{Name}->[0];
# 				}else{
# 					$element = $feature->{_gsf_tag_hash}->{ID}->[0];
# 				}
# 			}
# 			else{
			if ($feature->seq_id =~/MUSAC\d+/){				
				$element = $feature_pol->{_gsf_tag_hash}->{Derives_from}->[0];
				# my @count1 = grep {/$element/i} @list_gene;
# 				my $count1 = @count1;
# 				if ($count1 == 0){
# 					$element = $feature_pol->{_gsf_tag_hash}->{Name}->[0];
# 				}				
			}
			elsif ($feature->seq_id =~/COFCA\d+/){
				$element = $feature->{_gsf_tag_hash}->{Name}->[0];				
			}
			else{
				if (defined $feature_pol){
					if ($feature_pol->has_tag("Name")){
						if ($feature_pol->{_gsf_tag_hash}->{Name}->[0] =~/(.+)-.+/){
							$element = $1;
						}
						else{
							$element = $feature_pol->{_gsf_tag_hash}->{Name}->[0];
						}					
					}else{
						$element = $feature_pol->{_gsf_tag_hash}->{ID}->[0];
					}
				}
				else{
					if ($feature->has_tag("Name")){
						$element = $feature->{_gsf_tag_hash}->{Name}->[0];
					}else{
						$element = $feature->{_gsf_tag_hash}->{ID}->[0];
					}
				}
			}
			#}
			
			my @count = grep {/$element/i} @list_gene;
			my $count = @count;
			#for my $str (@count) {
       		#	print "$str\n";
   			#}
			
			if ($count >=1) {
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
				#Si région 5’ et brin - ou région 3’ et brin + et non transposable_element_gene
				elsif (($flank =~/prom/) && ($feature->primary_tag !~/transposable_element_gene/)
				|| ($feature->strand=~/^1/) && ($flank =~/term/) && ($feature->primary_tag !~/transposable_element_gene/)){				
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
					#Cas fin chromosome
					if(not exists $gene{$seq_id}{$gene+1}){					
						if ((defined $feature_cds)&&(defined $gene{$seq_id}{$gene-1})){				
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
my ($man, $help, $debug);

# parse options and print usage if there is a syntax error.
GetOptions("help|?"   => \$help,
           "man"      => \$man,
           "debug"    => \$debug, # a flag
           "begin=i" => \$begin,
           "end=i" => \$end)
#           "length=i" => \$length, # numeric
#           "file=s"   => \$data) # a string
    or pod2usage(2);
if ($help) {pod2usage(0);}
if ($man) {pod2usage(-verbose => 2);}


#print "Looking for files in $gff3_file\n";
#my @files = $gff3_file;
#print "Found " . @files . " files to process...\n";

flanking_region();

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
