#!/usr/bin/perl

use Bio::SeqIO;
use Bio::Tools::SeqStats;

my $file_input = shift;
my $file_output = shift;
my $seq = Bio::SeqIO->new(
			  -file => $file_input,
			  -format => "fasta"
			  )->next_seq();

my $out_genomic_seq = Bio::SeqIO->new(
				      -file   => ">$file_output",
				      -format => 'fasta' 
				      );
my $seq_stats  =  Bio::Tools::SeqStats->new(
					    -seq => $seq
                                            ); 
$seq->seq(uc($seq->seq));    
$out_genomic_seq->write_seq($seq);
my $hash_ref = $seq_stats->count_monomers();
my $clone_name = $seq->display_id;
my $length = $seq->length;
my $gc = ($hash_ref->{'G'} + $hash_ref->{'C'}) / $length * 100;
my $gc_content = sprintf("%.2f",$gc);
print "seq_name ", $seq->display_id, "\n";
print "seq_length ", $seq->length, "\n";
print "gc_content ", $gc_content, "\n";
foreach my $base (sort keys %{$hash_ref}) {
    print $base," ",$hash_ref->{$base},"\n";
}
    
