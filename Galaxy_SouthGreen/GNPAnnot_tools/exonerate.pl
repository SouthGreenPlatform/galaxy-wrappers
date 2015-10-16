#!/usr/bin/perl

use Bio::SeqIO; 

my $input_exonerate  = shift;
my $source   		 = shift;
my $database         = shift;
my $file_fasta       = shift;
my $output_exonerate = shift;

my $in = new Bio::SeqIO(
	-file => $file_fasta,
	-format => 'fasta'
);
while(my $seqobj = $in->next_seq) {
    my $display_id = (split(/\-/,$seqobj->display_id))[0]; 
    $seq_region{$display_id} = $seqobj;
}
$in->close;
my %result; 

open(IN,$input_exonerate);

while(<IN>) {
    chomp;
    my ($gene_id,$est_id,$frame,$description);
    my $file_est;
    my $file_est_transeq;
    my $file_est_transeq1;
    my $file_region;
    if ($source eq "nuc") {
		($gene_id,$est_id,$frame,$description) = (split(/\t/,$_));
		if (exists $exist{$gene_id}{$est_id}){next;}
		else {
			$exist{$gene_id}{$est_id} = 1;
			$file_est = $gene_id .".mrna.fna";
			$file_est_transeq = $gene_id .".mrna.faa";
			$file_est_transeq1 = $gene_id .".mrna1.faa";
			$file_region = $gene_id.".region.fna";
			system("fasta-fetch $database '$est_id' > $file_est");
			system("transeq -sequence $file_est -outseq $file_est_transeq -frame $frame -alternative true 2>/dev/null");
			my $seq = Bio::SeqIO->new(
				-file => $file_est_transeq,
				-format => "fasta"
			)->next_seq();
	
			$seq->display_id($est_id);
			my $out = Bio::SeqIO->new(
				-file => ">$file_est_transeq1",
				-format => "fasta"
			);
			$seq->desc('');
			$out->write_seq($seq);
		}
    }
    else {
		($gene_id,$est_id,$frame,$description) = (split(/\t/,$_));
		if (exists $exist{$gene_id}{$est_id}){next;}
		else {
			$exist{$gene_id}{$est_id} = 1;
			$gene_id = (split(/-/,$gene_id))[0];	
			$file_est = $gene_id .".mrna.fna";
			$file_est_transeq = $gene_id .".faa";
			$file_est_transeq1 = $gene_id .".1.faa";
			$file_region = $gene_id.".region.fna";
			system("fasta-fetch  $database '$est_id' > $file_est_transeq1");
			 
		}
    }
    if (defined $seq_region{$gene_id}) {
		my $seq_out = new Bio::SeqIO(
			-file => ">$file_region",
			-format=>"fasta"
		);
		$seq_out->write_seq($seq_region{$gene_id});
		system("exonerate --model protein2genome $file_est_transeq1 $file_region >> $output_exonerate");
		unlink $file_est if -e $file_est;
		unlink $file_est_transeq if -e $file_est_transeq;
		unlink $file_est_transeq1 if -e $file_est_transeq1;
		unlink $file_region if -e $file_region;
    }
}
close IN;
