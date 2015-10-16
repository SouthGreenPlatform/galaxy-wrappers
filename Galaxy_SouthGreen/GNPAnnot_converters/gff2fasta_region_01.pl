#!/usr/bin/perl 
use Carp qw (cluck confess croak); 
use Pod::Usage;
use Error;
use Fatal qw (open close); 
use Getopt::Long;
use Bio::Tools::GFF;
use Bio::SeqIO;


my $usage = q/
gff2fasta_region.pl [--gff <file.gff>] [--fasta <file_fasta>] [--output <outfile>] [--verbose] [--help]

Parameters
--gff        	file gff [required]
--fasta         fasta file [required]
--output        output file. Default file.gff.fna 
--verbose
--help
/;

my $file_gff;
my $file_fasta;
my $file_fasta_region;
my $help;
my $verbose;

GetOptions(
           'gff=s'    => \$file_gff, 
           'fasta=s'  => \$file_fasta,  
           'output=s' => \$file_fasta_region,	
           'verbose'  => \$verbose,
           'help'     => \$help
) or die $usage;

if($help){
    print $usage;
    exit 0;
} 
if ($file_gff eq "") {
    print "Warn :: --gff is empty\nPlease specify a valid gff file\n";
    print $usage;
    exit 0;
}
if ($file_fasta eq "") {
    print "Warn :: --fasta is empty\nPlease specify a valid gff file\n";
    print $usage;
    exit 0;
}
my $seqin = new Bio::SeqIO(
	-file   => $file_fasta,
	-format => 'fasta'
);
my $gffin = new Bio::Tools::GFF(
	-file   	 => $file_gff,
	-gff_version => 3
);

$file_fasta_region = "$file_gff.fna" if $file_fasta_region eq "";
my $seqout = Bio::SeqIO->new(
	-file   => ">$file_fasta_region",
	-format => 'fasta'
); 

my %sequence;
print "Loading sequence information\n" if $verbose;
while( my $seqobj = $seqin->next_seq) { 
	$sequence{$seqobj->display_id} = $seqobj;
}
$seqin->close; 

my %count;
my %gff;
print "Loading gff information\n" if $verbose;
while (my $feature = $gffin->next_feature){
	if ($feature->primary_tag() eq "gene") {	
		push @{$gff{$feature->seq_id}},$feature;
		$count{$feature->seq_id}++;
	}
}
$gffin->close;

print "Extract sequence\n" if $verbose;

foreach my $seq_id (sort {$a cmp $b} keys %gff){ 
	if ($sequence{$seq_id}){
		my $seqobj = $sequence{$seq_id};
		my $end_chr = $seqobj->length(); 
		my $overlap   = 5000;
		my $num_gene = $count{$seq_id};
		print "Run on $seq_id\n" if $verbose;
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
			my $query = join( "-", $locus_tag, $start_of_region, $end_of_region );
			my $seq_region = $seqobj->subseq( $start_of_region, $end_of_region );
			my $seq_region_obj    = Bio::PrimarySeq->new(
				-display_id => $query,
				-seq        => $seq_region
			);
			$seqout->write_seq($seq_region_obj); 
   		#	print progress_bar( $cpt, $num_gene, 25, '=' );
		}
	}
   	print "\n";
}
$seqout->close;

print "Done...\n" if $verbose;
print "Result in : $file_fasta_region\n";



sub progress_bar {
    my ( $got, $total, $width, $char ) = @_;
    $width ||= 25; $char ||= '=';
    my $num_width = length $total;
    sprintf "|%-${width}s| Done %${num_width}s genes of %s (%.2f%%)\r",
    $char x (($width-1)*$got/$total). '>',
    $got, $total, 100*$got/+$total;
}
