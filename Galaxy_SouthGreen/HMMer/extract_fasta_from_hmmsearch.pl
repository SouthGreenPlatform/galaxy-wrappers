#!/usr/bin/perl
#use lib '/usr/lib/perl5/vendor_perl/5.8.8/';
use Bio::SearchIO;  
use Getopt::Long;
my $hmmer_outfile;
my $reference_fastafile;
my $output;
my $type;
my $value;
GetOptions(
           'hmm=s'         => \$hmmer_outfile,
           'fasta=s'       => \$reference_fastafile,  
           'type=s'          => \$type, 
           'value=s'        => \$value
          )  or pod2usage(1);  
my $in = new Bio::SearchIO(
			   -format => 'hmmer3',
			   -file   => $hmmer_outfile
			  );

my @list_id;
my $cpt = 0;
while( my $result = $in->next_result ) {
    my $query_name = $result->query_accession(); 
    while(my $hit = $result->next_hit()) { 
    	if ($type eq "number") {
    		$cpt++;
    		push @list_id , $hit->name(); 
    		last if $cpt == $value;
    	}
    	else {
			if ($hit->significance() < $value){ 
	    		push @list_id , $hit->name(); 
    		}
    	}
	}
}
my %hashTemp = map { $_ => 1 } @list_id;
@list_id = sort keys %hashTemp;



my $out = $hmmer_outfile.'.txt';
my $out_fasta = $hmmer_outfile.'.fa';
open (handle, ">$out");
print handle join ("\n", @list_id), "\n";
close handle;
my $path =  $ENV{HOME};
system("$path/galaxy/tools/SouthGreen/HMMer/fasta-make-index $reference_fastafile -f");
system("$path/galaxy/tools/SouthGreen/HMMer/fasta-fetch $reference_fastafile -f $out ");





