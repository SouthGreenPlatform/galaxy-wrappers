#!/usr/bin/perl
use FindBin;                 
use lib "$FindBin::Bin/lib";
use Bio::SearchIO;
use Data::Dumper;
my $file = shift;
my $output = shift;

my $result = new Bio::SearchIO(
	-file => $file,
	-format => "blast"
)->next_result;
  
my @T;
my $cpt = 0;
my $query_length = $result->query_length;
while(my $hit = $result->next_hit()){
    my $hit_name = $hit->name; 
    my $hit_length = $hit->length();
    my $description = $hit->description();
    while(my $hsp = $hit->next_hsp()) {
		my $query_start    = $hsp->start('query');
		my $query_end      = $hsp->end('query');
		my $hit_start      = $hsp->start('hit');
		my $hit_end        = $hsp->end('hit');
		my $score          = $hsp->score();
		my $evalue         = $hsp->evalue();
		my $frac_identical = sprintf("%.f",$hsp->frac_identical() * 100);
		my $frac_conserved = sprintf("%.f",$hsp->frac_conserved() * 100);
		my $delta          = ($query_end > $query_start) ? $query_end - $query_start : $query_start - $query_end;
		my $pc_s           = int (( $hit_end - $hit_start + 1) / $query_length * 100);
		my $pc_q           = int (( $delta + 1) / $query_length * 100);
		my $strand         = $hsp->strand('query') == 1 ? "+" : "-";
		my $frame          = $strand . ($hsp->frame('query') + 1) ;
		$evalue            =~ /^e/ ? "1".$evalue : $evalue;
		($query_start,$query_end) = ($query_end,$query_start) if $strand eq "-";
		$T[$cpt++] = join(" ",$query_start, $query_end, $score, $evalue, $frame, $hit_name, $hit_start, $hit_end); 
    }
}

my @Tsort = sort sort4eugene(@T);
open(OUT,">$output");
foreach my $line (@Tsort) {
    print OUT $line ,"\n";
}
close OUT;


sub sort4eugene {
    my(@A) = split(' ',$a);
    my(@B) = split(' ',$b);
    if ( $A[5] ne $B[5] ) {
		return $A[5] cmp $B[5];
    }
    else {
		return $A[6] <=> $B[6];
    }
}

