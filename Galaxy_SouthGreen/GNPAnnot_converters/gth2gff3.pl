#!/usr/bin/perl

use Bio::Tools::GFF;
use File::Basename;
my $file_gth = shift; 
my $tool_data_path = shift;
my $pci      = shift;
my $pcs      = shift;
my $output_gff = shift;
my $output_gth_eugene = shift;

my $gffio = new Bio::Tools::GFF(
	-file => ">$output_gff",
	-gff_version => 3
);
my $tag;
open(GENESEQER,$file_gth) or die "ERROR $0 >$file_gth<\n";
open(EST,">$output_gth_eugene") or die "ERROR $0 >$file_out<\n";
my $exonok = 0;
my $seqok = 0;
my $match_len = 0;
my $curr_pci = 0;
my $curr_pcs = 0;
my @a_exon=();
my $strand_cdna=0;
my $strand_gen=0;
my $mrna_id = "";
my $est_len=0;
my %h_mrna=();
my $bank;
my %exist;
my $count_hit;
my $database;
while(my $lign=<GENESEQER>)  {
    chomp($lign);
    next if ( $lign =~ /^\s*$/ );
    next if ( $lign =~ /Intron/ );
    if ( $lign =~ /^Predicted gene structure/ || $lign =~ /^Genomic/) {
		($exonok,$seqok) = (1,0);
		next;
    }
    if ( $lign =~ /^MATCH\s+/ ) {
		my @a_F = split(' ',$lign);
		$curr_pci = $a_F[3] * 100;
		$curr_pcs = $a_F[5] * 100;
		$match_len = $a_F[4];
		$scaffold = $a_F[1];
		$scaffold =~ s/[\-\+]$//;
		($gene_name,$chr,$extend_start,$extend_end) = (split(/\-/,$scaffold));
		$scaffold = $chr;
		if ( $mrna_id eq '' ) {
	    	$mrna_id = $a_F[2];
	    	$mrna_id =~ s/[\-\+]$//;
		}
		($strand_cdna) = ( $a_F[2] =~ /\-$/ )  ? 1 : 0;
		($strand_gen) = ( $a_F[1] =~ /\-$/ )  ? 1 : 0;
    }
    if ( $lign =~ /^Alignment.+:\s*$/ ) {
		if ( ($curr_pci >= $pci) && ($curr_pcs >= $pcs)  && $#a_exon != -1 ) {
	    	my ($lib,$name) = (split(/\|/,$mrna_id))[0,1];
	    	my $est_name = join("_",$lib,$name);
	    	my $so = "EST_match";
	    	my ($min,$max,@a_exon) = &FilterReformatA6($match_len,$strand_gen,$mrna_id, $strand_cdna,$est_len,$tag,@a_exon);
	    	$min = $min  + $extend_start - 1;
	    	$max =  $max + $extend_start - 1;
	    	$h_mrna{$mrna_id} = "$match_len $strand_cdna $strand_gen $min $max";
			my $name1 = join("_",$est_name,$gene_name);
			$count_hit = 1 + $exist{$name1}++; 
			my $gene_id = join(	 "_",  $name1, "gene",  sprintf('%04d', $count_hit) );
			my $mrna_id = join(	 "_",    $name1, "mrna",  sprintf('%04d', $count_hit) ); 
	    	my $id = join("_", $database, $est_name);
	    	my $target = join( "+",$est_name,$min,$max );
	    	my $strand = $strand_gen == 0 ? "+" : "-";
	    	my $match =  new Bio::SeqFeature::Generic(
				-seq_id      => $gene_name,
				-source_tag  => $database,
				-primary_tag => "EST_match",
				-start       => $min,
				-end         => $max,
				-strand      => $strand,
				-tag         => {
					ID => $gene_id,
					Name => $gene_id,
					Target => $target 
				}
			);
	    	$match->add_tag_value("Alias",$alias) if $alias;
	   		$gffio->write_feature($match); 	
	    	for (my $i=0;$i<=$#a_exon;$i++) {
				my ($gen_begin,$gen_end,$real_strand,$mrna_id1,$cdna_begin,$cdna_end) = (split(/ /,$a_exon[$i]));	
	    		$gen_begin = $gen_begin + $extend_start - 1;
	    		$gen_end = $gen_end  + $extend_start - 1;
				my $match_part = new Bio::SeqFeature::Generic(
					-seq_id      => $gene_name,
					-source_tag  => $database,
					-primary_tag => 'match_part',
					-start       => $gen_begin,	
					-end         => $gen_end,
					-strand      => $strand,
					-tag         => {
						Parent => $gene_id	
					}
				);
				$gffio->write_feature($match_part);
	    	}
	    	print EST @a_exon;				
		}
		$match_len = 0;
		@a_exon=();
		$exonok = 0;
		$strand_gen=0;
		$strand_cdna=0;
		$mrna_id='';
		$est_len=0;
		next;
    }
    if ( $seqok )  {
		my ($seq) = $lign;
		$seq =~ s/[^A-Za-z]//g;
		$est_len+=length($seq);
    }
    if ( $lign =~ /^EST sequence/i ) {
		($mrna_id) = $lign =~ /description=(\S+)/;
		($bank) = $lign =~ /file=(\S+),/;
		$database = (split(/\//,$bank))[-1];
		$database = (split(/\./,$database))[0];
		$seqok = 1;
		next;
    }
    next if ( ! $exonok );
    push @a_exon, $lign if ( $lign =~ /^\s*Exon/ );
}
close(GENESEQER);
close(EST);
$gffio->close;

sub FilterReformatA6 {
    my ($match_len,$strand_gen,$mrna_id,$strand_cdna,$est_len,$tag,@a_exon) = @_;  
    my $min=99999999999999999999999;
    my $max=-1;
    my ($real_strand)  = $strand_gen ; 
    for(my $i=0;$i<=$#a_exon;$i++) {
        my ($gen_begin,$gen_end) = $a_exon[$i] =~ /Exon\s+\d+\s+(\d+)\s+(\d+)/;
        my ($cdna_begin,$cdna_end) = $a_exon[$i] =~ /cDNA\s+(\d+)\s+(\d+)/;
        ($gen_begin,$gen_end)   = ( $strand_gen == 0 ) ? ($gen_begin,$gen_end) : ($gen_end,$gen_begin);
        ($cdna_begin,$cdna_end)= (($est_len-$cdna_end+1),($est_len-$cdna_begin+1)) if ( $strand_gen == 1 );
        ($cdna_begin , $cdna_end) = ( $cdna_begin > $cdna_end ) ? ($cdna_end,$cdna_begin) : ($cdna_begin , $cdna_end); # A=6 format
        $a_exon[$i] = "$gen_begin $gen_end $match_len 0 $real_strand $tag|$mrna_id $cdna_begin $cdna_end\n";
		$a_exon_tab[$i] = join("\t",$gen_begin,$gen_end, $match_len, 0, $real_strand,$mrna_id, $cdna_begin, $cdna_end);
        $min = $gen_begin if  ( $gen_begin < $min );
        $max = $gen_end if  ( $gen_end > $max );
    }
    @a_exon = reverse(@a_exon) if ( $strand_gen == 1 );
    return ($min,$max,@a_exon);
}

