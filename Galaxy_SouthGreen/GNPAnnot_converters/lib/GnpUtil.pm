package GnpUtil;

use lib '.';
use XML::Simple qw(:strict);
use Getopt::Long;
use Bio::SearchIO;
use Bio::SeqIO;
use Bio::Tools::SeqStats;
use Bio::Index::Fasta;
use Bio::SeqFeature::Generic;
use Config::General;
use Data::Dumper;
use Net::SMTP;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(eval_prediction format_blast format_gth filter_gff format_netgene2 format_splicepredictor FilterReformatRegion FilterReformatA6 format_est format_repeatmasker format_sequence format_sim4 reverse_seq gc_content blast2tab encod decod sort_uniq adjust_coord read_db_conf get_default_conf read_embl_conf format_repseek get_sequence_models get_seq_obj  parse_interpro );

sub new {
    my $this = shift;
    my %arg = @_;
    my $self = {
                %arg
                };
    my $class = ref($this) || $this;
    bless $self, $class;
    return $self;
}

sub format_blast {
    my ($self,$file_in,$type)  = @_;
    my $file_tmp = $file_in ."_tmp";  
    my $file_out = $type eq "sp" ? $self->{target} .".blast0" : $self->{target}.".blast1";
    my $searchio= new Bio::SearchIO(
                                    -format => 'blast',
                                    -file   => $file_in
                                    );
    open(OUT,">".$file_tmp);
    while (my $result = $searchio->next_result) {
        while (my $hit = $result->next_hit) {
            while (my $hsp=$hit->next_hsp) {
                my $strand;
                if ($hsp->strand <0){$strand = "-";}
                else {$strand = "+";}
                my $frame = $hsp->query->frame+1;
                my $sens = $strand . $frame;
                my $score = $hsp->score;
                print OUT join("\t",$hsp->query->start,$hsp->query->end,$hsp->score, $hit->significance, $sens,$hit->name, $hsp->hit->start, $hsp->hit->end), "\n";
            }
        }
    }
    close OUT;
    system("sort -n -k 1,1 $file_tmp | sort -s -k 6,6  > $file_out");
    system("rm $file_tmp");
    return $file_out;
}

sub format_repseek {
    my ($self,$nb,$scaffold) = @_;
    my $file = $self->{output}->{repseek};
    open(IN,"sort -k $nb -n $file|");
    my %result = undef;
    my $cpt = 1;
    my $block = 1;
    my ($min_start_copy,$max_end_copy);
    while(<IN>) {
        chomp;
        my ($type,$first_copy,$second_copy,$length_first_copy,$length_second_copy,$identity,$score) = (split(/\t/,$_))[0,1,2,3,4,7,8];
        my ($type,$strand) = (split(/\./,$type));
        $strand = $strand eq "Dir" ? "+" : "-";
        my ($start_copy,$end_copy,$start_copy1,$end_copy1);
        if ($nb == 2) {
            $start_copy = $first_copy;
            $end_copy = $first_copy + $length_first_copy;
            $start_copy1 = $second_copy;
            $end_copy1 = $second_copy + $length_second_copy;
        }
        else {
            $start_copy = $second_copy;
            $end_copy = $second_copy + $length_second_copy;
            $start_copy1 = $first_copy;
            $end_copy1 = $first_copy + $length_first_copy;
        }
        if ($max_end_copy) {
            if ($start_copy <= $max_end_copy) {
                $id = join("_","repseek",$scaffold,$nb,$block);
                $min_start_copy = $start_copy if $min_start_copy > $start_copy;
                $max_end_copy = $end_copy if  $max_end_copy < $end_copy;
                push @{$result{$id}{start}},$start_copy;
                push @{$result{$id}{end}},$end_copy;
                push @{$result{$id}{identity}},$identity;
                push @{$result{$id}{score}},$score;
                push @{$result{$id}{strand}},$strand;
                push @{$result{$id}{type}},$type;
                push @{$result{$id}{start1}},$start_copy1;
                push @{$result{$id}{end1}},$end_copy1;	
                $result{$id}{min} = $min_start_copy;
                $result{$id}{max} = $max_end_copy;
            }
            else { 
                $block++;
                $id = join("_","repseek",$scaffold,$nb,$block);
                $min_start_copy = $start_copy;
                $max_end_copy = $end_copy;
                push @{$result{$id}{start}},$start_copy;
                push @{$result{$id}{end}},$end_copy;
                push @{$result{$id}{identity}},$identity;
                push @{$result{$id}{score}},$score;
                push @{$result{$id}{strand}},$strand;
                push @{$result{$id}{type}},$type;
                push @{$result{$id}{start1}},$start_copy1;
                push @{$result{$id}{end1}},$end_copy1;	
                $result{$id}{min} = $min_start_copy;
                $result{$id}{max} = $max_end_copy;
            }
        }
        else {
            $min_start_copy = $start_copy; #NFRbeg=
            $max_end_copy = $end_copy;     #NFRend=
            $id = join("_","repseek",$scaffold,$nb,$block);
            push @{$result{$id}{start}},$start_copy;
            push @{$result{$id}{end}},$end_copy;
            push @{$result{$id}{identity}},$identity;
            push @{$result{$id}{score}},$score;
            push @{$result{$id}{strand}},$strand;
            push @{$result{$id}{type}},$type;
            push @{$result{$id}{start1}},$start_copy1;
            push @{$result{$id}{end1}},$end_copy1;	
            $result{$id}{min} = $min_start_copy;
            $result{$id}{max} = $max_end_copy;
        }
        $cpt++;
    }
    close IN;
    return \%result;
}

sub format_repeatmasker {
    my ($self,$file_in) = @_;
    my $file_out = $file_in;
    $file_out =~ s/out/ig/;
    my $ostart = 0;
    my $ostop = 0;
    open(REPEAT,$file_in);
    open(IG,">".$file_out);
    while (<REPEAT>) {
        my @data = (split(/\s+/,$_));
        last if /~They/;
        if ($#data == 14) {
            $start = (split(/\s+/,$_))[5];
            $stop =  (split(/\s+/,$_))[6];
        }
        else {
            $start = (split(/\s+/,$_))[6];
            $stop =  (split(/\s+/,$_))[7];
        }
        if ($start <= $ostop) {
            if ($stop > $ostop) {
                $ostop = $stop;
            }
            next;
        }
        if ($ostop != 0){
            print IG $ostart, "\t", $ostop, "\n";
        }
        $ostart=$start;
        $ostop=$stop;
    }
    print IG $ostart, "\t", $ostop, "\n" if ($ostart != 0 && $ostop !=0);
    close IG;
    close REPEAT;
    return $file_out;
}

sub format_est {
    my ($self, $file_tmp, $tag, $file_out) = @_;
    open(IN, $file_tmp);
    open(OUT,">>".$file_out);
    while(my $line = <IN> ) {
        chop($line);
        my @F = split(' ',$line);
        $F[5] = join("|",$tag,$F[5]);
        my $result  = join(' ',@F);
        print OUT $result ,"\n";
    }
    close IN;
    close OUT;
}

sub format_gth {
    my ($self, $file_tmp, $tag, $file_out) = @_;    
    open(GENESEQER,"$file_tmp") or die "ERROR $0 >$file_tmp<\n";
    open(EST,">>$file_out") or die "ERROR $0 >$file_out<\n";
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
            if ( $mrna_id eq '' ) {
                $mrna_id = $a_F[2];
                $mrna_id =~ s/[\-\+]$//;
            }
            ($strand_cdna) = ( $a_F[2] =~ /\-$/ )  ? 1 : 0;
            ($strand_gen) = ( $a_F[1] =~ /\-$/ )  ? 1 : 0;
        }
        if ( $lign =~ /^Alignment.+:\s*$/ ) {
            if ( ($curr_pci >= $pci) && ($curr_pcs >= $pcs)  && $#a_exon != -1 ) {
                my ($min,$max,@a_exon) = $self->FilterReformatA6($match_len,$strand_gen,$mrna_id, $strand_cdna,$est_len,$tag,@a_exon);
                $h_mrna{$mrna_id} = "$match_len $strand_cdna $strand_gen $min $max";
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
            $seqok = 1;
            next;
        }
        next if ( ! $exonok );
        push @a_exon, $lign if ( $lign =~ /^\s*Exon/ );
    }
    close(GENESEQER);
    close(EST);
}

sub format_sim4 {
    my ($self,$file,$file_slim) = @_;
    open(IN,$file);
    open(OUT,">$file_slim");
    my ($seq1,$seq2);	
    my $cpt = 0;
    my %hash;
    while(my $line = <IN>) {
        next if( $line =~ /^\s+$/); chomp($line);
        if ($line =~ /seq1/) {
            $seq1 = $line;
        }
        else {
            if ($line =~ /seq2/) {
                $seq2 = $line;	
                $cpt = 0;
            }
            else {
                if ($line) {
                    if($line =~ /^\(?(\d+)\-(\d+)\)?\s+\(?(\d+)\-(\d+)\)?\s+(\d+)/ ) {
                        push @{$hash{$seq1}{$seq2}} , $line;
                    }
                    elsif ($line =~ /\s+/) {
                        $cpt++;
                        my $value = ($cpt%4) ;
                        push @{$hash{$seq1}{$seq2}} ,"\n" if $value == 1;
                        push @{$hash{$seq1}{$seq2}} , $line;
                    }
                    else {
                        push @{$hash{$seq1}{$seq2}} , $line;
                    }
                }
            }
        }
    }
    close IN;
    foreach my $seq1 (keys %hash) {
        foreach my $seq2 (keys %{$hash{$seq1}}) {
            print OUT "\n",$seq1 ,"\n";
            print OUT $seq2 ,"\n\n\n";
            print OUT join("\n",@{$hash{$seq1}{$seq2}}) ,"\n";
        }
    }
    close OUT;
    return $file_slim;
}

sub blast2tab {
    my ($self,$file,$start_region) = @_;
    my $in = new Bio::SearchIO(
                                -format => 'blast',
                                -file   => $file
                            );
    my $out_tab = $file .".tab";
    open(OUT,">".$out_tab);
    while( my $result = $in->next_result ) {
        while( my $hit = $result->next_hit ) {
            while( my $hsp = $hit->next_hsp ) {
                my $percent_identity = sprintf("%.2f",$hsp->percent_identity);
                my $strand_query = $hsp->strand('query');
                my $strand_hit   = $hsp->strand('hit');
                my $start_query  = $hsp->start('query') + $start_region - 1;
                my $end_query  = $hsp->end('query') + $start_region - 1;
                my $start_hit  = $hsp->start('hit');
                my $end_hit  = $hsp->end('hit');
                ($start_query,$end_query) = ($end_query,$start_query) if $strand_query == -1;
                ($start_hit,$end_hit) = ($end_hit,$start_hit) if $strand_hit == -1;
                my $homology_string = $hsp->homology_string;
                my $query_string =  $hsp->query_string;
                my $hit_string = $hsp->hit_string;
                my @homology_string = split(//,$homology_string);
                my @query_string = split (//,$query_string);
                my @hit_string  = split(//,$hit_string);
                my $cpt = 0;
                for (my $i=0;$i<=$#homology_string;$i++) {
                    if ($homology_string[$i] eq " " || $homology_string[$i] eq "+") {
                        $cpt++;
                    } 
                }
                my $cptx = 0;
                for (my $i=0;$i<=$#query_string;$i++) {
                    if ($query_string[$i] eq "X") {
                        $cptx++;
                    }
                }
                for (my $i=0;$i<=$#hit_string;$i++) {
                    if ($hit_string[$i] eq "X") {
                        $cptx++;
                    }
                }
                my $mismatch = $cpt - $cptx;
                print OUT join("\t",$result->query_name ,$hit->name,$percent_identity,$hsp->hsp_length,$mismatch,$hsp->gaps,$start_query,$end_query,$start_hit,$end_hit,$hsp->expect,$hsp->bits), "\n";
            }
        }
    }
    close OUT;
    return $out_tab;
}

sub read_db_conf {
    my ($self,$conf) = @_;
    my %conf;
    open CONF, $conf or die "Unable to open $conf : $!\n";
    while (<CONF>) {
        next if /^\#/;
        if (/(\w+)\s*=\s*(\S.*)$/) {
            $conf{$1} = $2;
        }
    }
    close CONF;
    return \%conf;
}

sub adjust_coord {
    my ($self,$hit_name,$query_name,$hit_desc,$hsps) = @_;
    my ($min_start,$max_end,$min_qstart,$max_qend,$gene_id,$cpt_exon,%result,%count_exon,$strand);
    my $cpt_exon;
    my $cpt_gene = 1;
    my $count_hsp = 0;
    my @hsps = @$hsps;
    for (my $i=0;$i<=$#hsps;$i++) {
        if ($min_start) {
            if ($strand != $hsps[$i]->strand('query')) {
                $cpt_gene++;
                $gene_id = $hit_name .".".$cpt_gene;
                $cpt_exon = $count_exon{$gene_id}++;
                $result{$gene_id}[$cpt_exon]{'start_bac'} = $hsps[$i]->start('query');
                $result{$gene_id}[$cpt_exon]{'end_bac'}   = $hsps[$i]->end('query');
                $result{$gene_id}[$cpt_exon]{'start_est'} = $hsps[$i]->start('hit');
                $result{$gene_id}[$cpt_exon]{'end_est'}   = $hsps[$i]->end('hit');
                $result{$gene_id}[0]{'strand'} =  $hsps[$i]->strand('hit');
                $result{$gene_id}[0]{'chr'} =  $query_name;
                $result{$gene_id}[0]{'desc'} =  $hit_desc;
                $min_start = $hsps[$i]->start('hit') if $hsps[$i]->start('hit');
                $max_end = $hsps[$i]->end('hit');
                $min_qstart = $hsps[$i]->start('query');
                $max_qend = $hsps[$i]->end('query');
                $strand = $hsps[$i]->strand('query');
            }
            else {
                if ($hsps[$i]->strand('query') == 1) {
                    if ( $hsps[$i]->start('hit') + $int > $max_end) {
                        $gene_id = $hit_name .".".$cpt_gene;
                        $cpt_exon = $count_exon{$gene_id}++;
                        $result{$gene_id}[$cpt_exon]{'start_bac'} = $hsps[$i]->start('query');
                        $result{$gene_id}[$cpt_exon]{'end_bac'}   = $hsps[$i]->end('query');
                        $result{$gene_id}[$cpt_exon]{'start_est'} = $hsps[$i]->start('hit');
                        $result{$gene_id}[$cpt_exon]{'end_est'}   = $hsps[$i]->end('hit');
                        $min_start = $hsps[$i]->start('hit') if $hsps[$i]->start('hit') < $min_start;
                        $max_end = $hsps[$i]->end('hit') if $hsps[$i]->end('hit') > $max_end;
                        $min_qstart = $hsps[$i]->start('query') if $hsps[$i]->start('query') < $min_qstart;
                        $max_qend = $hsps[$i]->end('query') if $hsps[$i]->end('query') > $max_qend;
                        $strand = $hsps[$i]->strand('query');
                    }
                    else {
                        $cpt_gene++;
                        $gene_id = $hit_name .".".$cpt_gene;
                        $cpt_exon = $count_exon{$gene_id}++;
                        $result{$gene_id}[$cpt_exon]{'start_bac'} = $hsps[$i]->start('query');
                        $result{$gene_id}[$cpt_exon]{'end_bac'}   = $hsps[$i]->end('query');
                        $result{$gene_id}[$cpt_exon]{'start_est'} = $hsps[$i]->start('hit');
                        $result{$gene_id}[$cpt_exon]{'end_est'}   = $hsps[$i]->end('hit');
                        $result{$gene_id}[0]{'strand'} =  $hsps[$i]->strand('hit');
                        $result{$gene_id}[0]{'chr'} =  $query_name;
                        $result{$gene_id}[0]{'desc'} =  $hit_desc;
                        $min_start = $hsps[$i]->start('hit') if $hsps[$i]->start('hit') < $min_start;
                        $max_end = $hsps[$i]->end('hit') if $hsps[$i]->end('hit') > $max_end;
                        $min_qstart = $hsps[$i]->start('query') if $hsps[$i]->start('query') < $min_qstart;
                        $max_qend = $hsps[$i]->end('query') if $hsps[$i]->end('query') > $max_qend;
                        $strand = $hsps[$i]->strand('query');
                    }
                }
                else {
                    if ($hsps[$i]->end('hit') - $int < $min_start) {
                        $gene_id = $hit_name .".".$cpt_gene;
                        $cpt_exon = $count_exon{$gene_id}++;
                        $result{$gene_id}[$cpt_exon]{'start_bac'} = $hsps[$i]->start('query');
                        $result{$gene_id}[$cpt_exon]{'end_bac'}   = $hsps[$i]->end('query');
                        $result{$gene_id}[$cpt_exon]{'start_est'} = $hsps[$i]->start('hit');
                        $result{$gene_id}[$cpt_exon]{'end_est'}   = $hsps[$i]->end('hit');
                        $min_start = $hsps[$i]->start('hit') if $hsps[$i]->start('hit') < $min_start;
                        $max_end = $hsps[$i]->end('hit') if $hsps[$i]->end('hit') > $max_end;
                        $min_qstart = $hsps[$i]->start('query') if $hsps[$i]->start('query') < $min_qstart;
                        $max_qend = $hsps[$i]->end('query') if $hsps[$i]->end('query') > $max_qend;
                        $strand = $hsps[$i]->strand('query');
                    }
                    else {
                        $cpt_gene++;
                        $gene_id = $hit_name .".".$cpt_gene;
                        $cpt_exon = $count_exon{$gene_id}++;
                        $result{$gene_id}[$cpt_exon]{'start_bac'} = $hsps[$i]->start('query');
                        $result{$gene_id}[$cpt_exon]{'end_bac'}   = $hsps[$i]->end('query');
                        $result{$gene_id}[$cpt_exon]{'start_est'} = $hsps[$i]->start('hit');
                        $result{$gene_id}[$cpt_exon]{'end_est'}   = $hsps[$i]->end('hit');
                        $result{$gene_id}[0]{'strand'} =  $hsps[$i]->strand('query');
                        $result{$gene_id}[0]{'chr'} =  $query_name;
                        $result{$gene_id}[0]{'desc'} =  $hit_desc;
                        $min_start = $hsps[$i]->start('hit') if $hsps[$i]->start('hit') < $min_start;
                        $max_end = $hsps[$i]->end('hit') if $hsps[$i]->end('hit') > $max_end;
                        $min_qstart = $hsps[$i]->start('query') if $hsps[$i]->start('query') < $min_qstart;
                        $max_qend = $hsps[$i]->end('query') if $hsps[$i]->end('query') > $max_qend;
                        $strand = $hsps[$i]->strand('query');
                    }
                }
            }
        }
        else {
            $min_start = $hsps[$i]->start('hit');
            $max_end = $hsps[$i]->end('hit');
            $min_qstart = $hsps[$i]->start('query');
            $max_qend = $hsps[$i]->end('query');
            $strand = $hsps[$i]->strand('query');
            $gene_id = $hit_name .".".$cpt_gene;
            $cpt_exon = $count_exon{$gene_id}++;
            $result{$gene_id}[$cpt_exon]{'start_bac'} = $hsps[$i]->start('query');
            $result{$gene_id}[$cpt_exon]{'end_bac'}   = $hsps[$i]->end('query');
            $result{$gene_id}[$cpt_exon]{'start_est'} = $hsps[$i]->start('hit');
            $result{$gene_id}[$cpt_exon]{'end_est'}   = $hsps[$i]->end('hit');
            $result{$gene_id}[0]{'strand'} =  $hsps[$i]->strand('query');
            $result{$gene_id}[0]{'chr'} =  $query_name;
            $result{$gene_id}[0]{'desc'} =  $hit_desc;
        }
        $count_hsp++;
        last if $count_hsp == $nb_hsp;
    }
    return \%result;
}

sub format_sequence {
    my ($self,$sequence) = @_;
    $sequence =~ s/\s+//g;
    my $lenght_line = 60;
    my $length_seq = length($sequence);
    my $nb = int $length_seq / $lenght_line;
    my @formated_sequence;
    if($nb < 1){
        push @formated_sequence, substr($sequence, 0, $lenght_line);
    }
    else{
        for(my $j = 0 ; $j <= $nb ; $j++){
            push @formated_sequence, substr($sequence, $lenght_line * $j, $lenght_line);
        }
    }
    return (join("\n", @formated_sequence) . "\n");
}

sub format_netgene2 {
    my ($self,$file,$length) = @_;
    my $file_out = $file .".gff3";
    my $sequence_name   = $self->{scaffold};   
    open (GFF, ">$file_out") || die "ERROR, can't write >$file_out<\n";
    my $acc_cpt=0;
    my $don_cpt=0;
    my $strand = $file eq $self->{target} .".splices" ? "+" : "-";
	open (GS, $file) || die "ERROR, can't open >$file<\n";
	while (my $line = <GS>) {
		chomp($line);
		my @a_line=split (" ", $line);
		#227 A  0      0      0.006 2  1.043  0.641  0      0      -     - 
		#228 G  0.580  0      0.009 3  1.043  0.641  0.303  0      -     - 
 		my ($start,$end,$score)=(0,0,0);
		if ($#a_line >= 9 ) {
			my ($feature, $feature_number)=("",0) ;
			$a_line[9]= $a_line[3] if ($a_line[9] eq "-");
			$a_line[8]= $a_line[2] if ($a_line[8] eq "-");
			if ( $a_line[9] != 0 ) {	
				$feature ="SO:0000164" ;
				$acc_cpt++;
				$feature_number=$acc_cpt;
				($start,$end) = ($strand eq "+") ? ($a_line[0]-1,$a_line[0]) : ($length-$a_line[0]+1,$length-$a_line[0]+2);
				$score=$a_line[9];
				my $id=$feature;
				$id =~ s/ /_/ig;
				printf GFF "$sequence_name\tNG2\t$feature\t$start\t$end\t$score\t$strand\t.\tID=$id:$sequence_name.$feature_number\n";
			}			
			if ( $a_line[8] != 0 ) {
				$feature  = "SO:0000163";
				$don_cpt++;
				$feature_number=$don_cpt;
				($start,$end) = ($strand eq "+") ? ($a_line[0],$a_line[0]+1) : ($length-$a_line[0],$length-$a_line[0]+1);
				$score=$a_line[8];
				my $id=$feature;
				$id =~ s/ /_/ig;
				printf GFF "$sequence_name\tNG2\t$feature\t$start\t$end\t$score\t$strand\t.\tID=$id:$sequence_name.$feature_number\n";
			} 
		}
	}
	close GS;
    close GFF;
}
sub format_splicepredictor {
    my ($self,$file) = @_;
    my $file_out = $file .".gff3";   
    open (GFF, ">$file_out") || die "ERROR, can't write >$file_out<\n";
    my $acc_cpt=0;
    my $don_cpt=0;
    my $sequence_name = $self->{scaffold}; 
    my $strand = $file eq $self->{target}.".spliceP" ? "+" : "-";
	open (GS, $file) || die "ERROR, can't open >$file<\n";
	while (my $line = <GS>) {
		chomp($line);
		my @a_line=split (" ", $line);
		#A     <-     98 gcagacttgcagtAGat     0.001  0.000  0.000   3 (1 1 1)      IAE-E-EDI
		#D --->      101           gatGTacgt   0.257  0.092  0.000   8 (3 4 1)     IAEE-E-DIIA
 		my ($start,$end,$score) = (0,0,0);
		if ($#a_line >= 11 ) {
			my ($feature, $feature_number) = ("",0) ;
			$score=$a_line[4];
			if ( $a_line[0] eq "A" ) {	
				$feature ="SO:0000164" ;
				$acc_cpt++;
				$feature_number = $acc_cpt;
				($start,$end) = ($strand eq "+") ? ($a_line[2]-1,$a_line[2]) : ($a_line[2],$a_line[2]+1);				
			}			
			if ( $a_line[0] eq "D" ) {
				$feature  = "SO:0000163";
				$don_cpt++;
				$feature_number = $don_cpt;
				($start,$end) = ($strand eq "+") ? ($a_line[2],$a_line[2]+1) : ($a_line[2]-1,$a_line[2]);
			} 
			my $id=$feature;
			$id =~ s/ /_/ig;
            next if $start == 0;
			printf GFF "$sequence_name\tSPred\t$feature\t$start\t$end\t$score\t$strand\t.\tID=$id:$sequence_name.$feature_number\n";
		}
	}
	close GS;
    close GFF;
}

sub gc_content {
    my $self = shift;
    my $seq = Bio::SeqIO->new(
                                -file => $self->{target},
                                -format => "fasta"
                            )->next_seq();
    my $seq_stats  =  Bio::Tools::SeqStats->new(
                                                -seq => $seq
                                            ); 
    my $hash_ref = $seq_stats->count_monomers();
    my $clone_name = $seq->display_id;
    my $length = $seq->length;
    my $gc = ($hash_ref->{'G'} + $hash_ref->{'C'}) / $length * 100;
    my $gc_content = sprintf("%.2f",$gc);
    $self->{clone_name} = $seq->display_id;
    $self->{length}     = $seq->length;
    $self->{gc_content} = $gc;
    foreach my $base (sort keys %{$hash_ref}) {
        $self->{base}->{$base} = $hash_ref->{$base};
    }
}


sub get_default_conf {
    my $self = shift;
    my $conf = $self->{clone_name} .".conf";
    unless (-e $conf) {
        open(CONF,">$conf") or die $!;
        print CONF "
<species>
common_name = $self->{common_name}
classification = $self->{classification}
description = $self->{common_name} genomic DNA, BAC clone $self->{clone_name} , complete sequence
accession  =
</species>

<reference>
title    = $self->{common_name} genomic DNA, BAC clone $self->{clone_name} , complete sequence
location = Unpublished
authors  =
</reference>

<source>
organism = $self->{common_name}
cultivar =
mol_type = genomic DNA
clone    = $self->{clone_name}
db_xref  = taxon:$self->{taxon_id}
</source>
<comment>
text     =
</comment>";
        close CONF;
    }
    $self->{embl_conf} = $conf;
}   

sub read_embl_conf {
    my $self = shift;
    my $conf = new Config::General($self->{embl_conf});
    my %config    = $conf->getall;
    $self->{read_embl} = \%config;
}

sub get_swissprot {
    my ($self,$id, $dbtag) = @_;
    my @id = @$id;
    my $file = $dbtag .".txt";
    open(OUT,">$file");
    print OUT join("\n",@id);
    close OUT; 
    my $dir_prog = $self->{dir_prog};
    my $file_temp = $self->{pwd} ."/".$dbtag .".faa";
    my $gff = $self->{pwd}  ."/gff3/swissprot/".$dbtag.".gff";
    my $db_name = $dbtag eq "SwissProt" ? $self->{db_uniprot}->{sprot}->{fasta} : $self->{db_uniprot}->{trembl}->{fasta};
    my $organism;
    system("/usr/local/bioinfo/meme/bin/fasta-fetch $db_name -f $file > $file_temp");
    print $file,"\n",$file_temp,"\n";
    my $seq_in = new Bio::SeqIO(
			    -format => 'fasta',
			    -file  => $file_temp
			    );
    my $gffio = Bio::Tools::GFF->new(
				     -file => ">$gff",
				     -gff_version => 3
				     );
    while (my $seqobj = $seq_in->next_seq()) {
	my ($id,$alias) = (split(/\|/,$seqobj->display_id))[1,2];
        $target = join("+",$id,$start,$end);
	my $desc = $seqobj->desc();
	my $description;
	my $species;
	my $gene_name;
        if ($desc =~ /(.*)\sOS=(.*)\sGN=(\S+)\s.*/) {
            $description = $1;
            $species     = $2;
	    $gene_name =$3;
        }
        else {
		$description = "Hypothetical protein";
		$species = "Unknown";
	}
        my $feature =  new Bio::SeqFeature::Generic(
                                                    -seq_id      => $id,
                                                    -source_tag  => $dbtag,
                                                    -primary_tag => 'polypeptide',
                                                    -start       => 1,
                                                    -end         => $seqobj->length(),
                                                    -tag         =>
                                                    {
                                                        ID        => $id,
                                                        Dbxref    => "$dbtag:$id"
                                                    }
                                                );
        $feature->add_tag_value("Alias",$alias) if $alias;
        $feature->add_tag_value("organism",$species);
        $feature->add_tag_value("product",$description);
        $gffio->write_feature($feature);

    } 
    $gffio->close; 
    #system("rm $file_temp") if -e $file_temp;
    #system("rm $file") if -e $file;
    return $gff;
}

sub get_sequence_models {
    my ($self,$index_name,$id,$path) = @_;
    my $gff = $path ."/gff3/models/".$id .".gff";
    my $organism;
    my $dbxref;
    if ($index_name eq "sb_pep") {
        $organism = "Sorghum bicolor";
        $dbxref   = "JGI";
    }
    elsif ($index_name eq "os_pep") {
        $organism = "Oryza sativa";
        $dbxref   = "TIGR";
    }
    else {
        $organism = "Oryza sativa";
        $dbxref   = "RAP-DB";
    }
    #unless (-e $gff) {
    my $dir = "/bank/index";
    my $type = "SDBM_File";
    my $dbobj = Bio::Index::Abstract->new("$dir/$index_name");
    my $seq = $dbobj->get_Seq_by_id($id);
    my $length = $seq->length();
    my $dbxref_id = join(":",$dbxref,$id);
    my $description = $seq->description();
    my $gffio = GPI::GFF->new(
                                -gff_version => 3,
                                -file        => ">$gff"
                            );
    my $models_gff3 =  new Bio::SeqFeature::Generic(
                                                    -seq_id      => $id,
                                                    -source_tag  => $dbxref,
                                                    -primary_tag => 'polypeptide',
                                                    -start       => 1,
                                                    -end         => $length,
                                                    -tag         =>  {
                                                                        ID       => $id,
                                                                        Name     => $id,
                                                                        Dbxref   => $dbxref_id,
                                                                        organism => $organism
                                                                    }
                                                    );
    $models_gff3->add_tag_value("Note",$description) if $description;
    print $gffio->write_feature($models_gff3);
    $gffio->close();
    return $gff;
}


sub get_est {
    my ($self,$est_name) = @_;
    my $lib;
    my $name;
    my $id;
    my $species;
    if ($est_name =~ /^CL/ || $est_name =~ /^KZ/) {
        $db_est = "/bank/theobroma_cacao/contig_region.fasta";
        $name = $est_name;
        $id = $est_name;
        $db_name = "TC_contig";
		$species = "Theobroma cacao";
    }
    else {
        ($lib,$name) = (split(/\|/,$est_name))[0,1];
        if ($lib eq "Vitac_vitis") {
            $db_est = "/bank/est/VV_embl_mrnas.fsa";
            $id = join("_",$lib,$name);
            $db_name = "VV_mrnas";
			$species = "Vitis vinifera";
        }
        elsif ($lib eq "singleton") {
            $db_est = "/bank/musa_acuminata/singleton.fasta";
            $id = join("_",$lib,$name);
            $db_name = "MS_mrnas";
			$species = "Musa acuminata";
        }
        elsif ($lib eq "full_length") {
            $db_est = "/bank/musa_acuminata/full_length.fasta";
            $id = join("_",$lib,$name);
            $db_name = "MS_full_length";
			$species = "Musa acuminata";
		}
        elsif ($lib eq "contig") {
            $db_est = "/bank/musa_acuminata/contig.fasta";
            $id = join("_",$lib,$name);
            $db_name = "MS_contig";
			$species = "Musa acuminata";
        }
        elsif ($lib eq "Malva_Cacao") {
            $db_est = "/bank/theobroma_cacao/est_region.fasta";
            $id = join("_",$lib,$name);
            $db_name = "TC_mrnas";
			$species = "Theobroma cacao";
        }
        elsif ($lib eq "Sapin_Citrus") {
            $db_est = "/bank/est/C_embl_mrnas.fsa";
            $id = join("_",$lib,$name);
            $db_name = "C_mrnas";
			$species = "Citrus";
        }
        elsif ($lib eq "TC_454") {
            $db_est = "/bank/theobroma_cacao/consensus_tgicl_454_cacao.fasta";
            $id = $name;
            $db_name = "TC_454";
			$species = "Theobroma cacao";
        }
        elsif ($lib eq "Brass_Arabido") {
            $db_est = "/bank/est/AT_embl_mrnas.fsa";
            $id = join("_",$lib,$name);
            $db_name = "AT_mrnas";
			$species = "Arabidopsis thaliana";
        }
        elsif ($lib eq "Malva_Gossyp") {
            $db_est = "/bank/est/G_embl_mrnas.fsa";
            $id = join("_",$lib,$name);
            $db_name = "G_mrnas";
			$species = "Gossypium";
        }
        elsif ($lib eq "Panico_Saccha") {
            $db_est = "/bank/est/SC_embl_mrnas.fsa";
            $id = join("_",$lib,$name);
            $db_name = "SC_mrnas";
			$species = "Saccharum officinarum";
        }
        elsif ($lib eq "Pooida_Tritic") {
            $db_est = "/bank/est/WH_embl_mrnas.fsa";
            $id = join("_",$lib,$name);
            $db_name = "WH_mrnas";
			$species = "Triticum aestivum";
        }
        elsif ($lib eq "Panico_ZeaMai") {
            $db_est = "/bank/est/ZM_embl_mrnas.fsa";
            $id = join("_",$lib,$name);
            $db_name = "ZM_mrnas";
			$species = "Zea maize";
        }
        elsif ($lib eq "Panico_Sorghu") {
            $db_est = "/bank/est/SB_embl_mrnas.fsa";
            $id = join("_",$lib,$name);
            $db_name = "SB_mrnas";
			$species = "Sorghum bicolor";
        }
        elsif ($lib eq "Pooida_OryzaS") {
            $db_est = "/bank/est/OS_embl_mrnas.fsa";
            $id = join("_",$lib,$name);
            $db_name = "OS_mrnas";
			$species = "Oryza sativa";
        }
        elsif ($lib eq "Cocose_Cocosn" || $lib eq "Phoeni_Phoeni" || $lib eq "Cocose_Elaeis" || $lib eq "Cocose") {
            $db_est = "/bank/est/AR_embl_mrnas.fsa";
            $id = join("_",$lib,$name);
            $db_name = "AR_mrnas";
			$species = "Elaeis guineensis";
        }
        elsif ($lib eq "Pooida_Hordeu") {
            $db_est = "/bank/est/HV_embl_mrnas.fsa";
            $id = join("_",$lib,$name);
            $db_name = "HV_mrnas";
			$species  = "Hordeum vulgare";
        }
        elsif ($lib eq "Populus") {
            $db_est = "/bank/est/Populus_ncbi_mrnas.fsa";
            $id = join("_",$lib,$name);
            $db_name = "PT_mrnas";
			$species = "Populus trichocarpa";
        }
        elsif ($lib eq "Vitis") {
            $db_est = "/bank/est/VV_ncbi_mrnas.fsa";
            $id = join("_",$lib,$name);
            $db_name = "VV_mrnas";
			$species = "Vitis vinifera";
        }
        elsif ($lib =~ /Solana/) {
            $db_est = "/bank/est/Solanales_ncbi_mrnas.fsa";
            $id = join("_",$lib,$name);
            $db_name = "SOL_mrnas";
			if ($lib eq "Solana_NicTab"){$species ="Nicotiana tabacum";}
			if ($lib eq "Solana_SolLyc"){$species ="Solanum lycopersicum";}
			if ($lib eq "Solana_SolTub"){$species ="Solanum tuberosum";}
			if ($lib eq "Solana_CapAnn"){$species ="Capsicum annuum";}
			if ($lib eq "Solana_SolMel"){$species ="Solanum melongena";}
			if ($lib eq "Solana_IpoNil"){$species ="Ipomoea nil";}
			if ($lib eq "Solana_NicBen"){$species ="Nicotiana benthamiana";}
			if ($lib eq "Solana_PetAxi"){$species ="Petunia axillaris subsp. axillaris";}
			if ($lib eq "Solana_SolTor"){$species ="Solanum torvum";}
			if ($lib eq "Solana_SolHab"){$species ="Solanum habrochaites";}
			if ($lib eq "Solana_IpoBat"){$species ="Ipomoea batatas";}
			if ($lib eq "Solana_PetHyb"){$species ="Petunia x hybrida";}
			if ($lib eq "Solana_NicLxS"){$species ="Nicotiana langsdorffii x Nicotiana sanderae";}
			if ($lib eq "Solana_SolPen"){$species ="Solanum pennellii";}
			if ($lib eq "Solana_NicSyl"){$species ="Nicotiana sylvestris";}
			if ($lib eq "Solana_SolCha"){$species ="Solanum chacoense";}
			if ($lib eq "Solana_HyoNig"){$species ="Hyoscyamus niger";}
			if ($lib eq "Solana_SolPhu"){$species ="Solanum phureja";}
			if ($lib eq "Solana_IpoTri"){$species ="Ipomoea trifida";}
			if ($lib eq "Solana_SolLyc"){$species ="Solanum lycopersicum var. cerasiforme";}
		}
		
        elsif ($lib eq "CofCan" || $lib eq "CofAra") {
            $db_est = "/bank/est/Coffea_unigene.fna";
            $id = join("_",$lib,$name);
			$db_name = "CC_unigen";
			$species = $lib eq "CofCan" ? "Coffea canephora":"Coffea arabica";
        }
        elsif ($lib =~ /Lamial/) {
            $db_est = "/bank/est/Lamiales_ncbi_mrnas.fsa";
            $id = join("_",$lib,$name);
            $db_name = "LA_mrnas";
			if ($lib eq "Lamial_MimGut"){$species ="Mimulus guttatus";}
			if ($lib eq "Lamial_StrHer"){$species ="Striga hermonthica";}
			if ($lib eq "Lamial_TriPus"){$species ="Triphysaria pusilla";}
			if ($lib eq "Lamial_TriVer"){$species ="Triphysaria versicolor";}
			if ($lib eq "Lamial_MimGut"){$species ="Mimulus guttatus var. nasutus";}
			if ($lib eq "Lamial_AntMaj"){$species ="Antirrhinum majus";}
			if ($lib eq "Lamial_OciBas"){$species ="Ocimum basilicum";}
			if ($lib eq "Lamial_MimLew"){$species ="Mimulus lewisii";}
			if ($lib eq "Lamial_SalMil"){$species ="Salvia miltiorrhiza";}
			if ($lib eq "Lamial_OleEur"){$species ="Olea europaea";}
			if ($lib eq "Lamial_SesInd"){$species ="Sesamum indicum";}
			if ($lib eq "Lamial_AviMar"){$species ="Avicennia marina";}
			if ($lib eq "Lamial_SalFru"){$species ="Salvia fruticosa";}
			if ($lib eq "Lamial_MenPip"){$species ="Mentha x piperita";}
			if ($lib eq "Lamial_SteRug"){$species ="Stenogyne rugosa";}
			if ($lib eq "Lamial_PicKur"){$species ="Picrorhiza kurrooa";}
			if ($lib eq "Lamial_AcaEbr"){$species ="Acanthus ebracteatus";}
			if ($lib eq "Lamial_AgaRug"){$species ="Agastache rugosa";}
			if ($lib eq "Lamial_TorFou"){$species ="Torenia fournieri";}
			if ($lib eq "Lamial_MenCan"){$species ="Mentha canadensis";}
			
        }
		elsif ($lib eq "Gentia_CofCan" || $lib eq "Gentia_CofAra" || $lib eq "Gentia_CatRos" || $lib eq "Gentia_KadCen" || $lib eq "Gentia_HedTer" || $lib eq "Gentia_OldAff" ||$lib eq "Gentia_EusGra" || $lib eq "Gentia_CofCxC" || $lib eq "Gentia_CofAxC" || $lib eq "Gentia_MitSpe" || $lib eq "Gentia_CofHyb" || $lib eq "Gentia_GymSyl") {
            $db_est = "/bank/est/Gentianales_ncbi_mrnas.fsa";
            $id = join("_",$lib,$name);
            $db_name = "GE_mrnas";
			if ($lib eq "Gentia_CofCan"){$species ="Coffea canephora";}
			if ($lib eq "Gentia_CofAra"){$species = "Coffea arabica";}
			if ( $lib eq "Gentia_CatRos"){$species ="Catharanthus roseus";}
			if ($lib eq "Gentia_KadCen"){$species ="Kadua centranthoides";}
			if ($lib eq "Gentia_HedTer" ){$species ="Hedyotis terminalis";}
			if ($lib eq "Gentia_OldAff") {$species ="Oldenlandia affinis";}
			if ($lib eq "Gentia_EusGra"){$species ="Eustoma grandiflorum";}	
			if ($lib eq "Gentia_CofCxC") {$species ="Coffea canephora x Coffea congensis";}	
			if ($lib eq "Gentia_CofAxC"){$species ="Coffea arabica x Coffea canephora";}	
			if ($lib eq "Gentia_MitSpe"){$species ="Mitragyna speciosa";}	
			if ($lib eq "Gentia_CofHyb"){$species ="Coffea hybrid cultivar";}	
			if ($lib eq "Gentia_GymSyl"){$species ="Gymnema sylvestre";}	
        }
		elsif ($lib eq "Gentia_CofAra_Caturra" || $lib eq "Gentia_CofAra_IAPAR" || $lib eq "Gentia_CofAra_S795" || $lib eq "Gentia_CofEug"){
			$db_est = "/bank/est/coffea_454_mrna.fsa";
            $id = join("_",$lib,$name);
            $db_name = "CC_mrnas";
			$species = "Coffea canephora";
    	}
        else {
            $db_est = "/bank/arabidopsis_thaliana/cDNA_full_reading_071121.txt";
            ($name,$alias) = (split(/\|/,$est_name))[0,1];
            $id = $name;
			$species = "Arabidopsis thaliana";
            $db_name = "AT_full_length";
        }
    }
    my $gff = $self->{pwd} ."/gff3/est/".$id .".gff";
    unless (-e $gff) {
        my $db_est_idx = $db_est.".idx";
        my $dbxref = join(":",$db_name,$id);
        my $tmp_fna = $id."_est.fna";
        my $cmd = "fastafetch -f $db_est -i $db_est_idx -F -q '$est_name' > $tmp_fna";
		my $length=100;
		if (-e $tmp_fna){
			my $seq = Bio::SeqIO->new(
										-file => "$tmp_fna",
										-format => "fasta"
									)->next_seq();
			$length = $seq->length();
        }
		system("rm tmp_fna") if -e "tmp_fna";
        my $gffio = GPI::GFF->new(
								-gff_version => 3,
								-file        => ">$gff"
								);
        my $est_gff3 =  new Bio::SeqFeature::Generic(
													 -seq_id      => $id,
													 -source_tag  => $db_name,
													 -primary_tag => 'region',
													 -start       => 1,
													 -end         => $length,
													 -tag         =>  
													 { 
														 ID       => $id,
														 Dbxref   => $dbxref,
														 organism => $species
													 }
													 );
        print $gffio->write_feature($est_gff3);
        $gffio->close();
        unlink $tmp_fna;
    }
    return ($gff,$db_name,$species);
}


sub get_repeat {
    my ($self,$id,$path) = @_;
    my ($name,$repeat) = (split(/\#/,$id));
    if ($id =~ /gi/) {
        $name = (split(/\|/,$name))[3];
        $name =~ s/\.\d*//;
    }
    my ($repeat_class,$repeat_family) = (split(/\//,$repeat));
    $repeat_class = "Unknown" unless $repeat;
    $repeat_family = "Unknown" unless $repeat_family;
    my $gff = $path ."/gff3/repeat/".$name .".gff";
    unless (-e $gff ){
        my $fna = $path ."/gff3/repeat/".$name .".fna";
        my $gffio = GPI::GFF->new(
                                    -gff_version => 3,
                                    -file        => ">$gff"
                                );
        my $db = $self->{db_repeat_with_repet};
        my $idx = $self->{db_repeat_with_repet} .".idx";
        my $cmd = "fastafetch -f $db -i $idx -F -q '$id' > $fna";
        system($cmd);
        my $seq = Bio::SeqIO->new(
                                    -file => $fna,
                                    -format => "fasta"
                                )->next_seq();
        my $length = $seq->length();
        system("rm $fna");
        my $est_gff3 =  new Bio::SeqFeature::Generic(
                                                        -seq_id      => $id,
                                                        -source_tag  => "repeat",
                                                        -primary_tag => 'region',  
                                                        -start       => 1,
                                                        -end         => $length,
                                                        -tag         =>  {
                                                                            ID       => $id,
                                                                            Name     => $name,
                                                                            repeat_class=>$repeat_class,
                                                                            repeat_subclass=>$repeat_family
                                                                        }
                                                    );
        print $gffio->write_feature($est_gff3);
        $gffio->close();
    }
    return $gff;
}

sub eval_prediction {
	my ($self,$QCov,$SCov,$identity,$description) = @_;
	my $product;
	my $cds_type;
	if ($identity >= 0.5) {
		if ($QCov >= 0.8) {
			if ($SCov >= 0.8) {
				if($description  =~ /[Hh]ypothetical [Pp]rotein/){
					$product = "Conserved $description";
				}
				else{
					$product = sprintf("%s",$description);
				}
				$cds_type = "complete";
			}
			else {
				if($description =~ /[Hh]ypothetical [Pp]rotein/) {
					$product = sprintf("Conserved %s",$description);
				} 
				else {	
					$product = sprintf("%s",$description);
				}
				$cds_type = "fragment";
			}				
		}
		else {
			if ($SCov >= 0.8) {
				if($description =~ /[Hh]ypothetical [Pp]rotein/) {
					$product = sprintf("Conserved %s",$description);
				} 
				else {	
					$product = sprintf("%s",$description);
				}
				$cds_type = "modules";
			}
			else {  
				$product  = "Hypothetical protein";
				$cds_type = "to_fill";
			}		
		}
	}
	elsif ($identity >= 0.25) {
		if ($QCov >= 0.8) {
			if ($SCov >= 0.8) {
				if($description =~ /[Hh]ypothetical [Pp]rotein/){
					$product = "Putative $description";
				}
				elsif($description =~ /[Pp]utative/){
					$product = $description;
				}
				else {
					$product = "Putative $description";
				}
				$cds_type = "complete";
			}
			else {
				if($description =~ /[Hh]ypothetical [Pp]rotein/){
					$product = "Putative $description";
				}
				elsif($description =~ /[Pp]utative/){
					$product = $description;
				}
				else {
					$product = "Putative $description";
				}
				$cds_type = "fragment";
			}
		}
		else {
			if ($SCov >= 0.8) {
				if($description =~ /[Hh]ypothetical [Pp]rotein/) {
					$product = $description;
				} 
				elsif($description =~ /[Pp]utative/){
					$product = $description;
				}   
				else {	
					$product = "Putative $description";
				}
				$cds_type = "modules";
			}
			else {	
				$product  = "Hypothetical protein";
				$cds_type = "to_fill";
			}
		}   
	}
	else {
		$product  = "Hypothetical protein";
		$cds_type = "to_fill";
	}
	return ($product,$cds_type);
}

sub reverse_seq {
    my ($self,$sequence) = @_;
    $sequence = reverse $sequence;
    $sequence =~ tr/gatcGATC/ctagCTAG/;
    return $sequence;
}

sub encod {
    my ($self,$encod) = @_; 
    $encod =~ s/([^a-zA-Z0-9_. :?^*\(\)\[\]@!-])/uc sprintf("%%%02x",ord($1))/eg;
#    $encod =~ s/([,;=%&()'"])/sprintf "%%%02X", ord($1)/gei;
    return $encod;
}
sub decod {
    my ($self,$decod) = @_;
    $decod =~ s/%([0-9a-z]{2})/chr(hex($1))/eig;
    return $decod;
}

sub sort_uniq {
    my ($self,$array) = @_;
    my %seen = ( );
    my @array = grep { ! $seen{$_} ++ }  @{$array};
    return \@array;
}

sub get_seq_obj {
    my ($self,$index_name,$id) = @_;
    my $dir = "/bank/index";
    my $type = "SDBM_File";
    my $dbobj = Bio::Index::Abstract->new("$dir/$index_name");
    my $seq = $dbobj->get_Seq_by_id($id);
    return $seq;
}

sub FilterReformatA6 {
    my ($self,$match_len,$strand_gen,$mrna_id,$strand_cdna,$est_len,$tag,@a_exon) = @_;  
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
        $min = $gen_begin if  ( $gen_begin < $min );
        $max = $gen_end if  ( $gen_end > $max );
    }
    @a_exon = reverse(@a_exon) if ( $strand_gen == 1 );
    return ($min,$max,@a_exon);
}

sub FilterReformatRegion {
	my ($self,$match_len,$strand_gen,$mrna_id,$strand_cdna,$est_len,$start_region ,@a_exon) = @_;
	my $min=99999999999999999999999;
	my $max=-1;
	my $min_cdna=99999999999999999999999;
	my $max_cdna=-1;
	my ($real_strand)  = $strand_gen ;
	for(my $i=0;$i<=$#a_exon;$i++) {
		my ($gen_begin,$gen_end) = $a_exon[$i] =~ /Exon\s+\d+\s+(\d+)\s+(\d+)/;
		# strand = 1 (revcomp), 0 (direct)
		#print "seqlen $match_len mrna_len $est_len\n";
		my ($cdna_begin,$cdna_end) = $a_exon[$i] =~ /cDNA\s+(\d+)\s+(\d+)/;
		($gen_begin,$gen_end)   = ( $strand_gen == 0 ) ? ($gen_begin,$gen_end) : ($gen_end,$gen_begin);
		$gen_begin = $gen_begin + $start_region - 1;
		$gen_end   = $gen_end   + $start_region - 1;
		($cdna_begin,$cdna_end)= (($est_len-$cdna_end+1),($est_len-$cdna_begin+1)) if ( $strand_gen == 1 );
		($cdna_begin , $cdna_end) = ( $cdna_begin > $cdna_end ) ? ($cdna_end,$cdna_begin) : ($cdna_begin , $cdna_end); # A=6 format
		$a_exon[$i] = join("\t",$gen_begin,$gen_end,$real_strand,$mrna_id,$cdna_begin,$cdna_end);
		$min = $gen_begin if  ( $gen_begin < $min );
		$max = $gen_end if  ( $gen_end > $max );
		$min_cdna = $cdna_begin if  ( $cdna_begin < $min_cdna);
		$max_cdna = $cdna_end if  ( $cdna_end > $max_cdna );
	}
	@a_exon = reverse(@a_exon) if ( $strand_gen == 1 );
	return ($min,$max,$min_cdna,$max_cdna,@a_exon);
}

sub filter_gff {
    my ($self,$file) = @_;
    my $file_tmp = $file .".tmp";
    my $gffi = GPI::GFF->new(
                            -file => $file,
                            -gff_version => 3
                            );
    my $gffo = GPI::GFF->new(
                            -file => ">$file_tmp",
                            -gff_version => 3
                            );
    my $found = 0;
    my $name;
    my %est;
    my %est_part;		
    my $id_parent;
    my $id;		  
    while($feature = $gffi->next_feature()) {
        my $source = $feature->source_tag();
        my $primary_tag = $feature->primary_tag();
        ($name)  = $feature->get_tag_values("Name") if $feature->has_tag("Name");
        if ($primary_tag eq "EST_match" || $primary_tag eq "cDNA_match" || $primary_tag eq "protein_match") {
            my $start = $feature->start();
            my $end   = $feature->end();
            if ($est{$name}{$start}) {
                $found = 1;
            }
            else {
                $est{$name}{$start} = $feature;
                $found = 0;
                ($id)  = $feature->get_tag_values("ID") if $feature->has_tag("ID");
            }
        }
        else {
            ($id_parent)  = $feature->get_tag_values("Parent") if $feature->has_tag("Parent");
            if ($id_parent eq $id) {
                push @{$est_part{$id_parent}}, $feature;
            }
        }
    }
    $gffi->close();
    foreach my $est (keys %est){
        foreach my $start(keys %{$est{$est}}) {
            my $feature = $est{$est}{$start};
            print $gffo->write_feature($feature);
            my ($id)  = $feature->get_tag_values("ID") if $feature->has_tag("ID");
            foreach my $id_parent (keys %est_part) {
                if ($id eq $id_parent) {
                    foreach my $feature_part( @{$est_part{$id_parent}}) {
                        print $gffo->write_feature($feature_part);
                    }
                }
            }
        }
    }
    $gffo->close();
    system("mv $file_tmp $file");
    return $file;
}


sub parse_interpro {
    my ($self,$gnpdb,$file_xml,$organism_id,$cvterm_id,$cvterm_desc_id,$analysis_id,$cvterm_note_id,$part_of) = @_;
    if (-e $file_xml) {
	my $ref = XMLin(
			$file_xml, 
			forcearray=>1, 
			keeproot=>1, 
			forcecontent=>1, 
			keyattr=>{}
			);
	my $source = "InterPro";
	my $prot_ref = $ref->{interpro_matches}->[0]->{protein};
	my $count;
	foreach my $protein_ref(@$prot_ref) {
	    $count++;
	    my $polypeptide_id = $protein_ref->{id};
	    $polypeptide_id =~ s/_g/_p/;
	    my $length = $protein_ref->{length};
	    my ($srcfeature_id) = $gnpdb->feature($polypeptide_id);
	    if ($srcfeature_id) {
		my ($protein_ref) = $protein_ref->{interpro};
		foreach my $interpro_ref (@$protein_ref) {
		    my $ipr_id    =  $interpro_ref->{id};
		    next if $ipr_id eq "noIPR";
		    my $ipr_desc  =  $interpro_ref->{name};
		    my $type      =  $interpro_ref->{type};
		    my $min_start = 0;
		    my $max_end   = 0;
		    my ($dbxref_id,$ipr_version) = $gnpdb->dbxref($ipr_id);
		    my $uniquename = join("_",$polypeptide_id,$ipr_id);
		    if ($dbxref_id) {
			($feature_ipr_id) = $gnpdb->feature($uniquename);
			unless ($feature_ipr_id) {
			    $gnpdb->insert_feature_with_dbxref($dbxref_id,$uniquename,$organism_id,$cvterm_id);
			    ($feature_ipr_id) = $gnpdb->feature($uniquename);
			}
		    }
		    else {
			$ipr_version = 1;
			$gnpdb->insert_dbxref("InterPro",$ipr_id,$ipr_version,$name);
			my ($dbxref_id,$ipr_version) = $gnpdb->dbxref($ipr_id);
			$gnpdb->insert_feature_with_dbxref($dbxref_id,$uniquename,$organism_id,$cvterm_id);
			($feature_ipr_id) = $gnpdb->feature($uniquename);
		    }
		    my $match_ref = $interpro_ref->{match};
		    foreach my $match (@$match_ref) {
			my $db_id        = $match->{id};
			my $dbname       = $match->{dbname};
			my $db_desc      = $match->{name};
			my $score        = $match->{score};
			my $location_ref = $match->{location};
			my $number_of_match = scalar(@$location_ref);
			my ($dbxref_db_id,$db_version) = $gnpdb->dbxref($db_id);
			my $uniquename_db = join("_",$polypeptide_id,$db_id);
			my $rank = 0;
			my ($start,$end,$score);
			my @match_feature_id;
			if ($number_of_match == 1) {
			    my $feature_db_id;
			    if ($dbxref_db_id) {
				($feature_db_id) = $gnpdb->feature($uniquename_db);
				unless ($feature_db_id) {
				    $gnpdb->insert_feature_with_dbxref($dbxref_db_id,$uniquename_db,$organism_id,$cvterm_id);
				    ($feature_db_id) = $gnpdb->feature($uniquename_db);
				}
			    }
			    else {
				$db_version = 1;
				$gnpdb->insert_dbxref($dbname,$db_id,$db_version,$db_desc);
				my ($dbxref_db_id,$db_version) = $gnpdb->dbxref($db_id);
				$gnpdb->insert_feature_with_dbxref($dbxref_db_id,$uniquename_db,$organism_id,$cvterm_id);
				($feature_db_id) = $gnpdb->feature($uniquename_db);
			    }
			    foreach my $location (@$location_ref) {
				my $temp_start    = $location->{start} - 1;
				my $temp_end      = $location->{end};
				$gnpdb->insert_featureloc($feature_db_id,$srcfeature_id,$temp_start,$temp_end,"0",$rank);
				$rank++;
				if ($start) {
				    $start = $temp_start if $temp_start < $min_start;
				    $end   = $temp_end  if $temp_end > $max_end;
				    $score = $location->{score};
				}
				else {
				    $start = $temp_start;
				    $end   = $temp_end;
				    $score = $location->{score};
				}
			    }
			    $gnpdb->insert_featureprop($feature_db_id,$cvterm_desc_id,$db_desc,'1');
			    if ($score eq "NA") {$score = "0";}
			    $gnpdb->insert_analysisfeature($feature_db_id,$analysis_id,$score);
			    if ($min_start) {
				$min_start = $start if $start < $min_start;
				$max_end   = $end   if $end > $max_end;
			    }
			    else {
				$min_start = $start ;
				$max_end = $end;
			    }
			}
			else {
			    foreach my $location (@$location_ref) {
				my $temp_start    = $location->{start} - 1;
				my $temp_end      = $location->{end};
				my $uniquename4match = $uniquename_db . ".".$rank;
				($feature_db_id_match) = $gnpdb->feature($uniquename4match);
				unless ($feature_db_id_match) {
				    $gnpdb->insert_feature_with_dbxref($dbxref_db_id,$uniquename4match,$organism_id,$cvterm_id);
				    ($feature_db_id_match) = $gnpdb->feature($uniquename4match);
				    push @match_feature_id , $feature_db_id_match;
				}
				$gnpdb->insert_featureloc($feature_db_id_match,$srcfeature_id,$temp_start,$temp_end,"0",$rank);
				$rank++;
				if ($start) {
				    $start = $temp_start if $temp_start < $min_start;
				    $end   = $temp_end  if $temp_end > $max_end;
				    $score = $location->{score};
				}
				else {
				    $start = $temp_start;
				    $end   = $temp_end;
				    $score = $location->{score};
				}
			    }
			    my $feature_db_id;
			    if ($dbxref_db_id) {
				($feature_db_id) = $gnpdb->feature($uniquename_db);
				unless ($feature_db_id) {
				    $gnpdb->insert_feature_with_dbxref($dbxref_db_id,$uniquename_db,$organism_id,$cvterm_id);
				    ($feature_db_id) = $gnpdb->feature($uniquename_db);
				}
			    }
			    else {
				$db_version = 1;
				$gnpdb->insert_dbxref($dbname,$db_id,$db_version,$db_desc);
				my ($dbxref_db_id,$db_version) = $gnpdb->dbxref($db_id);
				$gnpdb->insert_feature_with_dbxref($dbxref_db_id,$uniquename_db,$organism_id,$cvterm_id);
				($feature_db_id) = $gnpdb->feature($uniquename_db);
			    }
			    $gnpdb->insert_featureloc($feature_db_id,$srcfeature_id,$temp_start,$temp_end,"0",$rank);
			    $gnpdb->insert_featureprop($feature_db_id,$cvterm_desc_id,$db_desc,'1');
			    if ($score eq "NA") {$score = "0";}
			    $gnpdb->insert_analysisfeature($feature_db_id,$analysis_id,$score);
			    foreach my $feature_match_id (@match_feature_id) {
				$gnpdb->insert_feature_relationship($feature_match_id,$feature_db_id,$part_of,'0');
			    }
			    if ($min_start) {
				$min_start = $start if $start < $min_start;
				$max_end   = $end   if $end > $max_end;
			    }
			    else {
				$min_start = $start ;
				$max_end = $end;
			    }		
			    
			}
		    }
		    if ($ipr_id) {
			$gnpdb->insert_featureloc($feature_ipr_id,$srcfeature_id,$min_start,$max_end,"0","0");
			$gnpdb->insert_featureprop($feature_ipr_id,$cvterm_desc_id,$ipr_desc,'1');
			$gnpdb->insert_featureprop($feature_ipr_id,$cvterm_note_id,$type,'1');
			$gnpdb->insert_analysisfeature($feature_ipr_id,$analysis_id,"0");
		    }
		}
	    }
	}
    }
}


1;


