#!/usr/bin/perl

use FindBin;                 
use lib "$FindBin::Bin/lib";
use Bio::SearchIO;
use Bio::Tools::GFF;
use File::Basename;
use GnpUtil;
my $util = new GnpUtil();
my $file = shift; 
my $tool_data_path = shift;
my $output_gff3 = shift;
my $output_exonerate= shift;


my $database_name = &get_bank($file,$tool_data_path."/blastdb.loc");
&parse_tblastn($file,$database_name ,$output_gff3,$output_exonerate);

sub get_bank {
	my ($file,$blastdb) = @_;
	open(IN,$blastdb); 
	my %bank;
	while(<IN>){
		chomp;
		next if $_ =~ /^#/;
		my ($alias,$path) = (split(/\t/,$_))[0,2];
		chop($path);
		my ($filename,$path,$suffix) = fileparse($path);
		$bank{$filename} = $alias;   
	}
	close IN; 
    my $searchio = new Bio::SearchIO(
		-format => 'blast',
		-file   => $file
	);
    my $result = $searchio->next_result(); 
    my $database_name = (split(/\//,$result->database_name))[-1];  
    my $alias  = $bank{$database_name}; 
	return $alias;
}

sub parse_tblastn {
    my ($file ,$bank ,$output_gff3,$output_exonerate) = @_;
    my $source = join("_","tblastn",$bank);   
    my $searchio = new Bio::SearchIO(
		-format => 'blast',
		-file   => $file
	);	
    my $gff_out = new Bio::Tools::GFF(
		-file => ">$output_gff3",
		-gff_version => 3
	);
    open(TXT,">$output_exonerate");
    while(my $result = $searchio->next_result()) {
		if ($result) {
	    	my $query_length = $result->query_length;
	    	my $query_name   = $result->query_name; 
	    	my $query_desc   = $result->query_description(); 
	    	my $locus_tag    = $query_desc;
	    	if ($locus_tag eq "") {
	    		$locus_tag =  $query_name;
	    	}
	    	my $num_hits     = $result->num_hits();
	    	if ($num_hits >= 0) {
				my $rank_hit = 0;
				while(my $hit = $result->next_hit()) { 
		    		$rank_hit++;
		 		   	last if $rank_hit == 5;
		    		my $hit_length = $hit->length();
		   			my ($lib,$name,$description) = (split(/\|/,$hit->name()));
		   	 		my $est_name = join("_",$lib,$name); 
		    		$description .= " " .$hit->description();
		    		$description = $util->encod($description); 
		    		my $hsp_query_cum_length;	
		    		my $cum_num_identical;
		    		my $hsp_cum_length;
		    		my $hsp_hit_cum_length;
		    		my (@start_query,@end_query,@start_subject,@end_subject);
		    		my $strand;
		    		my $sfr;
		    		while(my $hsp = $hit->next_hsp()){
						$hsp_query_cum_length += $hsp->length( 'query');
						$hsp_hit_cum_length += $hsp->length( 'hit');
						$hsp_cum_length += $hsp->length( 'total');
						$cum_num_identical += $hsp->num_identical();
						push @start_query,$hsp->start('query');
						push @end_query, $hsp->end('query');
						push @start_subject,$hsp->start('subject');
						push @end_subject, $hsp->end('subject');
						$strand = $hsp->strand('hit') == 0 ? "+" : "-";
						$sfr = ($hsp->hit->frame + 1) * $hsp->strand('hit');
		    		}
		   			my $query_start = _min(@start_query);
		    		my $query_end   = _max(@end_query);
		    		my $hit_start   = _min(@start_subject);
		    		my $hit_end     = _max(@end_subject);
		    		my $evalue      = $hit->significance();	
		    		my $qcov = sprintf('%.2f',($hsp_query_cum_length / $query_length));
		    		my $scov = sprintf('%.2f',($hsp_hit_cum_length / $hit_length));
		    		my $identity = sprintf('%.2f',($cum_num_identical / $hsp_cum_length));
		    		$qcov = "1.00" if $qcov > 1;
		    		$scov = "1.00" if $scov > 1;
		    		$identity = "1.00" if $identity > 1;
		    		my ($product,$cds_type) = $util->eval_prediction($qcov,$scov,$identity,$description);
		    		my $id = join("_",$source,	$query_name, $est_name);
		    		my $target = join("+",$est_name, $hit_start,	$hit_end );
		    		my $feature =  new Bio::SeqFeature::Generic(
						-seq_id      => $query_desc,
						-source_tag  => $source,
						-primary_tag => "match",
						-start       => $query_start,
						-end         => $query_end,
						-strand      => $strand,    
						-tag         =>  {
							ID           => $id,
							Name         => $est_name,
							Target       => $target,
							qcov         => $qcov,
							scov         => $scov,
							identity     => $identity,
							completeness => $cds_type, 
							evalue       => $evalue
						}
					);
		    		$feature->add_tag_value("Alias",$alias) if $alias;  
		    		$feature->add_tag_value("description",$product) if $product;                                            
		    		my $evidence_tag= $bank .":" .$est_name;
		    		$gff_out->write_feature($feature);
		    		if ($rank_hit == 1) {
						my $id = $hit->name();
						print TXT join("\t",$locus_tag,$id,$sfr,$product),"\n";
		    		}
				}
	    	}
		}
    }
    $gff_out->close;
    close TXT;	
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

1;
