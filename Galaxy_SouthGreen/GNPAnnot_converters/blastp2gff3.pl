#!/usr/bin/perl

use FindBin;                 
use lib "$FindBin::Bin/lib";
use Bio::SearchIO;
use File::Basename;
use Data::Dumper;
use GnpUtil qw(encod eval_prediction);
use Bio::Tools::GFF;
my $util = new GnpUtil(); 
my $file1              = shift;
my $file2              = shift;
my $file3              = shift;
my $file_gff3          = shift;
my $file_embl          = shift;
my $output1_gff3       = shift;
my $output2_gff3 	   = shift;
my $output3_gff3 	   = shift;
my $output_eugene_gff3 = shift;
my $output_eugene_embl = shift;	
my $product;
my $gene; 
($output1_gff3, $product,$gene) = &parse_blastp($file1,$product,$gene,$output1_gff3);
($output2_gff3,$product,$gene) = &parse_blastp($file2,$product,$gene,$output2_gff3);
($output3_gff3, $product,$gene) = &parse_blastp($file3,$product,$gene,$output3_gff3);

my $gff_in = new Bio::Tools::GFF(
	-file => $file_gff3,
	-gff_version => 3
);
my $out    = new Bio::Tools::GFF (
	-file => ">$output_eugene_gff3",
	-gff_version => 3
);

my $embl_in    = Bio::SeqIO->new(
	-file   => $file_embl,
	-format => 'EMBL' 
);
my $embl_out = new Bio::SeqIO(
	-file => ">$output_eugene_embl",
	-format => "EMBL"
);
my $seqobj = $embl_in->next_seq();
my %locus_info;
my %locus_gene;
foreach my $locus (sort keys %$gene) {
    foreach my $source (keys %{$gene->{$locus}}) {
		my @gene = @{$gene->{$locus}->{$source}};
		$locus_gene{$locus} = {
			gene    => \@gene
	    };  
    }
}
foreach my $locus (sort keys %$product) {
    my $function = "Hypothetical protein";
    my $cds      = "missing_functional_completeness";
    my $origin   = "N/A";
    my @gene;
	foreach my $source (keys %{$product->{$locus}}) {
     	my @match = @{$product->{$locus}->{$source}};
      	foreach my $match (@match) {
	  		if ($match->{product} ne "Hypothetical protein") {
	      		$function = $match->{product};
	      		$cds     = $match->{cds_type};
	      		$origin  = $source;
	      		last ;
	  		}
      	}
  	}
    $locus_info{$locus} = {
		product => $function,
		cds_type => $cds
    };  
}

while(my $feature = $gff_in->next_feature) {
    my ($id) = $feature->get_tag_values("ID") if $feature->has_tag("ID");
    if (defined $locus_info{$id}) {
		$feature->add_tag_value("Ontology_term","PRODUCT:".$locus_info{$id}{product});
		$feature->add_tag_value("Ontology_term","CC_functional_completeness:".$locus_info{$id}{cds_type}); 
		my $gene_name = "unknown_gene";
		if ((scalar(@{$locus_gene{$id}{gene}}))) {
	    	$gene_name = $locus_gene{$id}{gene}[0];
	    	foreach my $gene (@{$locus_gene{$id}{gene}}) {
				$gene = uc($gene);
				$feature->add_tag_value("Ontology_term","CC_gene:$gene") if $gene;
	    	}
		}
		my $note = join("~ ",$id,$locus_info{$id}{product},$gene_name,$locus_info{$id}{cds_type});
		$feature->add_tag_value("note",$note);
		$out->write_feature($feature);
	}
    else {
		$out->write_feature($feature);
    }
}
$out->close;

foreach my $feature ($seqobj->get_SeqFeatures) {
    $feature->remove_tag("ID") if $feature->has_tag("ID");
    $feature->remove_tag("Name") if $feature->has_tag("Name");
    $feature->remove_tag("Parent") if $feature->has_tag("Parent");
    my ($locus_tag) = $feature->get_tag_values("locus_tag") if $feature->has_tag("locus_tag");
    if ($feature->primary_tag() eq "gene") {
		$locus_tag =~ s/_g/_p/;
		my $gene_name = "unknown_gene";
		if (defined $locus_info{$locus_tag}) {
	    	$feature->add_tag_value("product",$locus_info{$locus_tag}{product});
		}
    }
    if ($feature->primary_tag() eq "CDS") {
		$locus_tag =~ s/_g/_p/;
		my $id=$locus_tag;
		if (defined $locus_info{$locus_tag}) {
	    	my $gene_name = "unknown_gene";
	    	if ((scalar(@{$locus_gene{$id}{gene}}))) {
				$gene_name = $locus_gene{$id}{gene}[0];
				foreach my $gene (@{$locus_gene{$id}{gene}}) {
		    		$gene = uc($gene);
			    	$feature->add_tag_value("Ontology_term","CC_gene:$gene") if $gene;
				}
	    	}
	    	my $note = join("~ ",$id,$locus_info{$id}{product},$gene_name,$locus_info{$id}{cds_type});
	    	$feature->add_tag_value("product",$locus_info{$locus_tag}{product});
	    	$feature->add_tag_value("note",$note);
		}
    }
}
$embl_out->write_seq($seqobj);
$embl_out->close;
	
sub parse_blastp {
    my ($file,$product,$gene_name,$file_gff) = @_;
    my $searchio = new Bio::SearchIO( 
		-format => 'blast',
		-file   => $file
	); 
    my $gff_out = new Bio::Tools::GFF(
		-file => ">$file_gff",
		-gff_version => 3
	);
    my %product = %$product;
    my %gene_name = %$gene_name;
    #my $bank;
    while(my $result = $searchio->next_result()) { 
		my $bank;  
    	if ($result) {
	    	my $num_hits = $result->num_hits();
	    	if ($num_hits >= 0) {
				my $query_length = $result->query_length;
				my $query_name   = $result->query_name; 
				my $query_desc   = $result->query_description; 
				my $locus_tag = $query_desc;
				my $rank_hit = 0;
				while(my $hit = $result->next_hit()) { 
		    		$rank_hit++;
		    		last if $rank_hit == 5;
		    		my $hit_length = $hit->length();
		   			my $hit_name = $hit->name();
		   			my $hit_name_split = (split(/\|/,$hit_name))[0];
		   			if ($hit_name_split eq "tr"){
		   				$bank = "TrEMBL";
		   			}
		   			elsif ($hit_name_split eq "sp"){
		   				$bank = "SwissProt";
		   			}
		   			elsif ($hit_name_split =~ /^AT\d+/){
		   				$bank = "ARATH";
		   			}
		   			elsif ($hit_name_split =~ /^Os\d+/){
		   				$bank = "ORYSJ";
		   			}
		   			elsif ($hit_name_split =~ /^Sobic/  ||  $hit_name_split =~ /^Sb\d+/){ 
		   				$bank = "SORBI";
		   			}
		   			else {
		   				$bank = "NR";
		   			}
					my $source = join("_","blastp",$bank); 	
    				my $have_product = join("_","have_product",$bank);
		    		my ($description,$species,$gene_name,$dbxref,$name,$alias);	
		    		if ($hit->description() =~ /(.*)\sOS=(.*)\sGN=(\S+)\s.*/) {
						$description = $1;
						$species = $2;
						$gene_name = $3;
						push @{$gene_name{$locus_tag}{$source}} , $gene_name ;
		    		}
		    		elsif ($hit->description() =~ /(.*)\sOS=(.*)\sPE=/) {
						$description = $1;
						$species = $2;
		    		} 
		    		else {
						if ($bank eq "ARATH") {
			    			$description = (split(/\|/,$hit->description()))[2];
						}   
						else {
			    			$description = $hit->description();
						}
						$species = "Arabidopsis thaliana";
		    		}
		    		$description = "N/A" if $description eq "";
		    		my $description = $util->encod($description);
		    		if ($bank eq "SwissProt" || $bank eq "TrEMBL") {
						($name,$alias) =(split(/\|/,$hit_name))[1,2];
						if ($bank eq "SwissProt") {
			    			$dbxref = "SwissProt:".$name;
						} 
						else {
			    			$dbxref = "TrEMBL:".$name;
						}
		    		}
		    		elsif ($bank eq "ORYSJ") {
						($name,$alias) =(split(/\|/,$hit_name))[0,1];
						my $tigr_id = $name;
						$tigr_id =~ s/\.\d//g;
						$dbxref = "MSU:".$name;
						$species = "Oryza sativa";
		   			}	
		    		elsif ($bank eq "SORBI") {
						$name = $hit_name;
						my $jgi_id = $name;
						$jgi_id =~ s/\.\d//g;
						push @{$self->{jgi_id}} , $jgi_id;
						$dbxref = "JGI:".$name;
						$species = "Sorghum bicolor";
		    		}
		    		else {
						$name = $hit_name;
		    		}
		    		my $hsp_query_cum_length;	
		    		my $cum_num_identical;
		    		my $hsp_cum_length;
		    		my $hsp_hit_cum_length;
		    		my (@start_query,@end_query,@start_subject,@end_subject);
		    		my $strand;
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
		    		}
		    		my $query_start = _min(@start_query);
		    		my $query_end   = _max(@end_query);
		    		my $hit_start   = _min(@start_subject);
		    		my $hit_end     = _max(@end_subject);
		    		my $evalue      = $hit->significance();
		    		my $qcov = sprintf('%.2f',($hsp_query_cum_length / $query_length));
		    		my $scov = sprintf('%.2f',($hsp_hit_cum_length / $hit_length));
		    		my $identity = sprintf('%.2f',($cum_num_identical / $hsp_cum_length));
		    		my ($product,$cds_type) = $util->eval_prediction($qcov,$scov,$identity,$description);
		    		my $id = join("_",$source,$query_name,$name);
		    		my $target = join("+",$name,$hit_start,$hit_end);
		    		my $evidence_tag = $bank .":" .$name;
		    		my $rename = join("_",$name ,"renamed");
		    		my $feature =  new Bio::SeqFeature::Generic(
						-seq_id      => $locus_tag,
						-source_tag  => $source,
						-primary_tag => 'protein_match',
						-start       => $query_start,
						-end         => $query_end,
						-strand      => $strand,
						-tag         =>  {  
							ID          => $id,
							Name        => $rename,
							Target      => $target,
							qcov        => $qcov,
							scov        => $scov,
							identity    => $identity,
							description => $product,
							completeness=> $cds_type,
							organism    => $species,
							evalue      => $evalue
						}
					);
		    		$gff_out->write_feature($feature);
					#	push @{$product{$locus_tag}{$source}} , {gene=>$gene_name} ;
		    		push @{$product{$locus_tag}{$source}} , {
						product => $product,
						cds_type => $cds_type
					};    
				}
		    }
		}
    }
    $gff_out->close;
    return ($file_gff,\%product,\%gene_name);
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
