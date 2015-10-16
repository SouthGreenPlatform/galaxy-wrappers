#!/usr/bin/perl
use FindBin;
use Bio::SearchIO;
use Bio::Tools::GFF;
use Data::Dumper; 
use lib "$FindBin::Bin/lib";  
use GnpUtil;
use GPI::GFF;
my $file    	= shift; 
my $output_gff3 = shift;
my $output_embl = shift;
my $util    = new GnpUtil();
 

my $db_name = &get_bank($file);
&parse_exonerate($file,$db_name,$output_gff3,$output_embl);

sub get_bank {
	my $file = shift; 
	open(IN,$file);
	my $command = <IN>; 
	my $prog = "cdna2genome";
	if ($command =~ /protein2genome/){
		$prog = "protein2genome";
	} 
	close IN; 
	my $db_name; 
	my $in = Bio::SearchIO->new(
		-file => $file,
		-format => 'exonerate'
	);
	my $result = $in->next_result();
	my $query_name = $result->query_name(); 
	if ($prog eq "protein2genome") {
		if ($query_name =~ /^Os\d+/){
			$bank = "ORYSJ";
		}
		elsif ($query_name =~ /^tr\|/){
			$bank = "TrEMBL";
		}
		elsif ($query_name =~ /^sp\|/){
			$bank = "SwissProt";
		}
		elsif ($query_name =~ /^AT\d+/){
			$bank = "ARATH";
		} 	
		elsif ($query_name =~ /^Solyc/){
			$bank = "SOLLC";
		}
		elsif ($query_name =~ /^Tc\d+/) {
			$bank = "THECC";
		}
		elsif ($query_name =~ /^Bradi/) {
			$bank = "BRADI";
		}
		elsif ($query_name =~ /^Sobic/  ||  $query_name =~ /^Sb\d+/){ 
			$bank = "SORBI";
		}
		elsif ($query_name =~ /\_FGP/ 	||  $query_name =~ /^GRM/) {
			$bank = "MAIZE";
		}
		elsif ($query_name =~ /^\d/){
			$bank = "PDB";
		}
		else {
			$bank = "NR";
		}
	}
	else {
	
	}
	$in->close;
	return $bank;
}
sub parse_exonerate {
    my ($file,$db_name,$output_gff3,$output_embl) = @_; 
    open(EMBL,">$output_embl");
    my $out = new GPI::GFF(
		-file => ">$output_gff3",
		-gff_version=> 3
	);
    my $source = join("_","exonerate",$db_name);  
    my %exist;
    my $in = Bio::SearchIO->new(
		-file => $file,
		-format => 'exonerate'
	) or die "Can't open $file\n";
    while(my $result = $in->next_result()) {
		my $color = $db_name eq "sp" ? 10 : 7;
		my $query_name = $result->query_name();
		my $query_description = $result->query_description(); 
		my $count_hit = 0; 
		while(my $hit = $result->next_hit()) {
			my $query_name = $result->query_name();
			if ($hit) {
				my $hit_name  = $hit->name;
				my ($query , $chr ,$start_region , $end_region) = (split(/\-/,$hit_name));
				my $scaffold =$chr;
				my ($lib,$name,$est_name,$match,$alias);
				if ($db_name =~ /mrnas/  || $db_name =~ /\_contig/){
		 			($lib,$name) = (split(/\|/,$query_name))[0,1];
		    		$est_name = join("_",$lib,$name);
		    		$match = "gene";
		    		$description = $name;
				}
				elsif ($db_name eq "OG_ngs") {
		    		$name = $query_name;
		  		  	$est_name = $name;
		    		$match = "gene";
		    		$description = $name;
				}
				elsif ($db_name eq "SORBI") {
		    		$name = $query_name;
		    		$description = join(" ",$name,$query_description);
		    		$match = "gene";
				}
				elsif ($db_name eq "ARATH" ) {
					$match = "gene";
					$name = $query_name;
					$description = (split(/\|/,$query_description))[2]; 
					$description =~ s/^\s(.*)\s$/$1/;
				}
				elsif ($db_name eq "SOLLC") {
					$match = "gene";
					$name = $query_name;
					$description = $query_description ;  
				}
				elsif ($db_name eq "ORYSJ") {
		    		$name = (split(/\|/,$query_name))[0];
		    		$match = "gene";
		    		$description = $query_description;
				}
				else {
		    		($name,$alias) = (split(/\|/,$query_name))[1,2];
		    		$description = join(" ",$name,$query_description);
		    		$match = "gene";
				}
				$description = $util->encod($description);
				my @hsps = sort { $a-> start('query') <=> $b->start('query') } $hit->hsps();
				if (@hsps) {  
		    		my $hit_length = $hit->length();
		   	 		my ($coverage,$identity,$conserved);
		    		for my $hsp (@hsps) {
						$coverage += $hsp->length('subject');
						$conserved += $hsp->num_conserved();
						$identity += $hsp->num_identical();
		    		}
		    		my $strand = $hsps[0]->strand('hit') == 1 ? "+" : "-"; 
		    		#$count_hit++;
		    		my ($new_start,$new_end) ;
		    		if ($strand eq "-") {
						$new_start = $hsps[$#hsps]->start('hit') + $start_region - 1 ;
						$new_end = $hsps[0]->end('hit') + $start_region - 1;
		    		} 
		    		else {
						$new_start = $hsps[0]->start('hit') + $start_region - 1 ;
						$new_end = $hsps[$#hsps]->end('hit') + $start_region - 1;
		    		}
		    		my $target_start = $hsps[0]->start('query');
		    		my $target_end = $hsps[$#hsps]->end('query'); 
		    		my $name1 = join("_",$query,$name);
		    		$count_hit = 1 + $exist{$name1}++; 
		    		my $gene_id = join(	 "_", $query,  $name, "gene",  sprintf('%04d', $count_hit) );
		    		my $mrna_id = join(	 "_", $query,  $name, "mrna",  sprintf('%04d', $count_hit) );
		    		my $target = join("+",$name,$target_start,$target_end);
		    		my $feature_gene =  new Bio::SeqFeature::Generic(
						-seq_id      => $scaffold,
						-source_tag  => $source,
						-primary_tag => "gene",
						-start       => $new_start,
						-end         => $new_end,
						-strand      => $strand,
						-tag         =>  {
							ID		=> $gene_id,	
							Name	=> $name,
							Target	=> $target 
						}
					);
					$feature_gene->add_tag_value("Note",$description) if $description;
					$feature_gene->add_tag_value("Alias",$alias) if $alias; 
					$out->write_feature($feature_gene);
		    		my $feature_mrna =  new Bio::SeqFeature::Generic(
						-seq_id      => $scaffold,
						-source_tag  => $source,
						-primary_tag => "mRNA",
						-start       => $new_start,
						-end         => $new_end,
						-strand      => $strand,
						-tag         =>  {
							ID		=> $mrna_id,
							Parent	=> $gene_id,	
							Name	=> $name,
							Target	=> $target 
						}
					);
					$feature_mrna->add_tag_value("Note",$description) if $description;
					$feature_mrna->add_tag_value("Alias",$alias) if $alias; 
					$out->write_feature($feature_mrna);
					my $count_hsp = 0;
					my @location;
					for my $hsp ( @hsps ) {
						$count_hsp++;
						my $id_hsp = join("_",$source,$query,$name,"match_part".sprintf('%04d', $count_hsp));	
						my $start_query = $hsp->start('hit') + $start_region - 1;
						my $end_query = $hsp->end('hit') + $start_region - 1;
						my $start_hit  = $hsp->start('query');
						my $end_hit = $hsp->end('query');
						my $target_part = join("+",$name,$start_hit,$end_hit);
						my $match_part = new Bio::SeqFeature::Generic(
							-seq_id      => $scaffold,
							-source_tag  => $source,	
							-primary_tag => 'exon',
							-start       => $start_query,
							-end         => $end_query,	
							-strand      => $strand,
							-tag         =>  {
								Parent => $mrna_id,
								Target => $target_part
							}
						);
						push @location,$start_query ."..".$end_query;
						$out->write_feature($match_part);
					}
					my $tag;
					if ($strand eq "+") {
						if ($#location==0) {
							$tag = join(',',@location);
						}
						else {
							$tag = "join(".join(',',@location).")";
						}
					}
					else {
						if ($#location==0) {
							$tag = "complement(".join(',',@location).")"; 
						}
						else {
							$tag = "complement(join(".join(',',@location)."))";
						}
					}
		    
					print EMBL "FT   BLASTCDS        $tag\n";
					print EMBL "FT                   /note=\"$description\"\n" if $description;
					print EMBL "FT                   /target_id=\"$name\"\n";
					print EMBL "FT                   /target_start=\"$target_start\"\n";
					print EMBL "FT                   /target_end=\"$target_end\"\n";
					print EMBL "FT                   /color=$color\n";
				}
			}
		}
    }
    close EMBL;
    $out->close; 
}
