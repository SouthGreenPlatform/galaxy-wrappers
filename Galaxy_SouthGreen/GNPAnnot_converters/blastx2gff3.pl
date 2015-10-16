#!/usr/bin/perl
use FindBin;                 
use lib "$FindBin::Bin/lib";
use Bio::SearchIO;
use File::Basename;

use Data::Dumper;
use GnpUtil qw(encod eval_prediction);
use Bio::Tools::GFF;
my $util = new GnpUtil();
my $file = shift;
my $tool_data_path   = shift; 
my $output_gff3      = shift;
my $output_embl      = shift;
my $output_exonerate = shift;

my $database_name = &get_bank($file,$tool_data_path."/blastdb_p.loc");

&parse_blastx($file,$database_name,$output_gff3,$output_embl,$output_exonerate);
sub get_bank {
	my ($file,$blastdb) = @_;
	open(IN,$blastdb);
	print $blastdb,"\n";
	my %bank;
	while(<IN>){
		chomp;
		next if $_ =~ /^#/;
		my ($alias,$path) = (split(/\t/,$_))[0,2];
		$path =~ s/\s*//g;
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
sub parse_blastx {
    my ($file,$bank,$output_gff3,$output_embl,$output_exonerate) = @_;
 
    my $in = new Bio::SearchIO(
		-file => $file,
		-format => 'blast'
	);
    
    my $gff_out = new Bio::Tools::GFF->new(
		-file => ">$output_gff3",
		-gff_version => 3
	); 
    open(EMBL,">$output_embl");
    open(TXT,">$output_exonerate");
    my %exist;
    while(my $result = $in->next_result()) {     
		my $source_tag = join("_","blastx",$bank);
		my $blast_type = "BLASTX_HIT";
		my $color = 7;
		if ($result) {
			my $cpt = 0;
			my $cpt_gene = 0;
			my $int = 10;
			my $count_hit = 0;
			my $query = $result->query_name();
			my ($query_name , $chr , $start_region , $end_region) = (split(/\-/,$query));
			my $scaffold = $chr;
			my %result;
			while(my $hit = $result->next_hit()) {
				my $hit_name = $hit->name();
				my $hit_desc = $hit->description();
				my $num_hsps = $hit->num_hsps();       
				my @hsps = sort {$a->start('query') <=> $b->start('query')} $hit->hsps;
				my ($min_start,$max_end,$min_qstart,$max_qend,$gene_id,$cpt_exon,%result,%count_exon,$strand);
				my $cpt_exon;
				my $cpt_gene = 1;
				my $count_hsp = 0;
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
							$result{$gene_id}[$cpt_exon]{'frame'}   = $hsps[$i]->frame('query');
							$result{$gene_id}[0]{'strand'} =  $hsps[$i]->strand('hit');
							$result{$gene_id}[0]{'chr'} =  $chr;
							$result{$gene_id}[0]{'query_name'} = $query_name;
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
									$result{$gene_id}[$cpt_exon]{'frame'}   = $hsps[$i]->frame('query');
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
									$result{$gene_id}[$cpt_exon]{'frame'}   = $hsps[$i]->frame('query');
									$result{$gene_id}[0]{'strand'} =  $hsps[$i]->strand('hit');
									$result{$gene_id}[0]{'chr'} =  $chr;
									$result{$gene_id}[0]{'query_name'} = $query_name;
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
									$result{$gene_id}[$cpt_exon]{'frame'}   = $hsps[$i]->frame('query');
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
									$result{$gene_id}[$cpt_exon]{'start_bac'} 	= $hsps[$i]->start('query');
									$result{$gene_id}[$cpt_exon]{'end_bac'}   	= $hsps[$i]->end('query');
									$result{$gene_id}[$cpt_exon]{'start_est'} 	= $hsps[$i]->start('hit');
									$result{$gene_id}[$cpt_exon]{'end_est'}   = $hsps[$i]->end('hit');
									$result{$gene_id}[$cpt_exon]{'frame'}   = $hsps[$i]->frame('query');
									$result{$gene_id}[0]{'strand'} =  $hsps[$i]->strand('query');
									$result{$gene_id}[0]{'chr'} =  $chr;
									$result{$gene_id}[0]{'query_name'} = $query_name;
									$result{$gene_id}[0]{'desc'} 				=  $hit_desc;
									$min_start 	= $hsps[$i]->start('hit') if $hsps[$i]->start('hit') < $min_start;
									$max_end 	= $hsps[$i]->end('hit') if $hsps[$i]->end('hit') > $max_end;
									$min_qstart = $hsps[$i]->start('query') if $hsps[$i]->start('query') < $min_qstart;
									$max_qend 	= $hsps[$i]->end('query') if $hsps[$i]->end('query') > $max_qend;
									$strand 	= $hsps[$i]->strand('query');
								}
							}
						}
					}
					else {
						$min_start	= $hsps[$i]->start('hit');
						$max_end 	= $hsps[$i]->end('hit');
						$min_qstart = $hsps[$i]->start('query');
						$max_qend 	= $hsps[$i]->end('query');
						$strand 	= $hsps[$i]->strand('query');
						$gene_id 	= $hit_name .".".$cpt_gene;
						$cpt_exon 	= $count_exon{$gene_id}++;
						$result{$gene_id}[$cpt_exon]{'start_bac'} = $hsps[$i]->start('query');
						$result{$gene_id}[$cpt_exon]{'end_bac'}   = $hsps[$i]->end('query');
						$result{$gene_id}[$cpt_exon]{'start_est'} = $hsps[$i]->start('hit');
						$result{$gene_id}[$cpt_exon]{'end_est'}   = $hsps[$i]->end('hit');
						$result{$gene_id}[$cpt_exon]{'frame'}   = $hsps[$i]->frame('query');
						$result{$gene_id}[0]{'strand'}			  =  $hsps[$i]->strand('query');
						$result{$gene_id}[0]{'chr'} 			  =  $chr;
						$result{$gene_id}[0]{'query_name'} 		  = $query_name;
						$result{$gene_id}[0]{'desc'} 			  =  $hit_desc;
					}
					$count_hsp++;
				}
				$cpt_gene++;  
				foreach my $gene (sort keys %result) {
					my @location;
					my @exon = @{$result{$gene}};
					my $strand = $result{$gene}[0]{'strand'} == 1 ? "+" : "-";
					my $chr = $result{$gene}[0]{'chr'};
					my $desc = $util->encod($result{$gene}[0]{'desc'});
					my ($min_start,$max_end,$max_qend,$min_qstart);
					for (my $j=0;$j<=$#exon;$j++) {
						if ($min_qstart) {
							$min_start = $exon[$j]{'start_bac'} if $exon[$j]{'start_bac'} < $min_start ;
							$max_end   = $exon[$j]{'end_bac'} if  $exon[$j]{'end_bac'} > $max_end  ;
							$min_qstart = $exon[$j]{'start_est'} if $exon[$j]{'start_est'} < $min_qstart;
							$max_qend   = $exon[$j]{'end_est'} if $exon[$j]{'end_est'} > $max_qend;
						}
						else {
							$min_start = $exon[$j]{'start_bac'};
							$max_end   = $exon[$j]{'end_bac'};
							$min_qstart = $exon[$j]{'start_est'};
							$max_qend	= $exon[$j]{'end_est'};
						}
					}
					my ($name,$alias);
                    if ($bank eq  "SwissProt" || $bank eq "TrEMBL") { 
                        ($name,$alias) =(split(/\|/,$hit_name))[1,2]; 
                    } 
                    elsif ($bank eq "ORYSJ") {
                        ($name,$alias) =(split(/\|/,$hit_name))[0,1];
                    }
                    elsif ($bank eq "ARATH") {
                         $name = $hit_name ;
                         $desc = (split(/\|/,$result{$gene}[0]{'desc'}))[2];
                         $desc =~ s/^\s(.*)\s$/$1/;
                    } 
                    elsif ($bank eq "SOLLC") {
                         $name = $hit_name ;
                         $desc =  $result{$gene}[0]{'desc'};
                    } 
                    elsif ($bank eq "SORBI") {
                        $name = $hit_name;
                    }
                    else {
                        $name = $hit_name;
                    }
					my $target = join("+",$name,$min_qstart,$max_qend);
					my $gene_id = join("_","gene",$gene,$query_name);
					my $mrna_id = join("_","mrna",$gene,$query_name);
					$min_start = $min_start + $start_region - 1;
					$max_end = $max_end + $start_region - 1;
					my $feature_gene  =  new Bio::SeqFeature::Generic(
						-seq_id      => $scaffold,
						-source_tag  => $source_tag,
						-primary_tag => 'gene',
						-start       => $min_start,
						-end         => $max_end,
						-strand      => $strand,
						-tag         => {
											ID 		=> $gene_id,
											Name 	=> $name,	
											Target	=> $target
								   		}
							       );
		    		$feature_gene->add_tag_value("Note",$desc) if $desc;
					$feature_gene->add_tag_value("Alias",$alias) if $alias;
					$gff_out->write_feature($feature_gene);
					my $feature_mrna  =  new Bio::SeqFeature::Generic(
						-seq_id      => $scaffold,
						-source_tag  => $source_tag,
						-primary_tag => 'mRNA',
						-start       => $min_start,
						-end         => $max_end,
						-strand      => $strand,
						-tag         => {
							ID 		=> $mrna_id,
							Parent	=> $gene_id,
							Name 	=> $name,	
							Target	=> $target
						}
					);
		    		$feature_mrna->add_tag_value("Note",$desc) if $desc;
					$feature_mrna->add_tag_value("Alias",$alias) if $alias;
					$gff_out->write_feature($feature_mrna);
					for (my $j = 0;$j<=$#exon;$j++) {
						my $cpt_exon = $j+1;
						my $target = join("+",$name,$exon[$j]{'start_est'},$exon[$j]{'end_est'});
						my $start_exon = $exon[$j]{'start_bac'}+ $start_region - 1;
						my $end_exon =$exon[$j]{'end_bac'} + $start_region - 1;
						push @location,$start_exon ."..".$end_exon;
						my $match_part  =  new Bio::SeqFeature::Generic(
							-seq_id      => $scaffold,
							-source_tag  => $source_tag,
							-primary_tag => 'exon',
							-start       => $start_exon,
							-end         => $end_exon,
							-frame  	=> $exon[$j]{'frame'},
							-strand      => $strand,
							-tag         =>  	{
								Parent => $mrna_id,
								Target => $target
							}
						);
						$gff_out->write_feature($match_part);
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
					print EMBL "FT   BLASTX_HIT      $tag\n";
					print EMBL "FT                   /note=\"$desc\"\n" if $desc;
					print EMBL "FT                   /color=$color\n";
					if ($count_hit == 1) {
						$gene =~ s/\.\d*//g;
						if ($bank eq "ORYSJ" || $bank eq "SwissProt" || $bank eq "TrEMBL") {
							if (defined $exist{$query}{$name}) {
								next;
							}
							else {
								$exist{$query}{$name} = 1;
								print TXT join("\t",$query,$name,$desc),"\n";
							}
						} 
						else {
							if (defined $exist{$query}{$hit_name}) {
								next;
							}
							else {
								$exist{$query}{$hit_name} = 1;
								print TXT join("\t",$query,$hit_name,$desc),"\n";
							} 
						}
					}
				}
				$count_hit++; 
			}
		} 
    }
    close EMBL;
    close TXT;
    $gff_out->close(); 
}
