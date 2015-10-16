#!/usr/bin/perl

use Bio::Tools::GFF;
use Bio::Tools::Fgenesh;
 
my $file_fgenesh  = shift;
my $output_gff3   = shift;
my $output_embl   = shift;
my $output_eugene = shift;
my $seqname;
# Reformat fgenesh for Eugene
&reformat_fgenesh($file_fgenesh,$output_eugene);

# Reformat fgenesh on GFF3 & EMBL
&fgenesh2gff($file_fgenesh,$output_gff3,$output_embl);

sub reformat_fgenesh {
    my ($file,$output_eugene) = @_;   
    open(IN,$file);
    open(OUTPUT,">".$output_eugene);
    while (my $line = <IN>) {
        my @a_line = split (/\s+/, $line);
        if ($#a_line == 12  &&  $a_line[4]=~ /^CDS/) {
            my $seqname = $a_line[1];
            my $strand  = $a_line[2];
            my $feature = $a_line[4];
            my $score   = $a_line[8];
            my $start   = $a_line[5];
            my $end     = $a_line[7];
            $feature = "E.Init" if $feature eq "CDSf";
            $feature = "E.Intr" if $feature eq "CDSi";
            $feature = "E.Term" if $feature eq "CDSl";
            $feature = "E.Sngl" if $feature eq "CDSo";
            print OUTPUT join("\t",$seqname,"fgenesh",$feature,$start,$end,$score,$strand,"."),"\n";
        }
    }
    close OUTPUT;
    close INPUT; 
}



sub fgenesh2gff {
    my ($file_fgenesh,$output_gff3,$output_embl) = @_; 
    my $source = "fgenesh";
    my $gffio = new Bio::Tools::GFF->new(
		-file => ">$output_gff3",
		-gff_version => 3
	);
    open(EMBL,">$output_embl");
    my $fgenesh =  Bio::Tools::Fgenesh->new(-file => $file_fgenesh);
    my $cpt = 1;
    while (my $prediction = $fgenesh->next_prediction()) {
		my $scaffold = $prediction->seq_id();
		my $transcript = sprintf ($source ."_"."%s_t%06d", $scaffold, $cpt*10);
		my $gene_id = sprintf ($source ."_"."%s_g%06d", $scaffold, $cpt*10);
		my $polypeptide_id = sprintf ($source ."_"."%s_p%06d", $scaffold, $cpt*10);
		$prediction->seq_id($scaffold);
		$prediction->add_tag_value('ID',$gene_id);
		$prediction->primary_tag('gene');
		$prediction->source_tag($source);
		$prediction->add_tag_value("Name",$gene_id);
		$gffio->write_feature($prediction);
		$prediction->remove_tag('ID');
		$prediction->remove_tag('Name');
		$prediction->add_tag_value('ID',$transcript);
		$prediction->primary_tag('mRNA');
		$prediction->add_tag_value("Name",$transcript);
		$prediction->add_tag_value("Parent",$gene_id);
		$gffio->write_feature($prediction);
		$prediction->remove_tag('ID');
		$prediction->add_tag_value('ID',$polypeptide_id);
		$prediction->primary_tag('polypeptide');
		$prediction->add_tag_value("Derives_from",$transcript);
		$gffio->write_feature($prediction);
		my @exons = $prediction->exons();
		my ($location,@location);
		foreach my $exon (@exons) {
	    	if($exon->is_coding()) {
				$exon->primary_tag('CDS');
	    	}
	    	$exon->seq_id($scaffold);
	    	$exon->source_tag('fgenesh');
	    	$exon->primary_tag('CDS');
	    	$exon->add_tag_value('Parent',$transcript);
	    	$exon->remove_tag("score");
	    	$gffio->write_feature($exon);
	    	$location = $exon->start ."..".$exon->end;
	    	push @location,$location;
		}
		my $strand = $prediction->strand() == 1 ? "+" : "-";
		if ($strand eq "+") {
	    	if ($#location==0) {
				print EMBL "FT   CDS             ".join(',',@location) ,"\n";
	    	}
	    	else {
				print EMBL "FT   CDS             join(".join(',',@location).")\n";
	    	}
		}
		else {
	    	if ($#location==0) {
				print EMBL "FT   CDS             complement(".join(',',@location).")\n";
	    	}
	    	else {
				print EMBL "FT   CDS             complement(join(".join(',',@location)."))\n";
	    	}
		}
		print EMBL "FT                   /locus_tag=\"$gene_id\"\n";
		$cpt++;
    }
    close EMBL;
    $gffio->close();
    $fgenesh->close(); 
}
