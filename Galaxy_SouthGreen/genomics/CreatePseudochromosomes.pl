#!/usr/bin/perl

use strict;

use Getopt::Long;

my $usage = qq~Usage:$0 <args> [<opts>]
where <args> are:
    -f, --fasta          <FASTA file of scaffolds/contigs>
    -o, --order          <scaffold/contig order and orientation>
    -e, --esp            <species name>
    
      <opts> are:
    -c, --chrom_out      <Chromosome output file. Default: pseudomolecules.fa>
    -a, --agp_out        <AGP output file. Default: chr_from_scaffolds.agp>
    -p, --pseudo_gff_out <GFF output file. Default: pseudomolecules.gff3>
    -s, --scaff_gff_out  <GFF output file. Default: scaffolds.gff3>
    -n, --n_number       <N number for gap. Default: 100>
    -i, --is_chrom       <add "chr" particule for FASTA file(true/false). Default: true>
~;
$usage .= "\n";

my ($scaffolds,$order,$species,$chrom_output,$agp_output,$pseudo_gff,$scaff_gff,$n_number,$add_chrom_particule);

$chrom_output = "pseudomolecules.fa";
$agp_output = "chr_from_scaffolds.agp";
$pseudo_gff = "pseudomolecules.gff3";
$scaff_gff = "scaffolds.gff3";

$n_number = 100;
$add_chrom_particule = "true";

GetOptions(
	"fasta=s"          => \$scaffolds,
	"order=s"          => \$order,
	"esp=s"            => \$species,
	"chrom_output=s"   => \$chrom_output,
	"agp_output=s"     => \$agp_output,
	"pseudo_gff_out=s" => \$pseudo_gff,
	"scaff_gff_out=s"  => \$scaff_gff,
	"n_number=s"       => \$n_number,
	"is_chrom=s"       => \$add_chrom_particule
);


die $usage
  if ( !$scaffolds || !$order || !$species);
  


########################################################
# parse scaffolds FASTA file
######################################################## 
my %scaffolds;
my $scaffold_name;

use Bio::SeqIO;
my $in  = Bio::SeqIO->new(-file => "$scaffolds" , '-format' => 'Fasta');
while ( my $seq = $in->next_seq() ) 
{
	my $id = $seq->id();
	my $seq = $seq->seq();
	$scaffolds{$id} = $seq;
}



########################################################
# parse anchoring results file
########################################################
my %anchors;
open(my $ANCHORING,$order);
<$ANCHORING>;
while(<$ANCHORING>)
{
	my $line = $_;
	$line =~s/\n//g;
	$line =~s/\r//g;
	chomp($line);
	my ($chrom,$ordervalue,$scaffold_num,$sens) = split(/\t/,$line);
	if ($ordervalue =~/\d/)
	{
		my $scaffold_name = $scaffold_num;
		if ($chrom eq "0")
		{
			$chrom = "Un_random";
		}
		$anchors{$chrom}{$ordervalue} = $scaffold_name . "," . $sens;
	}
}
close($ANCHORING);



########################################################
# create pseudochromosome FASTA file
########################################################

open(my $GFF_SCAFF,">$scaff_gff");
open(my $GFF_PSEUDO,">$pseudo_gff");
open(my $AGP,">$agp_output");
print $AGP "##agp-version	2.0\n";
print $GFF_SCAFF "##gff-version 3\n";
print $GFF_PSEUDO "##gff-version 3\n";
my $out  = Bio::SeqIO->new(-file => ">$chrom_output" , '-format' => 'Fasta');
foreach my $chrom(sort(keys(%anchors)))
{
	my $ref = $anchors{$chrom};
	my %hash = %$ref;
	my $scaffold_ok = 0;
	my $pos_start = 1;
	my $pos_end = 1;
	my $n = 1;
	my $chrom_name;
	my $chrom_sequence = "";
	foreach my $num(sort {$a <=> $b} keys(%hash))
	{
		my ($scaffold_name,$sens) = split(",",$anchors{$chrom}{$num});
		my $scaffold_seq = $scaffolds{$scaffold_name};
		
		if (!$scaffold_seq)
		{
			$scaffold_name = "scaffold" . $scaffold_name;
			$scaffold_seq = $scaffolds{$scaffold_name};
		}
		
		
		if ($scaffold_ok)
		{
			$chrom_sequence .= "N"x$n_number;
		}
		$scaffold_ok++;
			
		my $sens_val = "+";
		if ($sens eq "AS" or $sens eq "-")
		{
			$chrom_sequence .= reverseComplement($scaffold_seq);
			$sens_val = "-";
		}
		elsif ($sens eq "S" or $sens eq "nd" or $sens eq "+")
		{
			$chrom_sequence .= $scaffold_seq;
		}
		
		if ($add_chrom_particule eq "true" or $add_chrom_particule eq "True")
		{
			$chrom_name = "chr".$chrom;
		}
		else
		{
			$chrom_name = $chrom;
		}
		
		my $size = length($scaffold_seq);
		$pos_end = $pos_start + $size - 1;
		print $AGP "$chrom_name	$pos_start	$pos_end	$n	W	$scaffold_name	1	$size	$sens_val\n";
		print $GFF_SCAFF "$chrom_name	Assembly	contig	$pos_start	$pos_end	.	$sens_val	.	ID=$scaffold_name;Name=$scaffold_name\n";
		$pos_start = $pos_end + 1;
		
		
		
		$size = $n_number;
		$pos_end = $pos_start + $size - 1;
		print $AGP "$chrom_name	$pos_start	$pos_end	$n	N	$n_number	fragment	no\n";
		$pos_start = $pos_end + 1;
		
		
		$n++;
		
	}
	my $seq = new Bio::Seq(-id=>"$chrom_name",-seq=>$chrom_sequence);
	$out->write_seq($seq);
	
	my $length_chrom = length($chrom_sequence);
	print $GFF_PSEUDO "$chrom_name	$species	chromosome	1	$length_chrom	.	.	.	ID=$chrom_name;Name=$chrom_name\n";
	
}
close($AGP);
close($GFF_PSEUDO);
close($GFF_SCAFF);


sub reverseComplement($)
{
	my $sequence = $_[0];

	my @nts = split("",reverse($sequence));
	my $sequence_reverse = "";
	foreach my $nt(@nts)
	{
		if ($nt eq 'A'){$nt = 'T';}
		elsif ($nt eq 'T'){$nt = 'A';}
		elsif ($nt eq 'G'){$nt = 'C';}
		elsif ($nt eq 'C'){$nt = 'G';}
		$sequence_reverse .= $nt;
	}
	return $sequence_reverse;
}