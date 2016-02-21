#!/usr/bin/perl

use strict;
use Switch;
use Getopt::Long;
use Bio::SeqIO;

my $usage = qq~Usage:$0 <args> [<opts>]
where <args> are:
    -v, --vcf           <VCF input>
    -g, --gff           <GFF annotation>
    -o, --out           <output in VCF>
    -t, --type          <type: reference sequence correspond to CDS, mRNA>
~;
$usage .= "\n";

my ($vcf,$gff,$out,$type);


GetOptions(
	"vcf=s"      => \$vcf,
	"gff=s"      => \$gff,
	"out=s"      => \$out,
    "type=s"     => \$type
);


die $usage
  if ( !$vcf || !$gff || !$out || !$type);
  
my %correspondence;
my %exons;
my %chrom_of_gene;
my %strands;
open(my $GFF,$gff);
while(<$GFF>)
{
	my @infos = split(/\t/,$_);
	my $feature = $infos[2];
	my $name;
	if ($feature eq 'mRNA' && /ID=([^;]*);.*Name=([^;]*)/)
	{
		$correspondence{$1} = $2;
		$name = $2;
	}
	if ($feature eq $type && /Parent=(.*)/)
	{
		my $chrom = $infos[0];
		my $begin = $infos[3];
		my $end = $infos[4];
		my $strand = $infos[6];
		my $gene = $correspondence{$1};
		if (!$gene && $type eq 'mRNA'){$gene = $name;}
		$chrom_of_gene{$gene} = $chrom;
		$exons{$gene}{$begin}{"end"} = $end;
    	$strands{$gene} = $strand;
	}
}
close($GFF);


my %hsps;
foreach my $g(keys(%exons))
{
	my $ref_hash = $exons{$g};
	my %hash = %$ref_hash;
	my @sorted_exons = sort {$a <=> $b } keys(%hash);
	my $sens = $strands{$g};
    if ($sens eq '-')
    {
    	@sorted_exons = sort {$b <=> $a } keys(%hash);
    }
	my $chromosome_position = 0;
	my $cumulated_size = 0;
	my $begin_query = 1;
    my $end_query;
    foreach my $start(@sorted_exons)
    {
    	if (!$chromosome_position)
		{
			$chromosome_position = $start;
		}
    	my $end = $exons{$g}{"$start"}{"end"};
    	my $size = $end - $start + 1;
    	$cumulated_size += $size;
    	$end_query = $begin_query + $size - 1;
    	if ($sens eq '-')
    	{
    		$hsps{$g}{"$end-$start"}{"pos_query"} = "$begin_query-$end_query";
    		$hsps{$g}{"$end-$start"}{"gaps_on_query"} = "";
    		$hsps{$g}{"$end-$start"}{"gaps_on_subject"} = "";
    	}
    	elsif ($sens eq '+')
    	{
    		$hsps{$g}{"$start-$end"}{"pos_query"} = "$begin_query-$end_query";
	    	$hsps{$g}{"$start-$end"}{"gaps_on_query"} = "";
	    	$hsps{$g}{"$start-$end"}{"gaps_on_subject"} = "";
    	}
    			
    	$begin_query = $end_query + 1;
    }
	#$snp_mapping{$g}{"chromosome_position"} = $chromosome_position;
}

my %pairs = ("A" => "T","T" => "A","G" => "C","C" => "G");

open(my $OUT,">$out");
open(my $VCF,$vcf);
while(<$VCF>)
{
	my @infos = split(/\t/,$_);
	if (scalar @infos > 8 && !/#CHROM/)
	{
		my $gene = $infos[0];
		my $position = $infos[1];
		my $id = $infos[2];
		my $ref_allele = $infos[3];
		my $alt_allele = $infos[4];
		
		my $ref_hsps = $hsps{$gene};
		my %hsps2 = %$ref_hsps;
		
		my $HSP;
		my $start_hsp;
		foreach my $hsp(keys(%hsps2))
		{
			if ($hsps{$gene}{$hsp}{"pos_query"})
			{
				my @pos_query = split("-",$hsps{$gene}{$hsp}{"pos_query"});
				if (($position <= $pos_query[0] && $position >= $pos_query[1]) or ($position <= $pos_query[1] && $position >= $pos_query[0]))
				{
					$HSP = $hsp;
					$start_hsp = $pos_query[0];
				}
			}
		}
		my $strand = $strands{$gene};
		my $hsp_pos = 0;
		if ($HSP =~/^(\d+)-(\d+)$/)
		{
			$hsp_pos = $1;
		}
		my $chromosome_position;
		if ($strand eq "-")
		{
			$chromosome_position = $hsp_pos - ($position - $start_hsp);
			$ref_allele = $pairs{$ref_allele};
			$alt_allele = $pairs{$alt_allele};
		}
		elsif ($strand eq "+")
		{
			$chromosome_position = $hsp_pos + ($position - $start_hsp);
		}
		my $chrom = $chrom_of_gene{$gene};
		if (/$gene	$position	$id	\w+	\w+	(.*)$/)
		{
			my $end_of_line = $1;
			$end_of_line =~s/;MQ0=/;GenePos=$position;MQ0=/g;
			print $OUT "$chrom	$chromosome_position	$id	$ref_allele	$alt_allele	$end_of_line\n";
		}
	}
	else
	{
		print $OUT $_;
	}
}
close($OUT);
  

