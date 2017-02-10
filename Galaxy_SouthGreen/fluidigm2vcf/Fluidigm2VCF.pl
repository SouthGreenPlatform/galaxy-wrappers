#!/usr/bin/perl

use strict;

use Getopt::Long;
use Bio::SeqIO;

my $usage = qq~Usage:$0 <args> [<opts>]
where <args> are:
    -i, --input         <Fluidigm csv file>
    -s, --snp           <SNP file. Optional>
    -r, --reference     <reference Fasta file. Optional>
    -o, --output        <output filename>
~;
$usage .= "\n";

my ($infile,$reference,$outfile,$snpfile);


GetOptions(
        "input=s"    => \$infile,
        "reference=s"=> \$reference,
        "output=s"   => \$outfile,
        "snpfile=s"  => \$snpfile
);


die $usage
  if ( !$infile || !$outfile);


my %sequences;
if (-e $reference){
	my $in  = Bio::SeqIO->new(-file => $reference ,-format => 'Fasta');
	while ( my $seq = $in->next_seq() ) {
		my $id = $seq -> id();
		#if ($id =~/c(hr\d+)/){$id="C".$1;$id=~s/Chr0/Chr/g;}
		my $sequence = $seq->seq();
		$sequence =~s/\n//g;
		$sequence =~s/\r//g;
		my $nt = substr($sequence,999,1);
		$sequences{$id} = $sequence;
	}
}


my %complement = ("A"=>"T","T"=>"A","G"=>"C","C"=>"G");

#print "$nb\n";
#exit;
my $nb_st=0;
my $nline = 0;
my $ok=0;
my %hash;
my %hash2;
my %individuals;
open(H,$infile);
<H>;
while(<H>)
{
	if (/,Converted,/){$ok=1;next;}
	if ($ok == 0){
		next;
	}
	#my ($snp,$alleles,$chr,$pos,$strand) = split(/\t/,$_);
	my ($id,$snp,$ref_allele,$alt_allele,$ind,$Type,$Auto,$Confidence,$Final,$alleles) = split(/,/,$_);
	if ($ind eq 'H2O'){next;}
	$individuals{$ind} = 1;
	my $pos = "";
	my $chr = 1;
	$alleles =~s/://g;
	$hash{$snp}{$ind} = $alleles;
	$hash2{$snp}{"ref_allele"} = $ref_allele;
	$hash2{$snp}{"alt_allele"} = $alt_allele;
}
close(H);


my $pos = 0;
open(my $O,">$outfile");
print $O "##fileformat=VCFv4.1\n";
print $O "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	".join("\t",keys(%individuals))."\n";
foreach my $snp(keys(%hash)){
	$pos++;
        my $chr = 1;
	my $ref_allele = $hash2{$snp}{"ref_allele"};
	my $alt_allele = $hash2{$snp}{"alt_allele"};
	my $alleles = $ref_allele."/".$alt_allele;
	$alleles =~s/0/A/g;
        $alleles =~s/1/C/g;
        $alleles =~s/2/G/g; 
	my $seq = $sequences{$chr};
	my $indice = $pos-1;
	my $real_pos = $indice + 1;
	if ($seq){
		$ref_allele = uc(substr($seq,$indice,1));
	}

	if (!$ref_allele){
		my @i = split(/\t/,$_);
		my ($a1,$a2) = split(//,$i[12]);
		$ref_allele = $a1;	
	}
	my $strand = "+";
	if ($alleles !~/$ref_allele/)
        {
		$strand = "-";
		$nb_st++;#next;
	}
	#print "$ref_allele $alleles $strand\n";
	#$ref_allele = "G";

	#print "$snp	$alleles	$chr	$real_pos	+	NA	NA	NA	NA	NA	NA";

	my $multi;
	if (!$ref_allele)
	{
		print "Ref: $ref_allele alleles:$alleles\n";
		print "Error:\n";
		next;
	}
	my @al_val = split("/",$alleles);

	


	print $O "$chr	$real_pos	.	$ref_allele	$alt_allele	100	.	AC=8;AF=1.00;AN=8;BaseQRankSum=0.747;DP=60;Dels=0.00;FS=0.000	GT";
		
	foreach my $ind(keys(%individuals))
	{
		my $al = $hash{$snp}{$ind};
		
		if ($al eq "No Call"){
			$al = "..";
		}
		else{
			$al =~s/$ref_allele/0/g;
			$al =~s/$alt_allele/1/g;
		}
		
		print $O "	". join("/",split("",$al));
	}
	print $O "\n";
}	
close(H);
close($O);

