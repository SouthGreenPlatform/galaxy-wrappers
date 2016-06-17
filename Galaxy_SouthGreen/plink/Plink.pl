
#!/usr/bin/perl

use strict;
use Getopt::Long;
use Bio::SeqIO;

my $usage = qq~Usage:$0 <args> [<opts>]

where <args> are:

    -i, --input          <VCF input>
    -o, --out            <Output basename>
      
      <opts> are:

    -s, --samples        <Samples to be analyzed. Comma separated list>
    -c, --chromosomes    <List of chromosomes to be analyzed.>
    -e, --export         <Output format (VCF/freq/plink. Default: VCF>
    -f, --frequency      <Minimum MAF. Default: 0.001>
    -m, --max_freq       <Maximum MAF. Default: 0.5>
    -a, --allow_missing  <Allowed missing data proportion per site. Must be comprised between 0 and 1. Default: 1>
    -n, --nb_alleles     <Accepted number of alleles (min,max). Default: 2,4>
    -t, --type           <Type of polymorphisms to keep (ALL/SNP/INDEL). Default: ALL>
    -b, --bounds         <Lower bound and upper bound for a range of sites to be processed (start,end). Default: 1, 100000000>
    -r, --remove_filt    <Remove all sites with a FILTER flag other than PASS (true/false). Default: false>
    -d, --distance       <Thin sites so that no two sites are within the specified distance from one another. Default: 0>
~;
$usage .= "\n";

my ($input,$out);

my $PLINK_EXE = "/usr/local/bioinfo/plink/1.90b3k/plink";


#my $indel_size_max = 500;
#my $indel_size_min = 1;
my $frequency_max = 0.5;
my $frequency_min = 0.001;
my $pos_max = 100000000000;
my $pos_min = 0;
my $filter_snp_type = "all";
my $remove_filt = "False";

my $missing_data = 1;
my $export = "VCF";
my $type = "ALL";
my $nb_alleles;
my $bounds;
my $samples;
my $chromosomes;
my $thin;

GetOptions(
	"input=s"        => \$input,
	"out=s"          => \$out,
	"samples=s"      => \$samples,
	"chromosomes=s"  => \$chromosomes,
	"frequency=s"    => \$frequency_min,
	"max_freq=s"     => \$frequency_max,
	"allow_missing=s"=> \$missing_data,
	"export=s"       => \$export,
	"type=s"         => \$type,
	"nb_alleles=s"   => \$nb_alleles,
	"bounds=s"       => \$bounds,
	"remove_filt=s"  => \$remove_filt,
	"distance=s"         => \$thin
);


die $usage
  if ( !$input || !$out);

if ($samples && $samples =~/^([\w\,\-\.]+)\s*$/){
	$samples = $1;
}
elsif ($samples){
	die "Error: Samples must be a comma separated list of string\n";
}
if ($bounds && $bounds =~/^([\d\,]+)\s*$/){
	$bounds = $1;
}
elsif($bounds){
	die "Error: Bounds must be a comma separated list of integers\n";
}

my $minfreq_cmd = "";
if ($frequency_min && $frequency_min > 0 && $frequency_min =~/^([\d\.]+)\s*$/){
	$frequency_min = $1;
	$minfreq_cmd = "--maf $frequency_min";
}
elsif ($frequency_min){
	die "Error: frequency must be an integer\n";
}
if ($thin && $thin =~/^([\d\.]+)\s*$/){
        $thin = $1;
}
elsif ($thin){
        die "Error: frequency must be an integer\n";
}
my $maxfreq_cmd = "";
if ($frequency_max && $frequency_max =~/^([\d\.]+)\s*$/){
	$frequency_max = $1;
	if ($frequency_max < 0.5){
		$maxfreq_cmd = "--max-maf $frequency_max";		
	}
}
elsif($frequency_max){
	die "Error: frequency must be an integer\n";
}
if ($missing_data =~/^([\d\.]+)\s*$/){
	$missing_data = $1;
	$missing_data = 1 - $missing_data;	
}
elsif ($missing_data){
	die "Error: Missing data must be an integer\n";
}
if ($nb_alleles && $nb_alleles =~/^(\d+,\d+)\s*$/){
	$nb_alleles = $1;
}
elsif($nb_alleles){
	die "Error: Nb alleles must be an integer\n";
}
if ($export && $export =~/^([\w]+)\s*$/){
	$export = $1;
}
elsif($export){
	die "Error: Export must be a string\n";
}
if ($type && $type =~/^([\w]+)\s*$/){
	$type = $1;
}
elsif($type){
	die "Error: Type must be a string\n";
}


my @dnasamples;
if ($samples)
{
	@dnasamples = split(",",$samples);
}
my @nalleles;
if ($nb_alleles)
{
	@nalleles = split(",",$nb_alleles);
}
my @boundaries;
if ($bounds)
{
	@boundaries = split(",",$bounds);
}


my $experiment = "chromosomes";
my $table = "";
my %genes;
my @snp_ids;
my @snp_ids_and_positions;
my @snp_ids_and_positions_all;
my $gene;
my $snp_num = 0;
my %ref_sequences;
my %snps_of_gene;

my $indiv_cmd = "";
if (@dnasamples)
{
	if (scalar @dnasamples > 1)
	{
		open(my $S,">$out.samples");
		foreach my $samp(@dnasamples){
			print $S "$samp	$samp\n";
		}
		close($S);
		$indiv_cmd = "--keep $out.samples ";
	}
	else
	{
		$indiv_cmd = "--indv " . join(" --indv ",@dnasamples);
	}
}

my $chrom_cmd = "";
if ($chromosomes)
{
	$chrom_cmd = "--chr ".$chromosomes
}

my $export_cmd = "--recode vcf-iid";
if ($export eq "bcf"){
	$export_cmd = "--recode bcf";
}
if ($export eq "freq"){
	$export_cmd = "--freq";
}
if ($export eq "plink"){
	$export_cmd = "--make-bed";
}
if ($export eq "bed"){
        $export_cmd = "--make-bed";
} 


my $nb_alleles_cmd = "--min-alleles 1 --max-alleles 4";
if (@nalleles)
{
	$nb_alleles_cmd = "--min-alleles $nalleles[0] --max-alleles $nalleles[1]";
}
my $bounds_cmd = "";
if (@boundaries && $chrom_cmd !~/,/)
{
        #$bounds_cmd = "--from-bp $boundaries[0] --to-bp $boundaries[1]";
}


 
my $type_cmd = "";
if ($type eq "INDEL")
{
	$type_cmd = "--exclude-snp";
}
if ($type eq "SNP")
{
	$type_cmd = "--snps-only";
}

my $filt_cmd = "";
if ($remove_filt eq "true")
{
	$filt_cmd = "--remove-filtered-all";
}

my $thin_cmd = "";
if ($thin){
	$thin_cmd = "--bp-space $thin";
}

#my $bcf_input = $input;
#$bcf_input =~s/vcf/bcf/g;
my $bcf_input;
my $bed_input = $input;
$bed_input =~s/\.bed//g;

if (-e "$bed_input.bed"){
        system("$PLINK_EXE --bfile $bed_input --out $out $type_cmd $export_cmd $chrom_cmd $indiv_cmd $minfreq_cmd $maxfreq_cmd --geno $missing_data $thin_cmd $bounds_cmd --allow-extra-chr 1>$out.plink.stdout 2>$out.plink.stderr");
	# for first 1000 SNPs
	system("$PLINK_EXE --bfile $bed_input --out $out.recode $type_cmd --recode vcf-fid $chrom_cmd $indiv_cmd $minfreq_cmd $maxfreq_cmd --geno $missing_data $thin_cmd $bounds_cmd --allow-extra-chr --thin-count 800 1>$out.2.plink.stdout 2>$out.2.plink.stderr");
}
elsif (-e $bcf_input){
	system("$PLINK_EXE --bcf $bcf_input --out $out $type_cmd $export_cmd $chrom_cmd $indiv_cmd $minfreq_cmd $maxfreq_cmd --geno $missing_data $thin_cmd $bounds_cmd --allow-extra-chr 1>$out.plink.stdout 2>$out.plink.stderr");
}
else
{
	system("$PLINK_EXE --vcf $input --out $out $type_cmd $export_cmd $chrom_cmd $indiv_cmd $minfreq_cmd $maxfreq_cmd -geno $missing_data $thin_cmd $bounds_cmd --allow-extra-chr 1>$out.3.plink.stdout 2>$out.3.plink.stderr");

}




	
