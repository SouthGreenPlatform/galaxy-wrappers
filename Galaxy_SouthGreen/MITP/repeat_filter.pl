#!/usr/bin/perl
use strict;
use warnings;
no warnings qw(uninitialized);
use Cwd 'abs_path';
use Getopt::Long;
use File::Basename;
#Variable setup
my %opt;
use vars qw/%opt/;

my $version="1.1 version";
my $program=$0;
$program=abs_path $program;
my @programpath = split /\//, $program;
my $programname = pop @programpath;
my $programdir  = join '/', @programpath;

my $usage=<<USAGE; #******* Instruction of this sub program *********#
	Author: Wanfei Liu & Chengqi Xin
	Email: <liuwf\@big.ac.cn> & <xincq\@big.ac.cn>
	Date: Nov 8, 2016
	Version: $version

	Usage: perl $program --repeat_fasta|-rf <repeat.fasta> --mat_fa|-mf <*_mature.fa> --mat_gff|-mg <*_mature.gff> --hair_fa|-hf <*_hairpin.fa> --hair_gff|-hg <*_hairpin.gff>
		
	The following options are necessary.
	--repeat_fasta|-rf	<repeat.fasta>
	--mat_fa|-mf	<*_mature.fa>	:candidate mature sequence file in fasta format
	--mat_gff|-mg	<*_mature.gff>	:candidate mature gff file in gff2 format
	--hair_fa|-hf	<*_hairpin.fa>	:candidate hairpin sequence file in fasta format
	--hair_gff|-hg	<*_hairpin.gff>	:candidate hairpin gff file in gff2 format

	The following options are optional.
	--output_dir|-o		<output_dir>	:default is sample under the current directory, you should change it if you have more than one sample to distinguish different samples
	--prefix|-p		<prefix=MI>	:prefix for file and miRNA ID, if you want to compare different samples, you should use different prefix, the prefix should be end with letters, default is MI
	--alignsoft|-as	<blast|blat>	:in default, blast was used.
	--identity|-i		<identity=90>	:minimum identity percent(%) (in mapped regions), default is 90	
	--filter_rate|-fr		<filter_rate=90>	:minimum filter overlap rate percent(%) between miRNA and repeat sequence, default is 50
	--help|-h				:print the usage information

USAGE

#Gather input
&GetOptions(
	"repeat_fasta|rf=s"	=>\$opt{repeat_fasta},
	"mat_fa|mf=s"		=>\$opt{mat_fa},
	"mat_gff|mg=s"		=>\$opt{mat_gff},
	"hair_fa|hf=s"		=>\$opt{hair_fa},
	"hair_gff|hg=s"		=>\$opt{hair_gff},
	"output_dir|o=s"	=>\$opt{output_dir},
	"prefix|p=s"		=>\$opt{prefix},
	"alignsoft|as=s"	=>\$opt{alignsoft},
	"identity|i=i"		=>\$opt{identity},	
	"filter_rate|fr=i"	=>\$opt{filter_rate},
	"help|h"	=>\$opt{help},
);

#Verify input
#Verify input
if (!defined $opt{repeat_fasta} or !defined $opt{mat_fa} or !defined $opt{mat_gff} or !defined $opt{hair_fa} or !defined $opt{hair_gff} or defined $opt{help}) {
	die "$usage\n";
}

#Default parameters
$opt{prefix}	= "MI" unless defined $opt{prefix};
$opt{alignsoft}	= "blast" unless defined $opt{alignsoft};
$opt{identity}	= 90 unless defined $opt{identity};
$opt{filter_rate}	= 90 unless defined $opt{filter_rate};

# Absolute path
$opt{repeat_fasta}= abs_path $opt{repeat_fasta} if defined $opt{repeat_fasta};
$opt{mat_fa}      = abs_path $opt{mat_fa} if defined $opt{mat_fa};
$opt{mat_gff}     = abs_path $opt{mat_gff} if defined $opt{mat_gff};
$opt{hair_fa}     = abs_path $opt{hair_fa} if defined $opt{hair_fa};
$opt{hair_gff}    = abs_path $opt{hair_gff} if defined $opt{hair_gff};

#Output directory
$opt{output_dir} = $ENV{'PWD'}."/sample"  unless defined $opt{output_dir};
$opt{output_dir} = abs_path $opt{output_dir} if defined $opt{output_dir};
`mkdir -p $opt{output_dir}` unless -d $opt{output_dir};
chdir $opt{output_dir};

#Build commond line
my $cmdline = $program;
$cmdline .= " -rf $opt{repeat_fasta}"  if defined $opt{repeat_fasta};
$cmdline .= " -mf $opt{mat_fa}" if defined $opt{mat_fa};
$cmdline .= " -mg $opt{mat_gff}" if defined $opt{mat_gff};
$cmdline .= " -hf $opt{hair_fa}" if defined $opt{hair_fa};
$cmdline .= " -hg $opt{hair_gff}" if defined $opt{hair_gff};
$cmdline .= " -o $opt{output_dir}"  if defined $opt{output_dir};
$cmdline .= " -p $opt{prefix}" if defined $opt{prefix};
$cmdline .= " -as $opt{alignsoft}" if defined $opt{alignsoft};
$cmdline .= " -i $opt{identity}" if defined $opt{identity};
$cmdline .= " -fr $opt{filter_rate}"  if defined $opt{filter_rate};

print `date`, "\n>>> $program <<<\n\n";
print $cmdline,"\n\n";

my (%repeat,$num,$time,);

#alignment with repeat sequence
print "Alignment with repeat sequence.\n";

if ($opt{alignsoft} eq "blast") {
	system "formatdb -i $opt{repeat_fasta} -p F";
	system "blastall -p blastn -d $opt{repeat_fasta} -i $opt{hair_fa} -v 1000 -b 1000 -W 7 -o $opt{prefix}\_candidate_hairpin.blast";
	system "perl $programdir/EblastN.pl -i $opt{prefix}\_candidate_hairpin.blast -o $opt{prefix}\_candidate_hairpin.eblastn";

	open (IN,"<$opt{prefix}\_candidate_hairpin.eblastn")||die("fail to open $opt{prefix}\_candidate_hairpin.eblastn.\n");
	while (<IN>) {
		chomp;
		next if (/^Query/);
		#MI34    18      1       18      1       18      19      36.2    5e-06   18/18   100     MI21405 ptc-miR169f
		my @list=split /\t/,$_;
		my $ident=$list[10];
		my $arate1=($list[3]-$list[2]+1)/$list[1];
		my $arate2=($list[5]-$list[4]+1)/$list[6];	
		if ($ident>=$opt{identity} and ($arate1>=$opt{filter_rate} or $arate2>=$opt{filter_rate})) {
			$repeat{$list[0]}++;
		}
	}
	close IN;
}else {
	system "blat $opt{repeat_fasta} $opt{hair_fa} -tileSize=8 -oneOff=1 -minMatch=1 -minIdentity=80 -noHead $opt{prefix}\_candidate_hairpin.blat";

	open (IN,"<$opt{prefix}\_candidate_hairpin.blat")||die("Cannot read $opt{prefix}\_candidate_hairpin.blat.\n");

	while (<IN>) {
		chomp;
		my @list=split /\t/,$_;
		my $ident=$list[0]/($list[0]+$list[1])*100;
		my $arate1=($list[12]-$list[11])/$list[10];
		my $arate2=($list[16]-$list[15])/$list[14];
		if ($ident>=$opt{identity} and ($arate1>=$opt{filter_rate} or $arate2>=$opt{filter_rate})) {
			$repeat{$list[0]}++;
		}
	}
	close IN;
}

open (OUT,">$opt{prefix}\_candidate_hairpin.filter")||die("fail to open $opt{prefix}\_candidate_hairpin.filter.\n");
$num=0;
foreach my $key (keys %repeat) {
	print OUT "$key\n";
	$num++;
}
close OUT;

$time=localtime;
print "The candidate hairpin sequence was aligned to repeat sequence: $time.\n\n";

#Obtain repeat filtered mature and hairpin sequence and gff file for final candidate miRNA
print "Obtain repeat filtered mature and hairpin sequence and gff file for final candidate miRNA.\n";

system "perl $programdir/fish.pl -tb list -tf fa -c -b $opt{prefix}\_candidate_hairpin.filter -f $opt{hair_fa} -o $opt{prefix}\_repeatfilter_hairpin.fa";
system "perl $programdir/fish.pl -tb list -tf fa -c -b $opt{prefix}\_candidate_hairpin.filter -f $opt{mat_fa} -o $opt{prefix}\_repeatfilter_mature.fa";
system "perl $programdir/fish.pl -tb list -tf gff -c -b $opt{prefix}\_candidate_hairpin.filter -f $opt{hair_gff} -o $opt{prefix}\_repeatfilter_hairpin.gff";
system "perl $programdir/fish.pl -tb list -tf gff -c -b $opt{prefix}\_candidate_hairpin.filter -f $opt{mat_gff} -o $opt{prefix}\_repeatfilter_mature.gff";

$time=localtime;
print "The mature and hairpin sequence and gff file for repeat filtered candidate miRNA have been obtained: $time.\n\n";

open (STAT,">>$opt{prefix}.statistics")||die("Cannot write to $opt{prefix}.statistics.\n");
print STAT "Filter_repeat_region_miRNA\t$num\n\n";
close STAT;
