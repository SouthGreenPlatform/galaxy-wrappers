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

	Usage: perl $program --masked_genome|-rf <repeat_masked_genome.fasta> --mat_fa|-mf <*_mature.fa> --mat_gff|-mg <*_mature.gff> --hair_fa|-hf <*_hairpin.fa> --hair_gff|-hg <*_hairpin.gff>
		
	The following options are necessary.
	--masked_genome|-rf	<repeat_masked_genome.fasta>
	--mat_fa|-mf	<*_mature.fa>	:candidate mature sequence file in fasta format
	--mat_gff|-mg	<*_mature.gff>	:candidate mature gff file in gff2 format
	--hair_fa|-hf	<*_hairpin.fa>	:candidate hairpin sequence file in fasta format
	--hair_gff|-hg	<*_hairpin.gff>	:candidate hairpin gff file in gff2 format

	The following options are optional.
	--output_dir|-o		<output_dir>	:default is sample under the current directory, you should change it if you have more than one sample to distinguish different samples
	--prefix|-p		<prefix=MI>	:prefix for file and miRNA ID, if you want to compare different samples, you should use different prefix, the prefix should be end with letters, default is MI
	--help|-h				:print the usage information

USAGE

#Gather input
&GetOptions(
	"masked_genome|rf=s"	=>\$opt{masked_genome},
	"mat_fa|mf=s"		=>\$opt{mat_fa},
	"mat_gff|mg=s"		=>\$opt{mat_gff},
	"hair_fa|hf=s"		=>\$opt{hair_fa},
	"hair_gff|hg=s"		=>\$opt{hair_gff},
	"output_dir|o=s"	=>\$opt{output_dir},
	"prefix|p=s"		=>\$opt{prefix},
	"help|h"	=>\$opt{help},
);

#Verify input
#Verify input
if (!defined $opt{masked_genome} or !defined $opt{mat_fa} or !defined $opt{mat_gff} or !defined $opt{hair_fa} or !defined $opt{hair_gff} or defined $opt{help}) {
	die "$usage\n";
}

#Default parameters
$opt{prefix}	= "MI" unless defined $opt{prefix};

# Absolute path
$opt{masked_genome}= abs_path $opt{masked_genome} if defined $opt{masked_genome};
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
$cmdline .= " -rf $opt{masked_genome}"  if defined $opt{masked_genome};
$cmdline .= " -mf $opt{mat_fa}" if defined $opt{mat_fa};
$cmdline .= " -mg $opt{mat_gff}" if defined $opt{mat_gff};
$cmdline .= " -hf $opt{hair_fa}" if defined $opt{hair_fa};
$cmdline .= " -hg $opt{hair_gff}" if defined $opt{hair_gff};
$cmdline .= " -o $opt{output_dir}"  if defined $opt{output_dir};
$cmdline .= " -p $opt{prefix}" if defined $opt{prefix};

print `date`, "\n>>> $program <<<\n\n";
print $cmdline,"\n\n";

my (%chr,$chr,$num,$time,);

#Read repeat masked genome fasta file
print "Read repeat masked genome fasta file.\n";

open (IN,"<$opt{masked_genome}")||die("Cannot read $opt{masked_genome}.\n");
while (<IN>) {
	chomp;
	if (/^>(\S+)/) {
		$chr=$1;
		next;
	}
	$chr{$chr}.=$_;
}
close IN;

$num=0;
open (OUT,">$opt{prefix}.repeat_region")||die("fail to open $opt{prefix}.repeat_region.\n");
open (IN,"<$opt{mat_gff}")||die("fail to open $opt{mat_gff}.\n");
while (<IN>) {
	chomp;
	#Chr10	MITP	miRNA	159065	159135	.	+	.	ACC="MI5"; ID="297,+"; EXPRESS="9";
	my @list=split /\t/,$_;
	if ($list[8]=~/ACC\=\"([^\"]+)\"/) {
		my $id=$1;
		my $id_seq=substr ($chr{$list[0]},$list[3]-1,abs($list[4]-$list[3])+1);
		my $id_base=($id_seq=~tr/ATGC/ATGC/);
		if ($id_base<length($id_seq)) {
			print OUT "$id\n";
			$num++;
		}
	}
}
close IN;
close OUT;

$time=localtime;
print "miRNA in repeat region was identified: $time.\n\n";

#Obtain repeat filtered mature and hairpin sequence and gff file for final candidate miRNA
print "Obtain repeat filtered mature and hairpin sequence and gff file for final candidate miRNA.\n";

system "perl $programdir/fish.pl -tb list -tf fa -c -b $opt{prefix}.repeat_region -f $opt{hair_fa} -o $opt{prefix}\_repeatfilter_hairpin.fa";
system "perl $programdir/fish.pl -tb list -tf fa -c -b $opt{prefix}.repeat_region -f $opt{mat_fa} -o $opt{prefix}\_repeatfilter_mature.fa";
system "perl $programdir/fish.pl -tb list -tf gff -c -b $opt{prefix}.repeat_region -f $opt{hair_gff} -o $opt{prefix}\_repeatfilter_hairpin.gff";
system "perl $programdir/fish.pl -tb list -tf gff -c -b $opt{prefix}.repeat_region -f $opt{mat_gff} -o $opt{prefix}\_repeatfilter_mature.gff";

$time=localtime;
print "The mature and hairpin sequence and gff file for repeat filtered candidate miRNA have been obtained: $time.\n\n";

open (STAT,">>$opt{prefix}.statistics")||die("Cannot write to $opt{prefix}.statistics.\n");
print STAT "Filter_repeat_region_miRNA\t$num\n\n";
close STAT;
