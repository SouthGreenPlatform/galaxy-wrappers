#!/usr/bin/perl
use strict;
use warnings;
#no warnings qw(uninitialized);
use Cwd 'abs_path';
use Getopt::Long;
use File::Basename;
#Variable setup
my %opt;
use vars qw/%opt/;

my $version="1.1 version";

my $usage=<<USAGE; #******* Instruction of this sub program *********#
	Author: Wanfei Liu & Chengqi Xin
	Email: <liuwf\@big.ac.cn> & <xincq\@big.ac.cn>
	Date: Jun 8, 2013
	Version: $version

	Usage: perl extract_seq.pl --genome|-g <*.fa> --filtered_cluster|-fc <*_read.cluster_filter>
	
	The following options are necessary.
	--genome|-g	<*.fa>	:genome sequence file in fasta format
	--filtered_cluster|-fc	<*_read.cluster_filter>	:the filtered read cluster file

	The following options are optional.
	--output_dir|-o		<output_dir>	:default is sample under the current directory, you should change it if you have more than one sample to distinguish different samples
	--prefix|-p		<prefix=MI>	:prefix for file and miRNA ID, if you want to compare different samples, you should use different prefix, the prefix should be end with letters, default is MI
	--strand_specific|-ss			:if you assign this parameter, it means the data is strand specific, if not, as default, it means the data is not strand specific
	--mincov|-mc		<mincov=10>	:minimum miRNA mapping coverage, default is 10 (for conserve miRNA, it should be 1, for novel miRNA from expression, it should be minimum expressed read number)
	--help|-h				:print the usage information

USAGE

#Gather input
&GetOptions(
	"genome|g=s"		=>\$opt{genome},
	"filtered_cluster|fc=s"		=>\$opt{filtered_cluster},
	"output_dir|o=s"	=>\$opt{output_dir},
	"prefix|p=s"		=>\$opt{prefix},
	"strand_specific|ss"=>\$opt{strand_specific},
	"mincov|mc=s"		=>\$opt{mincov},
	"help|h"	=>\$opt{help},
);

#Verify input
if (!defined $opt{genome} or !defined $opt{filtered_cluster} or defined $opt{help}) {
	die "$usage\n";
}

#Default parameters
$opt{prefix}	= "MI" unless defined $opt{prefix};
$opt{mincov}	= 10 unless defined $opt{mincov};

# Absolute path
$opt{genome}			= abs_path $opt{genome} if defined $opt{genome};
$opt{filtered_cluster}	= abs_path $opt{filtered_cluster} if defined $opt{filtered_cluster};

#Output directory
$opt{output_dir} = $ENV{'PWD'}."/sample"  unless defined $opt{output_dir};
$opt{output_dir} = abs_path $opt{output_dir} if defined $opt{output_dir};
`mkdir -p $opt{output_dir}` unless -d $opt{output_dir};
chdir $opt{output_dir};

#Build commond line
my $cmdline = $0;
$cmdline .= " -g $opt{genome}" if defined $opt{genome};
$cmdline .= " -fc $opt{filtered_cluster}" if defined $opt{filtered_cluster};
$cmdline .= " -o $opt{output_dir}"  if defined $opt{output_dir};
$cmdline .= " -p $opt{prefix}" if defined $opt{prefix};
$cmdline .= " -ss $opt{strand_specific}" if defined $opt{strand_specific};
$cmdline .= " -mc $opt{mincov}" if defined $opt{mincov};

print `date`, "\n>>> extract_seq.pl <<<\n\n";
print $cmdline,"\n\n";

my ($chr,$string,%gen,%len,$time,);

#Obtain chromosome sequence and length
print "Obtain chromosome sequence and length.\n";

open (IN,"<$opt{genome}")||die("Cannot read $opt{genome}.\n");
while (<IN>) {
	chomp;
	if (/^\#/) {
		next;
	}elsif (/^>(\S+)/) {
		if (defined $chr) {
			$gen{$chr}=$string;
			$len{$chr}=length $string;
			$string=undef;
		}
		$chr=$1;
		next;
	}
	$string.=$_;
}
$gen{$chr}=$string;
$len{$chr}=length $string;
$string=undef;
close IN;

$time=localtime;
print "Chromosome sequence and length are obtained: $time.\n\n";

#Creat mature and precursor sequence and gff file for candidate sequence
print "Create mature and precursor sequence and gff files.\n";
open (OUT1,">$opt{prefix}\_pre.fa")||die("Cannot write to $opt{prefix}\_pre.fa file.\n");
open (OUT2,">$opt{prefix}\_mature.fa")||die("Cannot write to $opt{prefix}\_mature.fa file.\n");
open (OUT3,">$opt{prefix}\_pre.gff")||die("Cannot write to $opt{prefix}\_pre.gff file.\n");
open (OUT4,">$opt{prefix}\_mature.gff")||die("Cannot write to $opt{prefix}\_mature.gff file.\n");

my $num=1;
open (IN,"<$opt{filtered_cluster}")||die("Cannot read $opt{filtered_cluster}.\n");
while (<IN>) {
	chomp;
	#aau-miR160      21      1       21      1137925 1137905 2539026 34.2    0.15    20/21   95      S000014
	next if (/^\#/);
	my @list=split /\t/,$_;
	next if ($list[12]<$opt{mincov});
	if ($list[5] > $list[4]) {
		my $mat_start=$list[4];
		my $mat_end=$list[5];
		my $pre_start=();
		if ($list[4]<=100) {
			$pre_start=1;
		}else {
			$pre_start=$list[4]-100;
		}
		my $pre_end=$list[5]+100;
		my $mat_seq=substr ($gen{$list[11]},$mat_start-1,abs($mat_end-$mat_start)+1);
		my $pre_seq=substr ($gen{$list[11]},$pre_start-1,abs($pre_end-$pre_start)+1);
		if (defined $opt{strand_specific}) {
			print OUT1 ">$opt{prefix}$num $list[0]\n$pre_seq\n";
			print OUT2 ">$opt{prefix}$num $list[0]\n$mat_seq\n";
			print OUT3 "$list[11]\tMIP\tmiRNA\t$pre_start\t$pre_end\t\.\t+\t\.\tACC=\"$opt{prefix}$num\"\; ID=\"$list[0],+\"\; EXPRESS=\"$list[$#list]\"\;\n";
			print OUT4 "$list[11]\tMIP\tmiRNA\t$mat_start\t$mat_end\t\.\t+\t\.\tACC=\"$opt{prefix}$num\"\; ID=\"$list[0],+\"\; EXPRESS=\"$list[$#list]\"\;\n";
			$num++;
		}else {
			print OUT1 ">$opt{prefix}$num $list[0]\n$pre_seq\n";
			print OUT2 ">$opt{prefix}$num $list[0]\n$mat_seq\n";
			print OUT3 "$list[11]\tMIP\tmiRNA\t$pre_start\t$pre_end\t\.\t+\t\.\tACC=\"$opt{prefix}$num\"\; ID=\"$list[0],+\"\; EXPRESS=\"$list[$#list]\"\;\n";
			print OUT4 "$list[11]\tMIP\tmiRNA\t$mat_start\t$mat_end\t\.\t+\t\.\tACC=\"$opt{prefix}$num\"\; ID=\"$list[0],+\"\; EXPRESS=\"$list[$#list]\"\;\n";
			$num++;
			$mat_seq=reverse $mat_seq;
			$mat_seq=~tr/ACGTacgt/TGCAtgca/;
			$pre_seq=reverse $pre_seq;
			$pre_seq=~tr/ACGTacgt/TGCAtgca/;
			print OUT1 ">$opt{prefix}$num $list[0]\n$pre_seq\n";
			print OUT2 ">$opt{prefix}$num $list[0]\n$mat_seq\n";
			print OUT3 "$list[11]\tMIP\tmiRNA\t$pre_start\t$pre_end\t\.\t-\t\.\tACC=\"$opt{prefix}$num\"\; ID=\"$list[0],+\"\; EXPRESS=\"$list[$#list]\"\;\n";
			print OUT4 "$list[11]\tMIP\tmiRNA\t$mat_start\t$mat_end\t\.\t-\t\.\tACC=\"$opt{prefix}$num\"\; ID=\"$list[0],+\"\; EXPRESS=\"$list[$#list]\"\;\n";
			$num++;
		}
	}else {
		my $mat_start=$list[5];
		my $mat_end=$list[4];
		my $pre_start=();
		if ($list[5]<=100) {
			$pre_start=1;
		}else {
			$pre_start=$list[5]-100;
		}
		my $pre_end=$list[4]+100;
		my $mat_seq=substr ($gen{$list[11]},$mat_start-1,abs($mat_end-$mat_start)+1);
		my $pre_seq=substr ($gen{$list[11]},$pre_start-1,abs($pre_end-$pre_start)+1);
		if (defined $opt{strand_specific}) {
			$mat_seq=reverse $mat_seq;
			$mat_seq=~tr/ACGTacgt/TGCAtgca/;
			$pre_seq=reverse $pre_seq;
			$pre_seq=~tr/ACGTacgt/TGCAtgca/;
			print OUT1 ">$opt{prefix}$num $list[0]\n$pre_seq\n";
			print OUT2 ">$opt{prefix}$num $list[0]\n$mat_seq\n";
			print OUT3 "$list[11]\tMIP\tmiRNA\t$pre_start\t$pre_end\t\.\t-\t\.\tACC=\"$opt{prefix}$num\"\; ID=\"$list[0],-\"\; EXPRESS=\"$list[$#list]\"\;\n";
			print OUT4 "$list[11]\tMIP\tmiRNA\t$mat_start\t$mat_end\t\.\t-\t\.\tACC=\"$opt{prefix}$num\"\; ID=\"$list[0],-\"\; EXPRESS=\"$list[$#list]\"\;\n";
			$num++;
		}else {
			print OUT1 ">$opt{prefix}$num $list[0]\n$pre_seq\n";
			print OUT2 ">$opt{prefix}$num $list[0]\n$mat_seq\n";
			print OUT3 "$list[11]\tMIP\tmiRNA\t$pre_start\t$pre_end\t\.\t+\t\.\tACC=\"$opt{prefix}$num\"\; ID=\"$list[0],-\"\; EXPRESS=\"$list[$#list]\"\;\n";
			print OUT4 "$list[11]\tMIP\tmiRNA\t$mat_start\t$mat_end\t\.\t+\t\.\tACC=\"$opt{prefix}$num\"\; ID=\"$list[0],-\"\; EXPRESS=\"$list[$#list]\"\;\n";
			$num++;
			$mat_seq=reverse $mat_seq;
			$mat_seq=~tr/ACGTacgt/TGCAtgca/;
			$pre_seq=reverse $pre_seq;
			$pre_seq=~tr/ACGTacgt/TGCAtgca/;
			print OUT1 ">$opt{prefix}$num $list[0]\n$pre_seq\n";
			print OUT2 ">$opt{prefix}$num $list[0]\n$mat_seq\n";
			print OUT3 "$list[11]\tMIP\tmiRNA\t$pre_start\t$pre_end\t\.\t-\t\.\tACC=\"$opt{prefix}$num\"\; ID=\"$list[0],-\"\; EXPRESS=\"$list[$#list]\"\;\n";
			print OUT4 "$list[11]\tMIP\tmiRNA\t$mat_start\t$mat_end\t\.\t-\t\.\tACC=\"$opt{prefix}$num\"\; ID=\"$list[0],-\"\; EXPRESS=\"$list[$#list]\"\;\n";
			$num++;
		}
	}
}
close IN;
close OUT1;
close OUT2;
close OUT3;
close OUT4;

$num=$num-1;

open (STAT,">>$opt{prefix}.statistics")||die("Cannot write to $opt{prefix}.statistics.\n");
print STAT "Extracted_seq\t$num\n\n";
close STAT;

$time=localtime;
print "The mature and precursor sequence and gff files have been created: $time.\n\n";
