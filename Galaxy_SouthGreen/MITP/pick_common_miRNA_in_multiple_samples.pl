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

my $version="1.0 version";

my $usage=<<USAGE; #******* Instruction of this program *********#
	Author: Wanfei Liu & Chengqi Xin
	Email: <liuwf\@big.ac.cn> & <xincq\@big.ac.cn>
	Date: Jun 8, 2013
	Version: $version

	This program can pick common miRNA according to the genome postion file of miRNA hairpin or mature sequence in gff format. Before you run it, you must make sure that different sample has different accession number for miRNA. Then you should cat all gff files of candidate miRNAs in different samples into one single file as input file.

	Usage: perl pick_common_miRNA_in_multiple_samples.pl --gff|-g <*.gff> --output|-o <*>
	
	The following options are necessary.
	--gff|-g	<*.gff>	:the gff file containing all miRNA in different samples
	--output|-o	<*>	:the output file

	The following options are optional.
	--overpercent|-hop	<hairoverpercent=50>	:minimum overlap percent(%), default is 50
	--help|-h				:print the usage information

USAGE

#Gather input
&GetOptions(
	"gff|s=s"	=>\$opt{gff},
	"output|o=s"	=>\$opt{output},
	"overpercent|op=i"	=>\$opt{overpercent},
	"help|h"	=>\$opt{help},
);

#Verify input
if (!defined $opt{gff} or !defined $opt{output} or defined $opt{help}) {
	die "$usage\n";
}

#Default parameters
$opt{overpercent}	= 50 unless defined $opt{overpercent};

# Absolute path
$opt{gff}         = abs_path $opt{gff};
$opt{output}      = abs_path $opt{output};

#Build commond line
my $cmdline = $0;
$cmdline .= " -g $opt{gff}" if defined $opt{gff};
$cmdline .= " -o $opt{output}" if defined $opt{output};
$cmdline .= " -op $opt{overpercent}" if defined $opt{overpercent};

print `date`, "\n>>> pick_common_miRNA_in_multiple_samples.pl <<<\n\n";
print $cmdline,"\n\n";

my $infile=$opt{gff};
my $outfile=$opt{output};

print "Read the input file.\n";

my (%repeat,);
open (IN,"<$infile")||die("fail to open $infile.\n");
while (<IN>) {
	chomp;
	#S091595 .       miRNA   194     211     .       -       .       ACC="MI1"; ID="ath-miR5021";
	my @list=split /\t/,$_;
	if ($list[8]=~/ACC\=\"([^\"]+)\"/) {
		$list[8]=~s/ACC\=//;
		$list[8]=~s/ID\=//;
		$list[8]=~s/EXPRESS\=//;
		$list[8]=~s/\"//g;
		$list[8]=~s/\;//g;
		my @feature=split /\s+/,$list[8];
		$repeat{$list[0]}{$list[6]}{"$list[3]\t$list[4]"}[0]=$list[3];
		$repeat{$list[0]}{$list[6]}{"$list[3]\t$list[4]"}[1].=$feature[0].",";
	}
}
close IN;

my $time=localtime;
print "The input file has been readed: $time.\n\n";

print "Compare the miRNA according to the genome position.\n";

open (OUT,">$outfile")||die("fail to open $outfile.\n");
foreach my $chr (keys %repeat) {
	foreach my $strand (keys %{$repeat{$chr}}) {
		my @record=sort {$repeat{$chr}{$strand}{$a}[0]<=>$repeat{$chr}{$strand}{$b}[0]} keys %{$repeat{$chr}{$strand}};
		for (my $i=0;$i<@record-1;$i++) {
			my ($start1,$end1)=split /\t/,$record[$i];
			my ($start2,$end2)=split /\t/,$record[$i+1];
			if (($end1>=$start2 and $end1<=$end2) or ($end2>=$start1 and $end2<=$end1)) {
				my $start=$start1;
				my $end=$end1;
				if ($start<$start2) {
					$start=$start2;
				}
				if ($end>$end2) {
					$end=$end2;
				}
				if (($end-$start+1)/($end1-$start1+1)>=$opt{overpercent}/100 or ($end-$start+1)/($end2-$start2+1)>=$opt{overpercent}/100) {
					$repeat{$chr}{$strand}{$record[$i]}[1].=$repeat{$chr}{$strand}{$record[$i+1]}[1];
					delete $repeat{$chr}{$strand}{$record[$i+1]};
					$record[$i+1]=$record[$i];
				}
			}
		}
		@record=sort {$repeat{$chr}{$strand}{$a}[0]<=>$repeat{$chr}{$strand}{$b}[0]} keys %{$repeat{$chr}{$strand}};
		for (my $i=0;$i<@record;$i++) {
			my @id=split /\,+/,$repeat{$chr}{$strand}{$record[$i]}[1];
			my $idnum=@id;
			print OUT "$idnum\t$repeat{$chr}{$strand}{$record[$i]}[1]\n";
		}
	}
}
close OUT;

$time=localtime;
print "The miRNA comparison has been finished: $time.\n\n";
