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

	This program can convert blat result to sam file.
	
	Usage: perl blat2sam.pl --blat|-b <blat_result> --output|-o <*>  --statistics|-s <*>
	
	The following options are necessary.
	--blat|-b	<blat_result>	:blat result file
	--output|-o			:output file
	--statistics|-s			:statistics file

	The following options are optional.
	--mismatch|-m		<mismatch=3>	:maximum mismatch value for mapping result (including indels), default is 3
	--identity|-i		<identity=90>	:minimum identity percent(%) (in mapped regions), default is 90
	--minimum|-min		<minimum=18>	:minimum miRNA length, default is 18
	--maximum|-max		<maximum=25>	:maximum miRNA length, default is 25
	--maxmap|-mm		<maxmap=10>	:maximum mapping position for each sequence, default is 10
	--help|-h		:print the usage information

USAGE

#Gather input
&GetOptions(
	"blat|b=s"	=>\$opt{blat},
	"output|o=s"	=>\$opt{output},
	"statistics|s=s"	=>\$opt{statistics},
	"mismatch|m=i"	=>\$opt{mismatch},
	"identity|i=i"	=>\$opt{identity},
	"minimum|min=i"	=>\$opt{minimum},
	"maximum|max=i"	=>\$opt{maximum},
	"maxmap|mm=i"	=>\$opt{maxmap},
	"help|h"	=>\$opt{help},
);

#Verify input
if (!defined $opt{blat} or !defined $opt{output} or !defined $opt{statistics} or defined $opt{help}) {
	die "$usage\n";
}

# Absolute path
$opt{blat}        = abs_path $opt{blat};
$opt{output}      = abs_path $opt{output};
$opt{statistics}      = abs_path $opt{statistics};

#Build commond line
my $cmdline = $0;
$cmdline .= " -b $opt{blat}" if defined $opt{blat};
$cmdline .= " -o $opt{output}" if defined $opt{output};
$cmdline .= " -s $opt{statistics}" if defined $opt{statistics};
$cmdline .= " -m $opt{mismatch}" if defined $opt{mismatch};
$cmdline .= " -i $opt{identity}" if defined $opt{identity};
$cmdline .= " -min $opt{minimum}" if defined $opt{minimum};
$cmdline .= " -max $opt{maximum}" if defined $opt{maximum};
$cmdline .= " -mm $opt{maxmap}" if defined $opt{maxmap};

print `date`, "\n>>> blat2sam.pl <<<\n\n";
print $cmdline,"\n\n";

print "Start to convert blat file into sam file.\n";

open (STAT,">>$opt{statistics}")||die("Cannot write to $opt{statistics}.\n");
open (OUT,">$opt{output}")||die("Cannot write to $opt{output}.\n");
open (IN,"<$opt{blat}")||die("Cannot read $opt{blat}.\n");
my $total_read=0;
my $filtered_read=0;
my $map=0;
my $id=undef;
while (<IN>) {
	chomp;
	next if (/^[\D]+/);
	my @list=split /\t/,$_;
	if (defined $id and $id ne $list[9]) {
		$total_read++;
		$filtered_read++ if ($map>0);
		$id=$list[9];
		$map=0;
		my $identity=$list[0]/($list[0]+$list[1])*100;
		my $mis=$list[1]+$list[5]+$list[7];
		if ($list[12]-$list[11]>=$opt{minmum} and $list[12]-$list[11]<=$opt{maximum} and $identity>=$opt{identity} and $mis<=$opt{mismatch}) {
			$map++;
			my $flag=&flag($list[8]);
			my $start=$list[15]+1;
			my $sam=($list[12]-$list[11])."M";
			my $NM="NM:i:".$mis;
			print OUT "$list[9]\t$flag\t$list[13]\t$start\t255\t$sam\t*\t0\t0\t*\t*\t$NM\n";
		}
	}else {
		$id=$list[9];
		next if ($map>=$opt{maxmap});
		my $identity=$list[0]/($list[0]+$list[1])*100;
		my $mis=$list[1]+$list[5]+$list[7];
		if ($list[12]-$list[11]>=$opt{minmum} and $list[12]-$list[11]<=$opt{maximum} and $identity>=$opt{identity} and $mis<=$opt{mismatch}) {
			$map++;
			my $flag=&flag($list[8]);
			my $start=$list[15]+1;
			my $sam=($list[12]-$list[11])."M";
			my $NM="NM:i:".$mis;
			print OUT "$list[9]\t$flag\t$list[13]\t$start\t255\t$sam\t*\t0\t0\t*\t*\t$NM\n";
		}
	}
}
$total_read++;
$filtered_read++ if ($map>0);
close IN;
close OUT;

print STAT "Total_read\t$total_read\nFiltered_read\t$filtered_read\n\n";
close STAT;

my $time=localtime;
print "The blat file has been converted to sam file.\n\n";

sub flag {
	## Bit    Description                                                Comment                                Value
    ## 0x1    template having multiple segments in sequencing            0: single-end 1: paired end            value: 2^^0 (  1)
    ## 0x2    each segment properly aligned according to the aligner     true only for paired-end alignments    value: 2^^1 (  2)
    ## 0x4    segment unmapped                                           ---                                           ---
    ## 0x8    next segment in the template unmapped                      ---                                           ---
    ## 0x10   SEQ being reverse complemented                                                                    value: 2^^4 ( 16)
    ## 0x20   SEQ of the next segment in the template being reversed                                            value: 2^^5 ( 32)
    ## 0x40   the first segment in the template                          read 1                                 value: 2^^6 ( 64)
    ## 0x80   the last segment in the template                           read 2                                 value: 2^^7 (128)
    ## 0x100  secondary alignment                                        ---                                           ---
    ## 0x200  not passing quality controls                               ---                                           ---
    ## 0x400  PCR or optical duplicate                                   ---                                           ---

	my $strand=shift;
	my $flag=0;

	if ($strand eq "+") {
		$flag+=0;
	}else {
		$flag+=16;
	}
	return $flag;
}
