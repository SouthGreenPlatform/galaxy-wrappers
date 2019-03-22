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
my $program=$0;
$program=abs_path $program;
my @programpath = split /\//, $program;
my $programname = pop @programpath;
my $programdir  = join '/', @programpath;

my $usage=<<USAGE; #******* Instruction of this sub program *********#
	Author: Wanfei Liu & Chengqi Xin
	Email: <liuwf\@big.ac.cn> & <xincq\@big.ac.cn>
	Date: Jun 8, 2013
	Version: $version

	Note: for identification miRNA, the parameter for blast and blat should change comparison to long sequence mapping. In our opinion, the blast command should be as \"blastall -p blastn -d database -i query -v 1000 -b 1000 -W 7 -o outfile\" and the blat command should be as \"blat database query -tileSize=8 -oneOff=1 -minMatch=1 -minIdentity=80 -noHead outfile\".

	Usage: perl $program { --sam <*> | --blast <*> | --blat <*> } --genome|-g <*.fa>
	
	The following options are necessary.
	--sam	<*>	:mapping result file in sam format
	--blast	<*>	:mapping result file in blast default output format
	--blat	<*>	:mapping result file in blat default output format
	--genome|-g	<*.fa>	:genome sequence file in fasta format

	The following options are optional.
	--output_dir|-o		<output_dir>	:default is sample under the current directory, you should change it if you have more than one sample to distinguish different samples
	--prefix|-p		<prefix=MI>	:prefix for file and miRNA ID, if you want to compare different samples, you should use different prefix, the prefix should be end with letters, default is MI
	--strand_specific|-ss			:if you assign this parameter, it means the data is strand specific, if not, as default, it means the data is not strand specific
	--maxmap|-mm		<maxmap=10>	:maximum mapping position for each sequence, default is 10
	--mismatch|-m		<mismatch=3>	:maximum mismatch value for mapping result (including indels), default is 3
	--identity|-i		<identity=90>	:minimum identity percent(%) (in mapped regions), default is 90	
	--minimum|-min		<minimum=18>	:minimum miRNA length, default is 18
	--maximum|-max		<maximum=25>	:maximum miRNA length, default is 25
	--oversize|-os		<oversize=16>	:minimum overlap size for read cluster, default is 16
	--overrate|-or		<overrate=80>	:minimum overlap rate percent(%) for read cluster, default is 80
	--help|-h				:print the usage information

USAGE

#Gather input
&GetOptions(
	"sam=s"				=>\$opt{sam},
	"blast=s"			=>\$opt{blast},
	"blat=s"			=>\$opt{blat},
	"genome|g=s"		=>\$opt{genome},
	"output_dir|o=s"	=>\$opt{output_dir},
	"prefix|p=s"		=>\$opt{prefix},
	"strand_specific|ss"=>\$opt{strand_specific},
	"maxmap|mm=i"		=>\$opt{maxmap},
	"mismatch|m=i"		=>\$opt{mismatch},
	"identity|i=i"		=>\$opt{identity},	
	"minimum|min=i"		=>\$opt{minimum},
	"maximum|max=i"		=>\$opt{maximum},
	"oversize|os=i"		=>\$opt{oversize},
	"overrate|or=i"		=>\$opt{overrate},
	"help|h"	=>\$opt{help},
);

#Verify input
if ((!defined $opt{sam} and !defined $opt{blast} and !defined $opt{blat}) or !defined $opt{genome} or defined $opt{help}) {
	die "$usage\n";
}

#Default parameters
$opt{prefix}	= "MI" unless defined $opt{prefix};
$opt{maxmap}	= 10 unless defined $opt{maxmap};
$opt{mismatch}	= 3 unless defined $opt{mismatch};
$opt{identity}	= 90 unless defined $opt{identity};
$opt{minimum}	= 18 unless defined $opt{minimum};
$opt{maximum}	= 25 unless defined $opt{maximum};
$opt{oversize}	= 16 unless defined $opt{oversize};
$opt{overrate}	= 80 unless defined $opt{overrate};

# Absolute path
$opt{sam}         = abs_path $opt{sam} if defined $opt{sam};
$opt{blast}       = abs_path $opt{blast} if defined $opt{blast};
$opt{blat}        = abs_path $opt{blat} if defined $opt{blat};
$opt{genome}      = abs_path $opt{genome};

#Output directory
$opt{output_dir} = $ENV{'PWD'}."/sample"  unless defined $opt{output_dir};
$opt{output_dir} = abs_path $opt{output_dir} if defined $opt{output_dir};
`mkdir -p $opt{output_dir}` unless -d $opt{output_dir};
chdir $opt{output_dir};

#Build commond line
my $cmdline = $program;
$cmdline .= " --sam $opt{sam}" if defined $opt{sam};
$cmdline .= " --blast $opt{blast}" if defined $opt{blast};
$cmdline .= " --blat $opt{blat}" if defined $opt{blat};
$cmdline .= " -g $opt{genome}" if defined $opt{genome};
$cmdline .= " -o $opt{output_dir}"  if defined $opt{output_dir};
$cmdline .= " -p $opt{prefix}" if defined $opt{prefix};
$cmdline .= " -ss $opt{strand_specific}" if defined $opt{strand_specific};
$cmdline .= " -mm $opt{maxmap}" if defined $opt{maxmap};
$cmdline .= " -m $opt{mismatch}" if defined $opt{mismatch};
$cmdline .= " -i $opt{identity}" if defined $opt{identity};
$cmdline .= " -min $opt{minimum}" if defined $opt{minimum};
$cmdline .= " -max $opt{maximum}" if defined $opt{maximum};
$cmdline .= " -os $opt{oversize}" if defined $opt{oversize};
$cmdline .= " -or $opt{overrate}" if defined $opt{overrate};

print `date`, "\n>>> $program <<<\n\n";
print $cmdline,"\n\n";

my ($chr,$string,%gen,%len,%id,%hash,);
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

my $time=localtime;
print "Chromosome sequence and length are obtained: $time.\n\n";

%hash=();
#Read the mapping result file
print "Read the mapping result file.\n";

if (defined $opt{blast}) {
	system "perl $programdir/EblastN.pl -i $opt{blast} -o $opt{prefix}.eblastn";
	system "perl $programdir/eblastn2sam.pl -e $opt{prefix}.eblastn -o $opt{prefix}.sam -s $opt{prefix}.statistics -m $opt{mismatch} -i $opt{identity} -min $opt{minimum} -max $opt{maximum} -mm $opt{maxmap}";

	open (IN,"<$opt{prefix}.sam")||die("Cannot read $opt{prefix}.sam.\n");
	while (<IN>) {
		chomp;
		#SRR027895.5259200_1     99      chr1    15443   255     76M     =       15623   256     CTCCAGAGGCCTCAGGTCCAGTCTCTAAAAATATCTCAGGAAGCTGCAGTGGCTGACCATTGCCTTGGACCGCTCT    BCBCCBBABCBBBBBB=BBAB?BBBBBAABBBBBBBBBBA@BBBABBBB?@AAA>?<702>7:;1:+670:9;<==    NM:i:1  NH:i:1
		next if(/^\@/);
		my @list=split /\t/,$_;
		next if ($list[2] eq "*");
		my $strand=&flag($list[1]);
		my $chr=$list[2];
		my $start=$list[3];	
		my @match=split /([M|I|D|S|H]+)/,$list[5];
		my $mis=0;#include indel
		my $mis2=0;#in mapped read
		my $len=0;#in reference
		my $len2=0;#in mapped read
		for (my $i=0;$i<@match-1;$i+=2) {
			if ($match[$i+1] eq "M") {
				$len+=$match[$i];
				$len2+=$match[$i];
			}elsif ($match[$i+1] eq "D") {
				$len+=$match[$i];
				$mis+=$match[$i];
			}elsif ($match[$i+1] eq "I") {
				$mis+=$match[$i];
			}
		}
		if ($_=~/NM\:i\:(\d+)/) {
			$mis+=$1;
			$mis2=$1;
		}
		my $end=$start+$len-1;
		if (defined $opt{strand_specific}) {
			if (!exists $hash{$chr}{$strand}{$start."\t".$end}) {
				$hash{$chr}{$strand}{$start."\t".$end}[0]=$len-$mis;
				$hash{$chr}{$strand}{$start."\t".$end}[1]++;
				$hash{$chr}{$strand}{$start."\t".$end}[2]=$start;
				$hash{$chr}{$strand}{$start."\t".$end}[3]++;
			}else {
				$hash{$chr}{$strand}{$start."\t".$end}[1]++;
				$hash{$chr}{$strand}{$start."\t".$end}[2]=$start;
				$hash{$chr}{$strand}{$start."\t".$end}[3]++;
				if ($len-$mis>$hash{$chr}{$strand}{$start."\t".$end}[0]) {
					$hash{$chr}{$strand}{$start."\t".$end}[0]=$len-$mis;
				}
			}
		}else {
			if (!exists $hash{$chr}{$start."\t".$end}) {
				$hash{$chr}{$start."\t".$end}[0]=$len-$mis;
				$hash{$chr}{$start."\t".$end}[1]++;
				$hash{$chr}{$start."\t".$end}[2]=$start;
				$hash{$chr}{$start."\t".$end}[3]++;
			}else {
				$hash{$chr}{$start."\t".$end}[1]++;
				$hash{$chr}{$start."\t".$end}[2]=$start;
				$hash{$chr}{$start."\t".$end}[3]++;
				if ($len-$mis>$hash{$chr}{$start."\t".$end}[0]) {
					$hash{$chr}{$start."\t".$end}[0]=$len-$mis;
				}
			}
		}
	}
	close IN;
}elsif (defined $opt{blat}) {
	system "perl $programdir/blat2sam.pl -b $opt{blat} -o $opt{prefix}.sam -s $opt{prefix}.statistics -m $opt{mismatch} -i $opt{identity} -min $opt{minimum} -max $opt{maximum} -mm $opt{maxmap}";

	open (IN,"<$opt{prefix}.sam")||die("Cannot read $opt{prefix}.sam.\n");
	while (<IN>) {
		chomp;
		#SRR027895.5259200_1     99      chr1    15443   255     76M     =       15623   256     CTCCAGAGGCCTCAGGTCCAGTCTCTAAAAATATCTCAGGAAGCTGCAGTGGCTGACCATTGCCTTGGACCGCTCT    BCBCCBBABCBBBBBB=BBAB?BBBBBAABBBBBBBBBBA@BBBABBBB?@AAA>?<702>7:;1:+670:9;<==    NM:i:1  NH:i:1
		next if(/^\@/);
		my @list=split /\t/,$_;
		next if ($list[2] eq "*");
		my $strand=&flag($list[1]);
		my $chr=$list[2];
		my $start=$list[3];	
		my @match=split /([M|I|D|S|H]+)/,$list[5];
		my $mis=0;#include indel
		my $mis2=0;#in mapped read
		my $len=0;#in reference
		my $len2=0;#in mapped read
		for (my $i=0;$i<@match-1;$i+=2) {
			if ($match[$i+1] eq "M") {
				$len+=$match[$i];
				$len2+=$match[$i];
			}elsif ($match[$i+1] eq "D") {
				$len+=$match[$i];
				$mis+=$match[$i];
			}elsif ($match[$i+1] eq "I") {
				$mis+=$match[$i];
			}
		}
		if ($_=~/NM\:i\:(\d+)/) {
			$mis+=$1;
			$mis2=$1;
		}
		my $end=$start+$len-1;
		if (defined $opt{strand_specific}) {
			if (!exists $hash{$chr}{$strand}{$start."\t".$end}) {
				$hash{$chr}{$strand}{$start."\t".$end}[0]=$len-$mis;
				$hash{$chr}{$strand}{$start."\t".$end}[1]++;
				$hash{$chr}{$strand}{$start."\t".$end}[2]=$start;
				$hash{$chr}{$strand}{$start."\t".$end}[3]++;
			}else {
				$hash{$chr}{$strand}{$start."\t".$end}[1]++;
				$hash{$chr}{$strand}{$start."\t".$end}[2]=$start;
				$hash{$chr}{$strand}{$start."\t".$end}[3]++;
				if ($len-$mis>$hash{$chr}{$strand}{$start."\t".$end}[0]) {
					$hash{$chr}{$strand}{$start."\t".$end}[0]=$len-$mis;
				}
			}
		}else {
			if (!exists $hash{$chr}{$start."\t".$end}) {
				$hash{$chr}{$start."\t".$end}[0]=$len-$mis;
				$hash{$chr}{$start."\t".$end}[1]++;
				$hash{$chr}{$start."\t".$end}[2]=$start;
				$hash{$chr}{$start."\t".$end}[3]++;
			}else {
				$hash{$chr}{$start."\t".$end}[1]++;
				$hash{$chr}{$start."\t".$end}[2]=$start;
				$hash{$chr}{$start."\t".$end}[3]++;
				if ($len-$mis>$hash{$chr}{$start."\t".$end}[0]) {
					$hash{$chr}{$start."\t".$end}[0]=$len-$mis;
				}
			}
		}
	}
	close IN;
}elsif (defined $opt{sam}) {
	open (STAT,">>$opt{prefix}.statistics")||die("Cannot write to $opt{prefix}.statistics.\n");
	open (IN,"<$opt{sam}")||die("Cannot read $opt{sam}.\n");
	%id=();
	while (<IN>) {
		chomp;
		#SRR027895.5259200_1     99      chr1    15443   255     76M     =       15623   256     CTCCAGAGGCCTCAGGTCCAGTCTCTAAAAATATCTCAGGAAGCTGCAGTGGCTGACCATTGCCTTGGACCGCTCT    BCBCCBBABCBBBBBB=BBAB?BBBBBAABBBBBBBBBBA@BBBABBBB?@AAA>?<702>7:;1:+670:9;<==    NM:i:1  NH:i:1
		next if(/^\@/);
		my @list=split /\t/,$_;
		next if ($list[2] eq "*");
		if (!defined $id{$list[0]}) {
			$id{$list[0]}=0;
		}else {
			next if ($id{$list[0]}>=$opt{maxmap});
		}
		my $strand=&flag($list[1]);
		my $chr=$list[2];
		my $start=$list[3];	
		my @match=split /([M|I|D|S|H]+)/,$list[5];
		my $mis=0;#include indel
		my $mis2=0;#in mapped read
		my $len=0;#in reference
		my $len2=0;#in mapped read
		for (my $i=0;$i<@match-1;$i+=2) {
			if ($match[$i+1] eq "M") {
				$len+=$match[$i];
				$len2+=$match[$i];
			}elsif ($match[$i+1] eq "D") {
				$len+=$match[$i];
				$mis+=$match[$i];
			}elsif ($match[$i+1] eq "I") {
				$mis+=$match[$i];
			}
		}
		if ($_=~/NM\:i\:(\d+)/) {
			$mis+=$1;
			$mis2=$1;
		}
		my $end=$start+$len-1;
		my $identity=($len2-$mis2)/$len2*100;
		if ($len>=$opt{minimum} and $len<=$opt{maximum} and $identity>=$opt{identity} and $mis<=$opt{mismatch}) {
			$id{$list[0]}++;
			if (defined $opt{strand_specific}) {
				if (!exists $hash{$chr}{$strand}{$start."\t".$end}) {
					$hash{$chr}{$strand}{$start."\t".$end}[0]=$len-$mis;
					$hash{$chr}{$strand}{$start."\t".$end}[1]++;
					$hash{$chr}{$strand}{$start."\t".$end}[2]=$start;
					$hash{$chr}{$strand}{$start."\t".$end}[3]++;
				}else {
					$hash{$chr}{$strand}{$start."\t".$end}[1]++;
					$hash{$chr}{$strand}{$start."\t".$end}[2]=$start;
					$hash{$chr}{$strand}{$start."\t".$end}[3]++;
					if ($len-$mis>$hash{$chr}{$strand}{$start."\t".$end}[0]) {
						$hash{$chr}{$strand}{$start."\t".$end}[0]=$len-$mis;
					}
				}
			}else {
				if (!exists $hash{$chr}{$start."\t".$end}) {
					$hash{$chr}{$start."\t".$end}[0]=$len-$mis;
					$hash{$chr}{$start."\t".$end}[1]++;
					$hash{$chr}{$start."\t".$end}[2]=$start;
					$hash{$chr}{$start."\t".$end}[3]++;
				}else {
					$hash{$chr}{$start."\t".$end}[1]++;
					$hash{$chr}{$start."\t".$end}[2]=$start;
					$hash{$chr}{$start."\t".$end}[3]++;
					if ($len-$mis>$hash{$chr}{$start."\t".$end}[0]) {
						$hash{$chr}{$start."\t".$end}[0]=$len-$mis;
					}
				}
			}
		}
	}
	close IN;

	my $total_read=0;
	my $filtered_read=0;
	foreach my $id (keys %id) {
		if ($id{$id}>0) {
			$filtered_read++;
		}
		$total_read++;
	}
	%id=();

	print STAT "Total_read\t$total_read\nFiltered_read\t$filtered_read\n";
	close STAT;
}

$time=localtime;
print "The mapping result has been readed: $time.\n\n";

#Cluster mapping result
print "Cluster mapping result.\n";
foreach my $chr (keys %hash) {
	if (defined $opt{strand_specific}) {
		foreach my $strand (keys %{$hash{$chr}}) {
			my @record=sort {$hash{$chr}{$strand}{$a}[2]<=>$hash{$chr}{$strand}{$b}[2]} keys %{$hash{$chr}{$strand}};
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
					if ($end-$start+1>=$opt{oversize} and (($end-$start+1)/($end1-$start1+1)>=$opt{overrate}/100 or ($end-$start+1)/($end2-$start2+1)>=$opt{overrate}/100)) {
						my $match1=$hash{$chr}{$strand}{$record[$i]}[0];
						my $read1=$hash{$chr}{$strand}{$record[$i]}[1];
						my $match2=$hash{$chr}{$strand}{$record[$i+1]}[0];
						my $read2=$hash{$chr}{$strand}{$record[$i+1]}[1];
						if ($read1>$read2) {
							$hash{$chr}{$strand}{$record[$i]}[0]=$match1;
							$hash{$chr}{$strand}{$record[$i]}[1]=$read1;
							$hash{$chr}{$strand}{$record[$i]}[2]=$start1;
							$hash{$chr}{$strand}{$record[$i]}[3]+=$hash{$chr}{$strand}{$record[$i+1]}[3];
							@{$hash{$chr}{$strand}{$record[$i+1]}}=();
							delete $hash{$chr}{$strand}{$record[$i+1]};
							$record[$i+1]=$record[$i];
						}elsif ($read1==$read2) {
							if ($match1>=$match2) {
								$hash{$chr}{$strand}{$record[$i]}[0]=$match1;
								$hash{$chr}{$strand}{$record[$i]}[1]=$read1;
								$hash{$chr}{$strand}{$record[$i]}[2]=$start1;
								$hash{$chr}{$strand}{$record[$i]}[3]+=$hash{$chr}{$strand}{$record[$i+1]}[3];
								@{$hash{$chr}{$strand}{$record[$i+1]}}=();
								delete $hash{$chr}{$strand}{$record[$i+1]};
								$record[$i+1]=$record[$i];
							}else {								
								$hash{$chr}{$strand}{$record[$i+1]}[0]=$match2;
								$hash{$chr}{$strand}{$record[$i+1]}[1]=$read2;
								$hash{$chr}{$strand}{$record[$i+1]}[2]=$start2;
								$hash{$chr}{$strand}{$record[$i+1]}[3]+=$hash{$chr}{$strand}{$record[$i]}[3];
								@{$hash{$chr}{$strand}{$record[$i]}}=();
								delete $hash{$chr}{$strand}{$record[$i]};
							}
						}else {
							$hash{$chr}{$strand}{$record[$i+1]}[0]=$match2;
							$hash{$chr}{$strand}{$record[$i+1]}[1]=$read2;
							$hash{$chr}{$strand}{$record[$i+1]}[2]=$start2;
							$hash{$chr}{$strand}{$record[$i+1]}[3]+=$hash{$chr}{$strand}{$record[$i]}[3];
							@{$hash{$chr}{$strand}{$record[$i]}}=();
							delete $hash{$chr}{$strand}{$record[$i]};
						}
					}
				}
			}
		}
	}else {
		my @record=sort {$hash{$chr}{$a}[2]<=>$hash{$chr}{$b}[2]} keys %{$hash{$chr}};
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
				if ($end-$start+1>=$opt{oversize} and (($end-$start+1)/($end1-$start1+1)>=$opt{overrate}/100 or ($end-$start+1)/($end2-$start2+1)>=$opt{overrate}/100)) {
					my $match1=$hash{$chr}{$record[$i]}[0];
					my $read1=$hash{$chr}{$record[$i]}[1];
					my $match2=$hash{$chr}{$record[$i+1]}[0];
					my $read2=$hash{$chr}{$record[$i+1]}[1];
					if ($read1>$read2) {
						$hash{$chr}{$record[$i]}[0]=$match1;
						$hash{$chr}{$record[$i]}[1]=$read1;
						$hash{$chr}{$record[$i]}[2]=$start1;
						$hash{$chr}{$record[$i]}[3]+=$hash{$chr}{$record[$i+1]}[3];
						@{$hash{$chr}{$record[$i+1]}}=();
						delete $hash{$chr}{$record[$i+1]};
						$record[$i+1]=$record[$i];
					}elsif ($read1==$read2) {
						if ($match1>=$match2) {
							$hash{$chr}{$record[$i]}[0]=$match1;
							$hash{$chr}{$record[$i]}[1]=$read1;
							$hash{$chr}{$record[$i]}[2]=$start1;
							$hash{$chr}{$record[$i]}[3]+=$hash{$chr}{$record[$i+1]}[3];
							@{$hash{$chr}{$record[$i+1]}}=();
							delete $hash{$chr}{$record[$i+1]};
							$record[$i+1]=$record[$i];
						}else {								
							$hash{$chr}{$record[$i+1]}[0]=$match2;
							$hash{$chr}{$record[$i+1]}[1]=$read2;
							$hash{$chr}{$record[$i+1]}[2]=$start2;
							$hash{$chr}{$record[$i+1]}[3]+=$hash{$chr}{$record[$i]}[3];
							@{$hash{$chr}{$record[$i]}}=();
							delete $hash{$chr}{$record[$i]};
						}
					}else {
						$hash{$chr}{$record[$i+1]}[0]=$match2;
						$hash{$chr}{$record[$i+1]}[1]=$read2;
						$hash{$chr}{$record[$i+1]}[2]=$start2;
						$hash{$chr}{$record[$i+1]}[3]+=$hash{$chr}{$record[$i]}[3];
						@{$hash{$chr}{$record[$i]}}=();
						delete $hash{$chr}{$record[$i]};
					}
				}
			}
		}
	}
}

$time=localtime;
print "The mapping result has been clustered: $time.\n\n";

#Print cluster result
print "Print the clusters.\n";
open CLUSTER, ">$opt{prefix}\_read.cluster" or die "Cannot write to $opt{prefix}\_read.cluster file.\n";
print CLUSTER "#ID\tQLength\tQStart\tQend\tTStart\tTend\tLength\tScore\tE-value\tOverlap/Total\tIdentity\tSubject_Name\tRead_num\n";
my $id=1;
my $read_num=0;
foreach my $chr (keys %hash) {
	if (defined $opt{strand_specific}) {
		foreach my $strand (keys %{$hash{$chr}}) {
			my @record=sort {$hash{$chr}{$strand}{$a}[2]<=>$hash{$chr}{$strand}{$b}[2]} keys %{$hash{$chr}{$strand}};
			for (my $i=0;$i<@record;$i++) {
				my ($start,$end)=split /\t/,$record[$i];				
				my $match=$hash{$chr}{$strand}{$record[$i]}[0];
				my $read=$hash{$chr}{$strand}{$record[$i]}[3];
				my $len=$end-$start+1;
				my $identity=sprintf "%.0f",$match/$len*100;
				if ($strand eq "-") {
					print CLUSTER "$id\t$len\t1\t$len\t$end\t$start\t$len{$chr}\t\.\t\.\t$match\/$len\t$identity\t$chr\t$read\n";
				}else {
					print CLUSTER "$id\t$len\t1\t$len\t$start\t$end\t$len{$chr}\t\.\t\.\t$match\/$len\t$identity\t$chr\t$read\n";
				}
				$id++;
				$read_num+=$read;
			}
		}
	}else {
		my @record=sort {$hash{$chr}{$a}[2]<=>$hash{$chr}{$b}[2]} keys %{$hash{$chr}};
		for (my $i=0;$i<@record;$i++) {
			my ($start,$end)=split /\t/,$record[$i];
			my $match=$hash{$chr}{$record[$i]}[0];
			my $read=$hash{$chr}{$record[$i]}[3];
			my $len=$end-$start+1;
			my $identity=sprintf "%.0f",$match/$len*100;
			print CLUSTER "$id\t$len\t1\t$len\t$start\t$end\t$len{$chr}\t\.\t\.\t$match\/$len\t$identity\t$chr\t$read\n";
			$id++;
			$read_num+=$read;
		}
	}
}
close CLUSTER;
%hash=();
$id=$id-1;

open (STAT,">>$opt{prefix}.statistics")||die("Cannot write to $opt{prefix}.statistics.\n");
print STAT "Cluster_number\t$id\n";
print STAT "Cluster_read\t$read_num\n\n";
close STAT;

$time=localtime;
print "The clusters have been printed: $time.\n\n";

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

	my $flag=shift;
	my $bin=sprintf("%02b",$flag);
	$bin=sprintf "%011d",$bin;
	$bin=reverse $bin;
	my @bin=split //,$bin;
	my ($pair,$align,$map1,$map2,$strand1,$strand2,$read1,$read2,$multiple,$quality,$duplicate,);
	#if ($bin[0] eq "0") {
	#	$pair=0;
	#}else {
	#	$pair=1;
	#}
	#if ($bin[1] eq "0") {
	#	$align=0;
	#}else {
	#	$align=1;
	#}
	#if ($bin[2] eq "0") {
	#	$map1=1;
	#}else {
	#	$map1=0;
	#}
	#if ($bin[3] eq "0") {
	#	$map2=1;
	#}else {
	#	$map2=0;
	#}
	if ($bin[4] eq "0") {
		$strand1="+";
	}else {
		$strand1="-";
	}
	#if ($bin[5] eq "0") {
	#	$strand2="+";
	#}else {
	#	$strand2="-";
	#}
	#if ($bin[6] eq "0") {
	#	$read1=0;
	#}else {
	#	$read1=1;
	#}
	#if ($bin[7] eq "0") {
	#	$read2=0;
	#}else {
	#	$read2=1;
	#}
	#if ($bin[8] eq "0") {
	#	$multiple=0;
	#}else {
	#	$multiple=1;
	#}
	#if ($bin[9] eq "0") {
	#	$quality=1;
	#}else {
	#	$quality=0;
	#}
	#if ($bin[10] eq "0") {
	#	$duplicate=0;
	#}else {
	#	$duplicate=1;
	#}
	#return $pair,$align,$map1,$map2,$strand1,$strand2,$read1,$read2,$multiple,$quality,$duplicate;
	return $strand1;
}
