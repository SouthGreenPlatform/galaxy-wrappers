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

my $usage=<<USAGE; #******* Instruction of this program *********#

	Author: Wanfei Liu & Chengqi Xin
	Email: <liuwf\@big.ac.cn> & <xinchq\@big.ac.cn>
	Date: Jun 8, 2013
	Version: $version

	Usage: perl expression.pl --statistics <*.statistics> --candidate_fa|-cf <*_candidate_mature.fa|*_candidate_hairpin.fa> --pre_gff|-pg <*_pre.gff> --pre_filter|-pf <*_pre.filter> --read_cluster|-rc <*_read.cluseter_filter|*_read.cluster>
	
	The following options are necessary.
	--statistics	<*.statistics>	:the statistics file produced before
	--candidate_fa|-cf	<*_candidate_mature.fa>	:candidate miRNA fasta file, *_candidate_mature.fa or *_candidate_hairpin.fa
	--pre_gff|-pg	<*_pre.gff>	:candidate precursor gff file in gff2 format
	--pre_filter|-pf	<*_pre.filter>	:filtered precursor file according to precursor RNAfold result
	--read_cluster|-rc	<*_read.cluster>	:read cluster file,  *_read.cluseter_filter or *_read.cluster

	The following options are optional.
	--output_dir|-o		<output_dir>	:default is sample under the current directory, you should change it if you have more than one sample to distinguish different samples
	--prefix|-p		<prefix=MI>	:prefix for file and miRNA ID, if you want to compare different samples, you should use different prefix, the prefix should be end with letters, default is MI
	--strand_specific|-ss			:if you assign this parameter, it means the data is strand specific, if not, as default, it means the data is not strand specific
	--help|-h				:print the usage information

USAGE

#Gather input
&GetOptions(
	"statistics=s"		=>\$opt{statistics},
	"candidate_fa|cf=s"		=>\$opt{candidate_fa},
	"pre_gff|pg=s"		=>\$opt{pre_gff},
	"pre_filter|pf=s"		=>\$opt{pre_filter},
	"read_cluster|rc=s"	=>\$opt{read_cluster},
	"output_dir|o=s"	=>\$opt{output_dir},
	"prefix|p=s"		=>\$opt{prefix},
	"strand_specific|ss"=>\$opt{strand_specific},
	"help|h"	=>\$opt{help},
);

#Verify input
if (!defined $opt{statistics} or !defined $opt{candidate_fa} or !defined $opt{pre_gff} or !defined $opt{pre_filter} or !defined $opt{read_cluster} or defined $opt{help}) {
	die "$usage\n";
}

#Default parameters
$opt{prefix}	= "MI" unless defined $opt{prefix};

# Absolute path
$opt{statistics}  = abs_path $opt{statistics} if defined $opt{statistics};
$opt{candidate_fa}= abs_path $opt{candidate_fa} if defined $opt{candidate_fa};
$opt{pre_gff}     = abs_path $opt{pre_gff} if defined $opt{pre_gff};
$opt{pre_filter}  = abs_path $opt{pre_filter} if defined $opt{pre_filter};
$opt{read_cluster}= abs_path $opt{read_cluster} if defined $opt{read_cluster};

#Output directory
$opt{output_dir} = $ENV{'PWD'}."/sample"  unless defined $opt{output_dir};
$opt{output_dir} = abs_path $opt{output_dir} if defined $opt{output_dir};
`mkdir -p $opt{output_dir}` unless -d $opt{output_dir};
chdir $opt{output_dir};

#Build commond line
my $cmdline = $0;
$cmdline .= " --statistics $opt{statistics}" if defined $opt{statistics};
$cmdline .= " -cf $opt{candidate_fa}" if defined $opt{candidate_fa};
$cmdline .= " -pg $opt{pre_gff}" if defined $opt{pre_gff};
$cmdline .= " -pf $opt{pre_filter}" if defined $opt{pre_filter};
$cmdline .= " -rc $opt{read_cluster}" if defined $opt{read_cluster};
$cmdline .= " -o $opt{output_dir}"  if defined $opt{output_dir};
$cmdline .= " -p $opt{prefix}" if defined $opt{prefix};
$cmdline .= " -ss $opt{strand_specific}" if defined $opt{strand_specific};

print `date`, "\n>>> expression.pl <<<\n\n";
print $cmdline,"\n\n";

my ($time,);

#Read statistics file
print "Read statistics file.\n";
my $total_read=0;
open (IN,"<$opt{statistics}")||die("Cannot read $opt{statistics}.\n");
while (<IN>) {
	chomp;
	if (/Filtered_read\s+(\d+)/) {
		$total_read=$1;
		last;
	}	
}
close IN;

$time=localtime;
print "Statistics file has been readed: $time.\n\n";

my %id=();
#Read candidate miRNA fa file 
print "Read candidate miRNA fa file.\n";
open (IN,"<$opt{candidate_fa}")||die("Cannot read $opt{candidate_fa}.\n");
while (<IN>) {
	chomp;
	if (/^\#/) {
	}elsif (/^>(\S+)/) {
		$id{$1}=1;
	}
}
close IN;

$time=localtime;
print "The candidate miRNA fa file has been readed: $time.\n\n";

print "Read precursor gff file.\n";

my %pos=();
my %hash=();
open (IN,"<$opt{pre_gff}")||die("Cannot read $opt{pre_gff}.\n");
while (<IN>) {
	chomp;
	#S091595 .       miRNA   194     211     .       -       .       ACC="MI1"; ID="ath-miR5021";
	my @list=split /\t/,$_;
	if ($list[8]=~/ACC\=\"([^\"]+)\"/) {
		$pos{$1}=$_ if (defined $id{$1});
	}
}
close IN;

$time=localtime;
print "The precursor gff file has been readed: $time.\n\n";

print "Read $opt{read_cluster} file.\n";

open (IN,"<$opt{read_cluster}")||die("Cannot read $opt{read_cluster}.\n");
while (<IN>) {
	chomp;
	#aau-miR160      21      1       21      1137925 1137905 2539026 34.2    0.15    20/21   95      S000014
	next if (/^\#/);
	my @list=split /\t/,$_;
	if (defined $opt{strand_specific}) {
		if ($list[5] > $list[4]) {
			$hash{$list[11]}{"+"}{"$list[4]\t$list[5]"}[0]=$list[4];
			$hash{$list[11]}{"+"}{"$list[4]\t$list[5]"}[1]=$list[12];
		}else {
			$hash{$list[11]}{"-"}{"$list[5]\t$list[4]"}[0]=$list[5];
			$hash{$list[11]}{"-"}{"$list[5]\t$list[4]"}[1]=$list[12];
		}
	}else {
		if ($list[5] > $list[4]) {
			$hash{$list[11]}{"$list[4]\t$list[5]"}[0]=$list[4];
			$hash{$list[11]}{"$list[4]\t$list[5]"}[1]=$list[12];
		}else {
			$hash{$list[11]}{"$list[5]\t$list[4]"}[0]=$list[5];
			$hash{$list[11]}{"$list[5]\t$list[4]"}[1]=$list[12];
		}
	}
}
close IN;

$time=localtime;
print "The $opt{read_cluster} file has been readed: $time.\n\n";

print "Read $opt{pre_filter} file and calulate the expression for candidate miRNA.\n";

open (OUT,">$opt{prefix}\_candidate.expression")||die("Cannot write to $opt{prefix}\_candidate.expression file.\n");
print OUT "#miRNA\tchr\tstrand\thairpin_start\thairpin_end\thairpin_read\thairpin_express(read_number/total_read*1M)\tmature_start\tmature_end\tmature_read\tmature_express\tstar_start\tstar_end\tstar_read\tstar_express\n";
my %pre=();
open (IN,"<$opt{pre_filter}")||die("fail to open $opt{pre_filter}.\n");
while (<IN>) {
	chomp;
	#miRNA  mfe     pair    unpair  unpair* mat_start       mat_end mat_start*      mat_end*
	#MI304252        -84.20  19      4       6       101     123     64      88
	next if (/^\#/);
	my @list=split /\t/,$_;
	next if (!defined $id{$list[0]});
	if ($list[5]<$list[7]) {
		my @pos=split /\t/,$pos{$list[0]};
		if ($pos[6] eq "+") {
			my $mat_start=$pos[3]+$list[5]-1;
			my $mat_end=$pos[3]+$list[6]-1;
			my $mat_start2=$pos[3]+$list[7]-1;
			my $mat_end2=$pos[3]+$list[8]-1;
			my $pre_start=$pos[3]+$list[5]-1;
			my $pre_end=$pos[3]+$list[8]-1;
			$pre{$pos[0]}{$pos[6]}{$list[0]}[0]=$pre_start;
			$pre{$pos[0]}{$pos[6]}{$list[0]}[1]="$pre_start\t$pre_end\t$mat_start\t$mat_end\t$mat_start2\t$mat_end2";
		}else {
			my $mat_end=$pos[4]-$list[5]+1;
			my $mat_start=$pos[4]-$list[6]+1;
			my $mat_end2=$pos[4]-$list[7]+1;
			my $mat_start2=$pos[4]-$list[8]+1;
			my $pre_end=$pos[4]-$list[5]+1;
			my $pre_start=$pos[4]-$list[8]+1;
			$pre{$pos[0]}{$pos[6]}{$list[0]}[0]=$pre_start;
			$pre{$pos[0]}{$pos[6]}{$list[0]}[1]="$pre_start\t$pre_end\t$mat_start\t$mat_end\t$mat_start2\t$mat_end2";
		}
	}elsif ($list[5]>$list[7]) {
		my @pos=split /\t/,$pos{$list[0]};
		if ($pos[6] eq "+") {
			my $mat_start=$pos[3]+$list[5]-1;
			my $mat_end=$pos[3]+$list[6]-1;
			my $mat_start2=$pos[3]+$list[7]-1;
			my $mat_end2=$pos[3]+$list[8]-1;
			my $pre_start=$pos[3]+$list[7]-1;
			my $pre_end=$pos[3]+$list[6]-1;
			$pre{$pos[0]}{$pos[6]}{$list[0]}[0]=$pre_start;
			$pre{$pos[0]}{$pos[6]}{$list[0]}[1]="$pre_start\t$pre_end\t$mat_start\t$mat_end\t$mat_start2\t$mat_end2";
		}else {
			my $mat_end=$pos[4]-$list[5]+1;
			my $mat_start=$pos[4]-$list[6]+1;
			my $mat_end2=$pos[4]-$list[7]+1;
			my $mat_start2=$pos[4]-$list[8]+1;
			my $pre_end=$pos[4]-$list[7]+1;
			my $pre_start=$pos[4]-$list[6]+1;
			$pre{$pos[0]}{$pos[6]}{$list[0]}[0]=$pre_start;
			$pre{$pos[0]}{$pos[6]}{$list[0]}[1]="$pre_start\t$pre_end\t$mat_start\t$mat_end\t$mat_start2\t$mat_end2";
		}
	}
}
close IN;

my $candidate_read=0;
foreach my $chr (keys %pre) {
	foreach my $strand (keys %{$pre{$chr}}) {
		my @id=sort {$pre{$chr}{$strand}{$a}[0]<=>$pre{$chr}{$strand}{$b}[0]} keys %{$pre{$chr}{$strand}};
		my @hash=();
		if (defined $opt{strand_specific}) {
			@hash=sort {$hash{$chr}{$strand}{$a}[0]<=>$hash{$chr}{$strand}{$b}[0]} keys %{$hash{$chr}{$strand}};
		}else {
			@hash=sort {$hash{$chr}{$a}[0]<=>$hash{$chr}{$b}[0]} keys %{$hash{$chr}};
		}

		my $h=0;
		my $min_overlap_h=undef;
		my ($hstart,$hend)=split /\t/,$hash[$h];
		for (my $i=0;$i<@id;$i++) {
			my ($pstart,$pend,$mstart,$mend,$mstart2,$mend2)=split /\t/,$pre{$chr}{$strand}{$id[$i]}[1];
			my $pre_read=0;
			my $mat_read=0;
			my $mat_read2=0;
			COM:
			if ($hend<$pstart) {
				$h++;
				if (defined $hash[$h]) {
					($hstart,$hend)=split /\t/,$hash[$h];
					goto COM;
				}else {
					printf OUT "$id[$i]\t$chr\t$strand\t$pstart\t$pend\t$pre_read\t%.2f\t$mstart\t$mend\t$mat_read\t%.2f\t$mstart2\t$mend2\t$mat_read2\t%.2f\n",($pre_read/$total_read*1000000),($mat_read/$total_read*1000000),($mat_read2/$total_read*1000000);
					$candidate_read+=$pre_read;
					next;
				}
			}elsif (($hend>=$pstart and $hend<=$pend) or ($pend>=$hstart and $pend<=$hend)) {
				$min_overlap_h=$h if (!defined $min_overlap_h);
				if (defined $opt{strand_specific}) {
					$pre_read+=$hash{$chr}{$strand}{$hash[$h]}[1];
					if (($hend>=$mstart and $hend<=$mend) or ($mend>=$hstart and $mend<=$hend)) {
						$mat_read+=$hash{$chr}{$strand}{$hash[$h]}[1];
					}
					if (($hend>=$mstart2 and $hend<=$mend2) or ($mend2>=$hstart and $mend2<=$hend)) {
						$mat_read2+=$hash{$chr}{$strand}{$hash[$h]}[1];
					}
				}else {
					$pre_read+=$hash{$chr}{$hash[$h]}[1];
					if (($hend>=$mstart and $hend<=$mend) or ($mend>=$hstart and $mend<=$hend)) {
						$mat_read+=$hash{$chr}{$hash[$h]}[1];
					}
					if (($hend>=$mstart2 and $hend<=$mend2) or ($mend2>=$hstart and $mend2<=$hend)) {
						$mat_read2+=$hash{$chr}{$hash[$h]}[1];
					}
				}
				$h++;
				if (defined $hash[$h]) {
					($hstart,$hend)=split /\t/,$hash[$h];
					goto COM;
				}else {
					printf OUT "$id[$i]\t$chr\t$strand\t$pstart\t$pend\t$pre_read\t%.2f\t$mstart\t$mend\t$mat_read\t%.2f\t$mstart2\t$mend2\t$mat_read2\t%.2f\n",($pre_read/$total_read*1000000),($mat_read/$total_read*1000000),($mat_read2/$total_read*1000000);
					$candidate_read+=$pre_read;
					next;
				}
			}elsif ($hstart>$pend) {
				printf OUT "$id[$i]\t$chr\t$strand\t$pstart\t$pend\t$pre_read\t%.2f\t$mstart\t$mend\t$mat_read\t%.2f\t$mstart2\t$mend2\t$mat_read2\t%.2f\n",($pre_read/$total_read*1000000),($mat_read/$total_read*1000000),($mat_read2/$total_read*1000000);
				$candidate_read+=$pre_read;
				$h=$min_overlap_h if (defined $min_overlap_h);
				($hstart,$hend)=split /\t/,$hash[$h];
				$min_overlap_h=undef;
				next;
			}
		}
	}
}
close OUT;

open (STAT,">>$opt{prefix}.statistics")||die("Cannot write to $opt{prefix}.statistics.\n");
print STAT "Candidate_read\t$candidate_read\n\n";
close STAT;

%pos=();
%hash=();
%pre=();

$time=localtime;
print "miRNA expresssion values have been obtained: $time.\n\n";
