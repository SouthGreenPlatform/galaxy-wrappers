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

my $version="1.1.0 version";
my $program=$0;
$program=abs_path $program;
my @programpath = split /\//, $program;
my $programname = pop @programpath;
my $programdir  = join '/', @programpath;

my $usage=<<USAGE; #******* Instruction of this program *********#

	Author: Wanfei Liu & Chengqi Xin
	Email: <liuwf\@big.ac.cn> & <xinchq\@big.ac.cn>
	Date: Jun 18, 2013
	Version: $version

	Introduction: miRNA is a widely known small non-coding RNA which can mediate gene regulation of most important biological processes in plants and animals. Therefore, identification conserve and novel miRNA and their target genes in model and new sequenced species are inevitable. However, the associated tools are often inconvenient, multi-step and difficult to use, especially for biologists who are short for bioinformatics knowledge. MITP is designed to identify miRNA easily and faster based on sequence mapping result from any mapping software which producing sam format output result, blast result (default output result) or blat result (default output result). The program provide a step praramter (8 steps) which can allow running program from any step and finishing all remaining steps. You also can run step by step using each step program. Please run these step programs at the same directory for running main program MITP.pl. Some steps are optional (filter, expression, class and target). When you select the related parameters belong to these optional steps, the program will run these steps. Otherwise, it skip these steps.

	Need to Note: for identification miRNA, the parameter for blast and blat should change comparison to long sequence mapping. In our opinion, the blast command should be as \"blastall -p blastn -d database -i query -v 1000 -b 1000 -W 7 -o outfile\" and the blat command should be as \"blat database query -tileSize=8 -oneOff=1 -minMatch=1 -minIdentity=80 -noHead outfile\".
	
	Usage: perl $program --step <*> { --sam <*> | --blast <*> | --blat <*> } --genome|-g <*.fa>
	
	The following options are necessary.

	--step			:This parameter allow you run program from any step. The details about each step was explained bellow.
		cluster			Run all steps from cluster. In cluster step, the program will cluster reads to obtain candidate regions for downsteam process. If you only want to run this step, please run cluster.pl.
		filter			Run remaining steps from filter. In filter step, the program will filter all regions in filter file (gff2 format) from candiate regions and calculate expression for filtered regions (adding a new attribute EXPRESS=\"Read_number\" at the end of gff2 file named as *.gff.read). The filter file can including any non-intron region from mRNA, rRNA, tRNA, snoRNA and even known miRNA. Please do not include intron region in this file. Otherwise, it will remove all reads in introns. This step is optional. If you only want to run this step, please run filter.pl.
		extract_seq		Run remaining steps from extract_seq. In this step, the program will extract candidate sequence and create *.gff file for candidate sequence. If you only want to run this step, please run extract_seq.pl.
		candidate		Run remaining steps from candidate. In this step, the program can obtain candidate miRNA.
		candidate_filter	Run remaining steps from candidate_filter. In this step, the program will filter hairpin sequences according to their attribute values and create the final sequence and gff files for candidate miRNA.
		expression		Run remaining steps from expression. In this step, the program will calculate the expression for candidate miRNA. If you only want to run this step, please run expression.pl.
		class			Run remaining steps from class. In this step, the program will extract conserve miRNA from candidate miRNA. If you only want to run this step, please run class.pl.
		target			Run the last step which is target. In this step, the program will obtain target genes for candidate miRNA. If you only want to run this step, please run target.pl.

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
	--alignsoft|-as		<blast|blat>	:in default, blast was used.
	--help|-h				:print the usage information

	These parameters belong to cluster step.
	--oversize|-os		<oversize=16>	:minimum overlap size for read cluster, default is 16
	--overrate|-or		<overrate=80>	:minimum overlap rate percent(%) for read cluster, default is 80

	These parameters belong to filter step (This step is optional, it will run only when you assign --filter_gff|-fg or --filter_fa parameter). You can provide a gff file (--filter_gff|-fg) or a fasta file (--filter_fa|-fa) for filter.
	--filter_gff|-fg	<*.gff>		:the filter record file in gff2 format, all clusters overlap with these records will be removed
	--filter_fa|-ff		<*.fa>		:the filter record file in fasta format, all clusters mapped to these sequences will be removed
	--filter_rate|-fr	<filter_rate=50>	:minimum filter overlap rate percent(%) between filter region and read cluster, default is 50

	These parameters belong to extract_seq step.
	--mincov|-mc		<mincov=10>	:minimum miRNA mapping coverage, default is 10 (for conserve miRNA, it should be 1, for novel miRNA from expression, it should be minimum expressed read number)
	
	These parameters belong to candidate step.
	--minbasepair|-mbp	<minbasepair=16>	:minimum base-pairs between mature and star of miRNA comparison, default is 16
	--maxbasebulge|-mbb	<maxbasebulge=3>	:maximum base bulge in mature and star miRNA comparison, default is 3
	--maxunpairbase|-mub	<maxunpairbase=6>	:maximum unpair base number in mature or star miRNA region, default is 6
	--matoverrate|-mor	<matoverrate=60>	:minimum redundant mature sequence overlap rate percent(%), default is 60

	These parameters belong to candidate_filter step. You can filter candidate miRNA according to the candidate sequence and structure attribute values. In default, it include the bellowing parameters.
	--mfei			<mfei=-0.85>	:the mfei means the minimal folding free energy index (MFEI). In default, keep miRNA which mfei value is equal or lower than -0.85
	--mfe			<mfe=-25>	:the mfe means the minimal folding free energy (MFE). In default, keep miRNA which mfe value is equal or lower than -25
	--hairlen		<len=50>	:in default, keep miRNA which hairpin length is equal or larger than 50

	we also provide other filter parameters if you want to further filter candidate miRNA. In default, these parameters are not used

	--pairbase				:minimum pair base number
	--pairpercent				:minimum pair base percent(%)
	--amfe					:the amfe means adjusted mfe (AMFE) represented the mfe of 100 nucleotides. You can set the maximum amfe value 
	--minapercent & --maxapercent		:minimum and maximum A base percent(%)
	--minupercent & --maxapercent		:minimum and maximum U base percent(%)
	--mingpercent & --maxapercent		:minimum and maximum G base percent(%)
	--mincpercent & --maxapercent		:minimum and maximum C base percent(%)
	--minaupercent & --maxapercent		:minimum and maximum AU base percent(%)
	--mingcpercent & --maxapercent		:minimum and maximum GC base percent(%)

	we can do self alignment for candidate miRNA to filter candidate miRNA with same mature and hairpin sequence, but located in different genome region.
	--selfalign|-sa				:if you assign this parameter, it means the program will filter candidate miRNA according to self alignment for mature and hairpin sequence, as default, does not do anything.
	--selfidentity|-si	<selfidentity=100>	:minimum identity percent(%) for filter self alignment result, default is 100
	--selfoverrate|-sor	<overrate=80>	:minimum overlap rate percent(%) for read cluster, default is 80

	we can draw second structure for candidate miRNA.
	--figure		<yes|no>	:default is yes, the program will produce second structure figure in pdf format for all candidate miRNA in a subdirectory named $opt{prefix}\_figure at output directory.

	These parameters belong to expression step. (This step is optional, it will run only when you assign --expression|-e parameter).
	--expression|-e				:if you assign this parameter, it means the program will calculate the expression for candidate miRNA, if not, as default, it means the program will not calculate the expression for candidate miRNA
	
	These parameters belong to class step. (This step is optional, it will run only when you assign --conserveseq|-cs parameter).
	--conserveseq|-cs	<*.fa>		:The mature sequence of conserve miRNA. If you assign this parameter, the program will extract conserve miRNA from candidate miRNA. This file format must like the mature sequence file in miRBase
	--triletterabbr|-tla	<triletterabbr=new>	:three letter abbreviation for studied species, it used to compare miRNA conservation in multiple species, default is new

	These parameters belong to target step. (This step is optional, it will run only when you assign --target_seq|-ts parameter).
	--target_seq|-ts	<*.fa>		:The target sequence file. If you assign this parameter, the program will do target prediction for candidate miRNA.
	--target_tool|-tt	<*>		:The program for target prediction (blast, blat, targetfinder, miranda and RNAhybrid). In default, it will run using the program assigned by --alignsoft|-as parameter. For plant, you can identify target genes by our rule set according to blast or blat alignment (in default) or by targetfinder; for animal, you can use miranda or RNAhybrid.
	--dataset|-ds	<dataset=3utr_fly|3utr_worm|3utr_human>	:three data set name in RNAhybird program for target prediction, in default is 3utr_human. This parameter is only need when you use RNAhybrid to predict target

	Note: for conserve miRNA identification, --mincov|-mc should be 1; for novel miRNA identification, --mincov|-mc should be the minimum read number for expression. For --filter_gff|-fg or --filter_fa|-ff parameter, please only inculde region or sequence you want to removed (for example, do not include intron in animal).

USAGE

#Gather input
&GetOptions(
	"step=s"			=>\$opt{step},
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
	"alignsoft|as=s"	=>\$opt{alignsoft},

	"oversize|os=i"		=>\$opt{oversize},
	"overrate|or=i"		=>\$opt{overrate},

	"filter_gff|fg=s"	=>\$opt{filter_gff},
	"filter_fa|ff=s"	=>\$opt{filter_fa},
	"filter_rate|fr=i"	=>\$opt{filter_rate},
	
	"mincov|mc=s"		=>\$opt{mincov},

	"minbasepair|mbp=i"	=>\$opt{minbasepair},
	"maxbasebulge|mbb=i"=>\$opt{maxbasebulge},
	"maxunpairbase|mub=i"=>\$opt{maxunpairbase},
	"matoverrate|mor=i"	=>\$opt{matoverrate},

	"mfei=i"			=>\$opt{mfei},
	"mfe=i"				=>\$opt{mfe},
	"hairlen=i"			=>\$opt{hairlen},	
	"pairbase=i"		=>\$opt{pairbase},
	"pairpercent=i"		=>\$opt{pairpercent},
	"amfe=i"			=>\$opt{amfe},
	"minapercent=i"		=>\$opt{minapercent},
	"maxapercent=i"		=>\$opt{maxapercent},
	"minupercent=i"		=>\$opt{minupercent},
	"maxupercent=i"		=>\$opt{maxupercent},
	"mingpercent=i"		=>\$opt{mingpercent},
	"maxgpercent=i"		=>\$opt{maxgpercent},
	"mincpercent=i"		=>\$opt{mincpercent},
	"maxcpercent=i"		=>\$opt{maxcpercent},
	"minaupercent=i"	=>\$opt{minaupercent},
	"maxaupercent=i"	=>\$opt{maxaupercent},
	"mingcpercent=i"	=>\$opt{mingcpercent},
	"maxgcpercent=i"	=>\$opt{maxgcpercent},
	"selfalign|sa"		=>\$opt{selfalign},
	"selfidentity|si=i"	=>\$opt{selfidentity},
	"selfoverrate|sor=i"=>\$opt{selfoverrate},
	"figure=s"			=>\$opt{figure},

	"expression|e"	=>\$opt{expression},

	"conserveseq|cs=s"	=>\$opt{conserveseq},
	"triletterabbr|tla=s"	=>\$opt{triletterabbr},

	"target_seq|ts=s"	=>\$opt{target_seq},
	"target_tool|tt=s"	=>\$opt{target_tool},
	"dataset|ds=s"	=>\$opt{dataset},
	
	"help|h"	=>\$opt{help},
);

#Verify input
if ((!defined $opt{sam} and !defined $opt{blast} and !defined $opt{blat}) or !defined $opt{genome} or defined $opt{help}) {
	die "$usage\n";
}

#Default parameters
$opt{step}		="cluster" unless defined $opt{step};
$opt{prefix}	= "MI" unless defined $opt{prefix};
$opt{maxmap}	= 10 unless defined $opt{maxmap};
$opt{mismatch}	= 3 unless defined $opt{mismatch};
$opt{identity}	= 90 unless defined $opt{identity};
$opt{minimum}	= 18 unless defined $opt{minimum};
$opt{maximum}	= 25 unless defined $opt{maximum};
$opt{alignsoft}	= "blast" unless defined $opt{alignsoft};
$opt{oversize}	= 16 unless defined $opt{oversize};
$opt{overrate}	= 80 unless defined $opt{overrate};
$opt{filter_rate}	= 50 unless defined $opt{filter_rate};
$opt{mincov}	= 10 unless defined $opt{mincov};
$opt{minbasepair}	= 16 unless defined $opt{minbasepair};
$opt{maxbasebulge}	= 3 unless defined $opt{maxbasebulge};
$opt{maxunpairbase}	= 6 unless defined $opt{maxunpairbase};
$opt{matoverrate}	= 60 unless defined $opt{matoverrate};
$opt{mfei}	= -0.85 unless defined $opt{mfei};
$opt{mfe}	= -25 unless defined $opt{mfe};
$opt{hairlen}	= 50 unless defined $opt{hairlen};
$opt{selfidentity}	= 100 unless defined $opt{selfidentity};
$opt{selfoverrate}	= 80 unless defined $opt{selfoverrate};
$opt{figure}	= "yes" unless defined $opt{figure};
$opt{triletterabbr}	= "new" unless defined $opt{triletterabbr};
$opt{target_tool}	= $opt{alignsoft} unless defined $opt{target_tool};
$opt{dataset}	= "3utr_human" unless defined $opt{dataset};

# Absolute path
$opt{sam}         = abs_path $opt{sam} if defined $opt{sam};
$opt{blast}       = abs_path $opt{blast} if defined $opt{blast};
$opt{blat}        = abs_path $opt{blat} if defined $opt{blat};
$opt{genome}      = abs_path $opt{genome};
$opt{filter_gff}  = abs_path $opt{filter_gff} if defined $opt{filter_gff};
$opt{filter_fa}   = abs_path $opt{filter_fa} if defined $opt{filter_fa};
$opt{conserveseq} = abs_path $opt{conserveseq} if defined $opt{conserveseq};
$opt{target_seq}  = abs_path $opt{target_seq} if defined $opt{target_seq};

#Output directory
$opt{output_dir} = $ENV{'PWD'}."/sample"  unless defined $opt{output_dir};
$opt{output_dir} = abs_path $opt{output_dir} if defined $opt{output_dir};
`mkdir -p $opt{output_dir}` unless -d $opt{output_dir};
chdir $opt{output_dir};

#Build commond line
my $cmdline = $program;
$cmdline .= " --step $opt{step}" if defined $opt{step};
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
$cmdline .= " -as $opt{alignsoft}" if defined $opt{alignsoft};
$cmdline .= " -os $opt{oversize}" if defined $opt{oversize};
$cmdline .= " -or $opt{overrate}" if defined $opt{overrate};
$cmdline .= " -fg $opt{filter_gff}"  if defined $opt{filter_gff};
$cmdline .= " -ff $opt{filter_fa}"  if defined $opt{filter_fa};
$cmdline .= " -fr $opt{filter_rate}"  if defined $opt{filter_rate};
$cmdline .= " -mc $opt{mincov}" if defined $opt{mincov};
$cmdline .= " -mbp $opt{minbasepair}" if defined $opt{minbasepair};
$cmdline .= " -mbb $opt{maxbasebulge}" if defined $opt{maxbasebulge};
$cmdline .= " -mub $opt{maxunpairbase}" if defined $opt{maxunpairbase};
$cmdline .= " -mor $opt{matoverrate}" if defined $opt{matoverrate};
$cmdline .= " --mfei $opt{mfei}" if defined $opt{mfei};
$cmdline .= " --mfe $opt{mfe}" if defined $opt{mfe};
$cmdline .= " --hairlen $opt{hairlen}" if defined $opt{hairlen};
$cmdline .= " --pairbase $opt{pairbase}" if defined $opt{pairbase};
$cmdline .= " --pairpercent $opt{pairpercent}" if defined $opt{pairpercent};
$cmdline .= " --amfe $opt{amfe}" if defined $opt{amfe};
$cmdline .= " --minapercent $opt{minapercent}" if defined $opt{minapercent};
$cmdline .= " --maxapercent $opt{maxapercent}" if defined $opt{maxapercent};
$cmdline .= " --minupercent $opt{minupercent}" if defined $opt{minupercent};
$cmdline .= " --maxupercent $opt{maxupercent}" if defined $opt{maxupercent};
$cmdline .= " --mingpercent $opt{mingpercent}" if defined $opt{mingpercent};
$cmdline .= " --maxgpercent $opt{maxgpercent}" if defined $opt{maxgpercent};
$cmdline .= " --mincpercent $opt{mincpercent}" if defined $opt{mincpercent};
$cmdline .= " --maxcpercent $opt{maxcpercent}" if defined $opt{maxcpercent};
$cmdline .= " --minaupercent $opt{minaupercent}" if defined $opt{minaupercent};
$cmdline .= " --maxaupercent $opt{maxaupercent}" if defined $opt{maxaupercent};
$cmdline .= " --mingcpercent $opt{mingcpercent}" if defined $opt{mingcpercent};
$cmdline .= " --maxgcpercent $opt{maxgcpercent}" if defined $opt{maxgcpercent};
$cmdline .= " -sa $opt{selfalign}" if defined $opt{selfalign};
$cmdline .= " -si $opt{selfidentity}" if defined $opt{selfidentity};
$cmdline .= " -sor $opt{selfoverrate}" if defined $opt{selfoverrate};
$cmdline .= " --figure $opt{figure}" if defined $opt{figure};
$cmdline .= " -e $opt{expression}" if defined $opt{expression};
$cmdline .= " -cs $opt{conserveseq}" if defined $opt{conserveseq};
$cmdline .= " -tla $opt{triletterabbr}" if defined $opt{triletterabbr};
$cmdline .= " -ts $opt{target_seq}" if defined $opt{target_seq};
$cmdline .= " -tt $opt{target_tool}" if defined $opt{target_tool};
$cmdline .= " -ds $opt{dataset}" if defined $opt{dataset};

open LOG, ">log" or die "Cannot write to log file\n";
print LOG `date`, "\n>>> $program <<<\n\n";
print LOG $cmdline,"\n\n";

my ($time,$chr,$string,%gen,%len,%id,%hash,);
if ($opt{step} eq "cluster" or $opt{step} eq "extract_seq" or defined $opt{filter_fa}) {
	#Obtain chromosome sequence and length
	print LOG "Obtain chromosome sequence and length.\n";
	
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
	print LOG "Chromosome sequence and length are obtained: $time.\n\n";
}

%hash=();
if ($opt{step} eq "cluster") {
	#Read the mapping result file
	print LOG "Read the mapping result file.\n";
	
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
			next if ($list[5]=~/S|H/);
			if (!defined $id{$list[0]}) {
				$id{$list[0]}=0;
			}else {
				if ($id{$list[0]}>=$opt{maxmap}) {
					next;
				}
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

		print STAT "Total_read\t$total_read\nFiltered_read\t$filtered_read\n\n";
		close STAT;
	}

	$time=localtime;
	print LOG "The mapping result has been readed: $time.\n\n";

	#Cluster mapping result
	print LOG "Cluster mapping result.\n";
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
	print LOG "The mapping result has been clustered: $time.\n\n";

	#Print cluster result
	print LOG "Print the clusters.\n";
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
	print LOG "The clusters have been printed: $time.\n\n";
}

%hash=();

if ($opt{step} eq "cluster" or $opt{step} eq "filter") {
	if (defined $opt{filter_gff}) {
		#Read $opt{prefix}_read.cluster file.\n";
		print LOG "Read $opt{prefix}\_read.cluster file.\n";

		open (IN,"<$opt{prefix}\_read.cluster")||die("Cannot read $opt{prefix}\_read.cluster.\n");
		while (<IN>) {
			chomp;
			#aau-miR160      21      1       21      1137925 1137905 2539026 34.2    0.15    20/21   95      S000014
			next if (/^\#/);
			my @list=split /\t/,$_;
			if (defined $opt{strand_specific}) {
				if ($list[5] > $list[4]) {
					$hash{$list[11]}{"+"}{"$list[4]\t$list[5]"}[0]=$list[4];
					$hash{$list[11]}{"+"}{"$list[4]\t$list[5]"}[1]=$list[12];
					$hash{$list[11]}{"+"}{"$list[4]\t$list[5]"}[2]=$_;
				}else {
					$hash{$list[11]}{"-"}{"$list[5]\t$list[4]"}[0]=$list[5];
					$hash{$list[11]}{"-"}{"$list[5]\t$list[4]"}[1]=$list[12];
					$hash{$list[11]}{"-"}{"$list[5]\t$list[4]"}[2]=$_;
				}
			}else {
				if ($list[5] > $list[4]) {
					$hash{$list[11]}{"$list[4]\t$list[5]"}[0]=$list[4];
					$hash{$list[11]}{"$list[4]\t$list[5]"}[1]=$list[12];
					$hash{$list[11]}{"$list[4]\t$list[5]"}[2]=$_;
				}else {
					$hash{$list[11]}{"$list[5]\t$list[4]"}[0]=$list[5];
					$hash{$list[11]}{"$list[5]\t$list[4]"}[1]=$list[12];
					$hash{$list[11]}{"$list[5]\t$list[4]"}[2]=$_;
				}
			}
		}
		close IN;

		$time=localtime;
		print LOG "The $opt{prefix}_read.cluster file has been readed: $time.\n\n";

		#Filter clusters which belong to regions in filter file (gff2 format)
		print LOG "Filter clusters which belong to regions in filter file.\n";

		my %filter=();
		my $connector=undef;
		open (IN,"<$opt{filter_gff}")||die("Cannot read $opt{filter_gff}.\n");
		while (<IN>) {
			chomp;
			#S091595 .       miRNA   194     211     .       -       .       ACC="MI1"; ID="ath-miR5021";
			my @list=split /\t/,$_;
			if (!defined $connector) {
				if ($list[8]=~/\=/) {
					$connector="=";
				}else {
					$connector=" ";
				}
			}
			$filter{$list[0]}{$list[6]}{"$list[3]\t$list[4]"}[0]=$list[3];
			$filter{$list[0]}{$list[6]}{"$list[3]\t$list[4]"}[1]=$_;
		}
		close IN;

		my $out=basename($opt{filter_gff}).".read";

		open FILTER, ">$out" or die "Cannot write to $out file.\n";
		foreach my $chr (keys %filter) {
			foreach my $strand (keys %{$filter{$chr}}) {
				my @filter=sort {$filter{$chr}{$strand}{$a}[0]<=>$filter{$chr}{$strand}{$b}[0]} keys %{$filter{$chr}{$strand}};
				my @hash=();
				if (defined $opt{strand_specific}) {
					@hash=sort {$hash{$chr}{$strand}{$a}[0]<=>$hash{$chr}{$strand}{$b}[0]} keys %{$hash{$chr}{$strand}};
				}else {
					@hash=sort {$hash{$chr}{$a}[0]<=>$hash{$chr}{$b}[0]} keys %{$hash{$chr}};
				}

				my $h=0;
				my $min_overlap_h=undef;
				my ($hstart,$hend)=split /\t/,$hash[$h];
				for (my $f=0;$f<@filter;$f++) {
					my ($fstart,$fend)=split /\t/,$filter[$f];
					my $express=0;
					FCOM:
					if ($hend<$fstart) {
						$h++;
						if (defined $hash[$h]) {
							($hstart,$hend)=split /\t/,$hash[$h];
							goto FCOM;
						}else {
							print FILTER "$filter{$chr}{$strand}{$filter[$f]}[1] EXPRESS$connector\"$express\"\n";
							next;
						}
					}elsif (($hend>=$fstart and $hend<=$fend) or ($fend>=$hstart and $fend<=$hend)) {
						$min_overlap_h=$h if (!defined $min_overlap_h);
						my $start=$hstart;
						my $end=$hend;
						if ($start<$fstart) {
							$start=$fstart;
						}
						if ($end>$fend) {
							$end=$fend;
						}
						if (($end-$start+1)/($hend-$hstart+1)>=$opt{filter_rate}/100 or ($end-$start+1)/($fend-$fstart+1)>=$opt{filter_rate}/100) {
							if (defined $opt{strand_specific}) {
								$express+=$hash{$chr}{$strand}{$hash[$h]}[1];
								$hash{$chr}{$strand}{$hash[$h]}[3]=1;
							}else {
								$express+=$hash{$chr}{$hash[$h]}[1];
								$hash{$chr}{$hash[$h]}[3]=1;
							}
						}
						$h++;
						if (defined $hash[$h]) {
							($hstart,$hend)=split /\t/,$hash[$h];
							goto FCOM;
						}else {
							print FILTER "$filter{$chr}{$strand}{$filter[$f]}[1] EXPRESS$connector\"$express\"\n";
							next;
						}
					}elsif ($hstart>$fend) {
						print FILTER "$filter{$chr}{$strand}{$filter[$f]}[1] EXPRESS$connector\"$express\"\n";
						$h=$min_overlap_h if (defined $min_overlap_h);
						($hstart,$hend)=split /\t/,$hash[$h];
						$min_overlap_h=undef;
						next;
					}
				}
			}
		}
		close FILTER;

		$time=localtime;
		print LOG "The clusters which belong to regions in filter file have been filtered: $time.\n\n";

		#Print filtered cluster result
		print LOG "Print filtered cluster result.\n";
		open CLUSTER, ">$opt{prefix}\_read.cluster_filter" or die "Cannot write to $opt{prefix}\_read.cluster_filter file.\n";
		print CLUSTER "#ID\tQLength\tQStart\tQend\tTStart\tTend\tLength\tScore\tE-value\tOverlap/Total\tIdentity\tSubject_Name\tRead_num\n";
		my $filtered_cluster=0;
		my $filtered_read=0;
		foreach my $chr (keys %hash) {
			if (defined $opt{strand_specific}) {
				foreach my $strand (keys %{$hash{$chr}}) {
					my @record=sort {$hash{$chr}{$strand}{$a}[0]<=>$hash{$chr}{$strand}{$b}[0]} keys %{$hash{$chr}{$strand}};
					for (my $i=0;$i<@record;$i++) {
						if (!defined $hash{$chr}{$strand}{$record[$i]}[3]) {
							print CLUSTER "$hash{$chr}{$strand}{$record[$i]}[2]\n";
							$filtered_cluster++;
							$filtered_read+=$hash{$chr}{$strand}{$record[$i]}[1];
						}
						@{$hash{$chr}{$strand}{$record[$i]}}=();
						delete $hash{$chr}{$strand}{$record[$i]};
					}
				}
			}else {
				my @record=sort {$hash{$chr}{$a}[0]<=>$hash{$chr}{$b}[0]} keys %{$hash{$chr}};
				for (my $i=0;$i<@record;$i++) {
					if (!defined $hash{$chr}{$record[$i]}[3]) {
						print CLUSTER "$hash{$chr}{$record[$i]}[2]\n";
						$filtered_cluster++;
						$filtered_read+=$hash{$chr}{$record[$i]}[1];
					}
					@{$hash{$chr}{$record[$i]}}=();
					delete $hash{$chr}{$record[$i]};
				}
			}
		}
		close CLUSTER;

		open (STAT,">>$opt{prefix}.statistics")||die("Cannot write to $opt{prefix}.statistics.\n");
		print STAT "Cluster_number_after_filter\t$filtered_cluster\n";
		print STAT "Cluster_read_after_filter\t$filtered_read\n\n";
		close STAT;

		%hash=();

		$time=localtime;
		print LOG "The filtered clusters result have been printed: $time.\n\n";

	}elsif (defined $opt{filter_fa}) {

		#Read $opt{prefix}_read.cluster file.\n";
		print LOG "Read $opt{prefix}\_read.cluster file.\n";

		open (OUT,">$opt{prefix}\_read.cluster.fa")||die("Cannot write to $opt{prefix}\_read.cluster.fa file.\n");
		open (IN,"<$opt{prefix}_read.cluster")||die("Cannot read $opt{prefix}_read.cluster.\n");
		while (<IN>) {
			chomp;
			#aau-miR160      21      1       21      1137925 1137905 2539026 34.2    0.15    20/21   95      S000014
			next if (/^\#/);
			my @list=split /\t/,$_;
			if ($list[5] > $list[4]) {
				my $cluster_start=$list[4];
				my $cluster_end=$list[5];
				my $cluster_seq=substr ($gen{$list[11]},$cluster_start-1,abs($cluster_end-$cluster_start)+1);
				print OUT ">$list[0]\n$cluster_seq\n";
				$hash{$list[0]}=$_;
			}else {
				my $cluster_start=$list[5];
				my $cluster_end=$list[4];
				my $cluster_seq=substr ($gen{$list[11]},$cluster_start-1,abs($cluster_end-$cluster_start)+1);
				$cluster_seq=reverse $cluster_seq;
				$cluster_seq=~tr/ACGTacgt/TGCAtgca/;
				print OUT ">$list[0]\n$cluster_seq\n";
				$hash{$list[0]}=$_;
			}
		}
		close IN;
		close OUT;

		$time=localtime;
		print LOG "The $opt{prefix}_read.cluster file has been readed: $time.\n\n";

		#Filter clusters which belong to regions in filter file (gff2 format)
		print LOG "Filter clusters which belong to regions in filter file.\n";

		if ($opt{alignsoft} eq "blast") {
			system "formatdb -i $opt{filter_fa} -p F";
			system "blastall -p blastn -d $opt{filter_fa} -i $opt{prefix}\_read.cluster.fa -v 1000 -b 1000 -W 7 -o $opt{prefix}\_read.cluster.blast";
			system "perl $programdir/EblastN.pl -i $opt{prefix}\_read.cluster.blast -o $opt{prefix}\_read.cluster.eblastn";

			my %filter=();
			open (IN,"<$opt{prefix}\_read.cluster.eblastn")||die("Cannot read $opt{prefix}\_read.cluster.eblastn.\n");
			while (<IN>) {
				next if (/^Query/);
				my @list=split(/\t/,$_);
				my ($overlap,$total)=split /\//,$list[9];
				my $mis=$total-$overlap;
				if ($list[10]>=$opt{identity} and $mis<=$opt{mismatch} and (($list[3]-$list[2]+1)/$list[1]>=$opt{filter_rate}/100 or (abs($list[5]-$list[4])+1)/$list[6]>=$opt{filter_rate}/100)) {
					my @cluster=split /\t/,$hash{$list[0]};
					$filter{$list[11]}+=$cluster[$#cluster];
					delete $hash{$list[0]};
				}
			}
			close IN;

			my $out=basename($opt{filter_fa})."read";

			open FILTER, ">$out" or die "Cannot write to $out file.\n";
			open IN, "<$opt{filter_fa}" or die "Cannot read $opt{filter_fa} file.\n";
			while (<IN>) {
				if (/^>(\S+)/) {
					if (defined $filter{$1}) {
						print FILTER "$1\t$filter{$1}\n";
					}else {
						print FILTER "$1\t0\n";
					}
				}
			}
			close IN;
			close FILTER;

			$time=localtime;
			print LOG "The clusters which belong to regions in filter file have been filtered: $time.\n\n";

			#Print filtered cluster result
			print LOG "Print filtered cluster result.\n";
			open CLUSTER, ">$opt{prefix}\_read.cluster_filter" or die "Cannot write to $opt{prefix}\_read.cluster_filter file.\n";
			print CLUSTER "#ID\tQLength\tQStart\tQend\tTStart\tTend\tLength\tScore\tE-value\tOverlap/Total\tIdentity\tSubject_Name\tRead_num\n";
			my $filtered_cluster=0;
			my $filtered_read=0;
			foreach my $id (sort {$a<=>$b} keys %hash) {
				print CLUSTER "$hash{$id}\n";
				$filtered_cluster++;
				my @list=split /\t/,$hash{$id};
				$filtered_read+=$list[$#list];
			}
			close CLUSTER;

			open (STAT,">>$opt{prefix}.statistics")||die("Cannot write to $opt{prefix}.statistics.\n");
			print STAT "Cluster_number_after_filter\t$filtered_cluster\n";
			print STAT "Cluster_read_after_filter\t$filtered_read\n\n";
			close STAT;

			%hash=();

			$time=localtime;
			print LOG "The filtered clusters result have been printed: $time.\n\n";

		}elsif ($opt{alignsoft} eq "blat") {
			system "blat $opt{filter_fa} $opt{prefix}\_read.cluster.fa -tileSize=8 -oneOff=1 -minMatch=1 -minIdentity=80 -noHead $opt{prefix}\_read.cluster.blat";
			
			my %filter=();
			open (IN,"<$opt{prefix}\_read.cluster.blat")||die("Cannot read $opt{prefix}\_read.cluster.blat.\n");
			while (<IN>) {
				chomp;
				next if (/^[\D]+/);
				my @list=split /\t/,$_;
				my $identity=$list[0]/($list[0]+$list[1])*100;
				my $mis=$list[1]+$list[5]+$list[7];
				if ($identity>=$opt{identity} and $mis<=$opt{mismatch} and (($list[12]-$list[11])/$list[10]>=$opt{filter_rate}/100 or ($list[16]-$list[15])/$list[14]>=$opt{filter_rate}/100)) {
					my @cluster=split /\t/,$hash{$list[9]};
					$filter{$list[13]}+=$cluster[$#cluster];
					delete $hash{$list[9]};
				}
			}
			close IN;

			my $out=basename($opt{filter_fa}).".read";

			open FILTER, ">$out" or die "Cannot write to $out file.\n";
			open IN, ">$opt{filter_fa}" or die "Cannot read $opt{filter_fa} file.\n";
			while (<IN>) {
				if (/^>(\S+)/) {
					if (defined $filter{$1}) {
						print FILTER "$1\t$filter{$1}\n";
					}else {
						print FILTER "$1\t0\n";
					}
				}
			}
			close IN;
			close FILTER;

			$time=localtime;
			print LOG "The clusters which belong to regions in filter file have been filtered: $time.\n\n";

			#Print filtered cluster result
			print LOG "Print filtered cluster result.\n";
			open CLUSTER, ">$opt{prefix}\_read.cluster_filter" or die "Cannot write to $opt{prefix}\_read.cluster_filter file.\n";
			print CLUSTER "#ID\tQLength\tQStart\tQend\tTStart\tTend\tLength\tScore\tE-value\tOverlap/Total\tIdentity\tSubject_Name\tRead_num\n";
			my $filtered_cluster=0;
			my $filtered_read=0;
			foreach my $id (sort {$a<=>$b} keys %hash) {
				print CLUSTER "$hash{$id}\n";
				$filtered_cluster++;
				my @list=split /\t/,$hash{$id};
				$filtered_read+=$list[$#list];
			}
			close CLUSTER;

			open (STAT,">>$opt{prefix}.statistics")||die("Cannot write to $opt{prefix}.statistics.\n");
			print STAT "Cluster_number_after_filter\t$filtered_cluster\n";
			print STAT "Cluster_read_after_filter\t$filtered_read\n\n";
			close STAT;

			%hash=();

			$time=localtime;
			print LOG "The filtered clusters result have been printed: $time.\n\n";
		}
	}else {
		print LOG "Does not assign filter file.\n";
		system "cp $opt{prefix}_read.cluster $opt{prefix}_read.cluster_filter";
		$time=localtime;
		print LOG "Just do a copy named $opt{prefix}_read.cluster_filter for $opt{prefix}_read.cluster file: $time.\n\n";

		open (STAT,">>$opt{prefix}.statistics")||die("Cannot write to $opt{prefix}.statistics.\n");
		print STAT "Filtered_cluster_number\t0\n";
		print STAT "Filtered_cluster_read\t0\n\n";
		close STAT;
	}
}

if ($opt{step} eq "cluster" or $opt{step} eq "filter" or $opt{step} eq "extract_seq") {

	#Creat mature and precursor sequence and gff file for candidate sequence
	print LOG "Create mature and precursor sequence and gff files.\n";
	open (OUT1,">$opt{prefix}\_pre.fa")||die("Cannot write to $opt{prefix}\_pre.fa file.\n");
	open (OUT2,">$opt{prefix}\_mature.fa")||die("Cannot write to $opt{prefix}\_mature.fa file.\n");
	open (OUT3,">$opt{prefix}\_pre.gff")||die("Cannot write to $opt{prefix}\_pre.gff file.\n");
	open (OUT4,">$opt{prefix}\_mature.gff")||die("Cannot write to $opt{prefix}\_mature.gff file.\n");

	my $num=1;
	open (IN,"<$opt{prefix}_read.cluster_filter")||die("Cannot read $opt{prefix}_read.cluster_filter.\n");
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
			#filter low complexity miRNA
			my $type1=&low_complexity ($mat_seq);
			#my $type2=&low_complexity ($pre_seq);
			next if ($type1==0);
			if (defined $opt{strand_specific}) {
				print OUT1 ">$opt{prefix}$num $list[0]\n$pre_seq\n";
				print OUT2 ">$opt{prefix}$num $list[0]\n$mat_seq\n";
				print OUT3 "$list[11]\tMITP\tmiRNA\t$pre_start\t$pre_end\t\.\t+\t\.\tACC=\"$opt{prefix}$num\"\; ID=\"$list[0],+\"\; EXPRESS=\"$list[$#list]\"\;\n";
				print OUT4 "$list[11]\tMITP\tmiRNA\t$mat_start\t$mat_end\t\.\t+\t\.\tACC=\"$opt{prefix}$num\"\; ID=\"$list[0],+\"\; EXPRESS=\"$list[$#list]\"\;\n";
				$num++;
			}else {
				print OUT1 ">$opt{prefix}$num $list[0]\n$pre_seq\n";
				print OUT2 ">$opt{prefix}$num $list[0]\n$mat_seq\n";
				print OUT3 "$list[11]\tMITP\tmiRNA\t$pre_start\t$pre_end\t\.\t+\t\.\tACC=\"$opt{prefix}$num\"\; ID=\"$list[0],+\"\; EXPRESS=\"$list[$#list]\"\;\n";
				print OUT4 "$list[11]\tMITP\tmiRNA\t$mat_start\t$mat_end\t\.\t+\t\.\tACC=\"$opt{prefix}$num\"\; ID=\"$list[0],+\"\; EXPRESS=\"$list[$#list]\"\;\n";
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
			#filter low complexity miRNA
			my $type1=&low_complexity ($mat_seq);
			#my $type2=&low_complexity ($pre_seq);
			next if ($type1==0);
			if (defined $opt{strand_specific}) {
				$mat_seq=reverse $mat_seq;
				$mat_seq=~tr/ACGTacgt/TGCAtgca/;
				$pre_seq=reverse $pre_seq;
				$pre_seq=~tr/ACGTacgt/TGCAtgca/;
				print OUT1 ">$opt{prefix}$num $list[0]\n$pre_seq\n";
				print OUT2 ">$opt{prefix}$num $list[0]\n$mat_seq\n";
				print OUT3 "$list[11]\tMITP\tmiRNA\t$pre_start\t$pre_end\t\.\t-\t\.\tACC=\"$opt{prefix}$num\"\; ID=\"$list[0],-\"\; EXPRESS=\"$list[$#list]\"\;\n";
				print OUT4 "$list[11]\tMITP\tmiRNA\t$mat_start\t$mat_end\t\.\t-\t\.\tACC=\"$opt{prefix}$num\"\; ID=\"$list[0],-\"\; EXPRESS=\"$list[$#list]\"\;\n";
				$num++;
			}else {
				$mat_seq=reverse $mat_seq;
				$mat_seq=~tr/ACGTacgt/TGCAtgca/;
				$pre_seq=reverse $pre_seq;
				$pre_seq=~tr/ACGTacgt/TGCAtgca/;
				print OUT1 ">$opt{prefix}$num $list[0]\n$pre_seq\n";
				print OUT2 ">$opt{prefix}$num $list[0]\n$mat_seq\n";
				print OUT3 "$list[11]\tMITP\tmiRNA\t$pre_start\t$pre_end\t\.\t-\t\.\tACC=\"$opt{prefix}$num\"\; ID=\"$list[0],-\"\; EXPRESS=\"$list[$#list]\"\;\n";
				print OUT4 "$list[11]\tMITP\tmiRNA\t$mat_start\t$mat_end\t\.\t-\t\.\tACC=\"$opt{prefix}$num\"\; ID=\"$list[0],-\"\; EXPRESS=\"$list[$#list]\"\;\n";
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
	print LOG "The mature and precursor sequence and gff files have been created: $time.\n\n";
}

if ($opt{step} eq "cluster" or $opt{step} eq "filter" or $opt{step} eq "extract_seq" or $opt{step} eq "candidate") {

	#Do RNAfold for obtained precursor sequence
	print LOG "Do RNAfold for obtained precursor sequence.\n";
	system "less $opt{prefix}\_pre.fa | RNAfold --noconv --noPS > $opt{prefix}\_pre.RNAfold";

	$time=localtime;
	print LOG "The RNAfold for obtained precursor sequence has been done: $time.\n\n";

	#Filter precursor based on RNAfold result
	print LOG "Filter precursor based on RNAfold result.\n";
	my ($key,%mat,);
	open (MAT,"<$opt{prefix}\_mature.fa")||die("Cannot read $opt{prefix}\_mature.fa.\n");
	while (<MAT>) {
		chomp;
		if (/^>(\S+)/) {
			$key=$1;
		}else{
			$mat{$key}.=$_;
		}
	}
	close MAT;

	open (FILTER,">$opt{prefix}\_pre.filter")||die("fail to open $opt{prefix}\_pre.filter.\n");
	print FILTER "#miRNA\tMFE\tPair\tMature_unpair\tStar_unpair\tMature_start\tMature_end\tStar_start\tStar_end\n";
	open (ALL,">$opt{prefix}\_pre.all")||die("fail to open $opt{prefix}\_pre.all.\n");
	print ALL "#miRNA\tMFE\tPair\tMature_unpair\tStar_unpair\tMature_start\tMature_end\tStar_start\tStar_end\n";

	open (IN,"<$opt{prefix}\_pre.RNAfold")||die("Cannot read $opt{prefix}\_pre.RNAfold.\n");
	my $seq_num=0;
	while (!eof(IN)) {
		#>MI1 1
		#AAACCAGATTACAGGTCCAGATCAGATCTCGGGCCCAAACCAAACCAGCCCACAGATCCAATCAGACCCAACTCAGCCCAAGAACAGTTCCCTAGCCCAAACCAGATCTGATCAGACCAGATCAGATCAAATC
		#.....(((((...((((..(((..(((((.((((.............))))..)))))..))).)))).......((....(((...)))....))........(((((((((......))))))))).))))
		my ($key,$seq,$struct,);
		$key = <IN>; $seq = <IN>; $struct = <IN>;
		chomp ($key,$seq,$struct);
		$key=~s/^>(\S+).*/$1/;
		if ($seq =~ m/$mat{$key}/) {
			my $pos = index($seq,$mat{$key});
			my $mat_start1=$pos+1;
			my $mat_end1=$mat_start1+length($mat{$key})-1;
			my $mat_start2=undef;
			my $mat_end2=undef;
			my $pair=0;
			my @unpair1=();
			my @unpair2=();
			my $before=undef;
			my $align=undef;
			my $mfe=undef;
			if ( $struct=~/^(\S+)\s+\(([^\)]+)/ ) {
				$align=$1;$mfe=$2;
				my @align=split //,$align;
				my %align=();
				my %head=();
				my $head=0;
				for (my $i=0;$i<@align;$i++) {
					if ($align[$i] eq ".") {
						$align{$i+1}=0;
					}elsif ($align[$i] eq "(") {
						$head++;
						$head{$head}=$i+1;
					}elsif ($align[$i] eq ")") {
						$align{$head{$head}}=$i+1;
						$align{$i+1}=$head{$head};
						delete $head{$head};
						$head--;
					}
				}
				my @key=sort {$a<=>$b} keys %align;
				for (my $i=0;$i<@key;$i++) {
					if ($key[$i]>=$mat_start1 and $key[$i]<=$mat_end1) {
						if ($align{$key[$i]}==0) {
							if (!defined $before) {
								$before=1;
							}else {
								if ($before==0) {
									$before=1;
								}else {
									$before++;
								}
							}
						}elsif ($align{$key[$i]}!=0) {
							if (!defined $mat_end2) {
								if (!defined $before) {
									$mat_end2=$align{$key[$i]};
								}else {
									$mat_end2=$align{$key[$i]}+$before;
								}
								$mat_start2=$align{$key[$i]};
							}else {
								$mat_start2=$align{$key[$i]};
							}
							if (!defined $before) {
								$before=0;
							}else {
								if ($before!=0) {
									push (@unpair1,$before);
									$before=0;
								}
							}
							$pair++;
						}
					}
				}
				if (defined $before and $before!=0) {
					push (@unpair1,$before);
					$mat_start2-=$before if (defined $mat_start2);
				}
				$before=undef;
				for (my $i=0;$i<@key;$i++) {
					if (defined $mat_start2 and defined $mat_end2 and $key[$i]>=$mat_start2 and $key[$i]<=$mat_end2) {
						if ($align{$key[$i]}==0) {
							if (!defined $before) {
								$before=1;
							}else {
								if ($before==0) {
									$before=1;
								}else {
									$before++;
								}
							}
						}elsif ($align{$key[$i]}!=0) {
							if (!defined $before) {
								$before=0;
							}else {
								if ($before!=0) {
									push (@unpair2,$before);
									$before=0;
								}
							}
						}
					}
				}
				if (defined $before and $before!=0) {
					push (@unpair2,$before);
				}
				@unpair1=sort {$b<=>$a} @unpair1;
				@unpair2=sort {$b<=>$a} @unpair2;
				my $unpair1=0;
				my $unpair2=0;
				foreach my $num (@unpair1) {
					$unpair1+=$num;
				}
				foreach my $num (@unpair2) {
					$unpair2+=$num;
				}
				if ($mat_start2>$mat_end1 or $mat_start1>$mat_end2) {
					if (defined $unpair1[0] and defined $unpair2[0] and $unpair1[0]<=$opt{maxbasebulge} and $unpair2[0]<=$opt{maxbasebulge} and $pair>=$opt{minbasepair} and $unpair1<=$opt{maxunpairbase} and $unpair2<=$opt{maxunpairbase} and $mat_end2>$mat_start2) {
						print FILTER "$key\t$mfe\t$pair\t$unpair1\t$unpair2\t$mat_start1\t$mat_end1\t$mat_start2\t$mat_end2\n";
						$seq_num++;
					}elsif (defined $unpair1[0] and $unpair1[0]<=$opt{maxbasebulge} and $pair>=$opt{minbasepair} and $unpair1<=$opt{maxunpairbase} and $unpair2<=$opt{maxunpairbase} and $mat_end2>$mat_start2) {
						print FILTER "$key\t$mfe\t$pair\t$unpair1\t$unpair2\t$mat_start1\t$mat_end1\t$mat_start2\t$mat_end2\n";
						$seq_num++;
					}elsif (defined $unpair2[0] and $unpair2[0]<=$opt{maxbasebulge} and $pair>=$opt{minbasepair} and $unpair1<=$opt{maxunpairbase} and $unpair2<=$opt{maxunpairbase} and $mat_end2>$mat_start2) {
						print FILTER "$key\t$mfe\t$pair\t$unpair1\t$unpair2\t$mat_start1\t$mat_end1\t$mat_start2\t$mat_end2\n";
						$seq_num++;
					}elsif ($pair>=$opt{minbasepair} and $unpair1<=$opt{maxunpairbase} and $unpair2<=$opt{maxunpairbase} and $mat_end2>$mat_start2) {
						print FILTER "$key\t$mfe\t$pair\t$unpair1\t$unpair2\t$mat_start1\t$mat_end1\t$mat_start2\t$mat_end2\n";
						$seq_num++;
					}
				}
				print ALL "$key\t$mfe\t$pair\t$unpair1\t$unpair2\t$mat_start1\t$mat_end1\t$mat_start2\t$mat_end2\n";
			}
		}
	}
	close IN;
	close FILTER;
	close ALL;

	open (STAT,">>$opt{prefix}.statistics")||die("Cannot write to $opt{prefix}.statistics.\n");
	print STAT "Extracted_seq_after_filter\t$seq_num\n";
	close STAT;

	$time=localtime;
	print LOG "The RNAfold result of precursor has been filtered: $time.\n\n";

	#Extract hairpin sequence and create hairpin gff file
	print LOG "Extract hairpin sequence and create hairpin gff file for filtered precursor.\n";
	my (%pre,%pos,%id,);
	open (PRE,"<$opt{prefix}\_pre.fa")||die("Cannot read $opt{prefix}\_pre.fa.\n");
	while (<PRE>) {
		chomp;
		if (/^>(\S+)/) {
			$key=$1;
		}else{
			$pre{$key}.=$_;
		}
	}
	close PRE;

	open (IN,"<$opt{prefix}\_pre.gff")||die("Cannot read $opt{prefix}\_pre.gff.\n");
	while (<IN>) {
		chomp;
		#S091595 .       miRNA   194     211     .       -       .       ACC="MI1"; ID="ath-miR5021";
		my @list=split /\t/,$_;
		if ($list[8]=~/ACC\=\"([^\"]+)\"/) {
			$pos{$1}=$_;
		}
	}
	close IN;

	open (FASTA,">$opt{prefix}\_hairpin.fa")||die("Cannot write to $opt{prefix}\_hairpin.fa file.\n");
	open (GFF,">$opt{prefix}\_hairpin.gff")||die("Cannot write to $opt{prefix}\_hairpin.gff.\n");
	open (FILTER,"<$opt{prefix}\_pre.filter")||die("fail to open $opt{prefix}\_pre.filter.\n");
	while (<FILTER>) {
		chomp;
		#miRNA  mfe     pair    unpair  unpair* mat_start       mat_end mat_start*      mat_end*
		#MI304252        -84.20  19      4       6       101     123     64      88
		next if (/^\#/);
		my @list=split /\t/,$_;
		$id{$list[0]}=0;
		my @pos=split /\t/,$pos{$list[0]};
		$pos[8]=~s/ACC\=*//;
		$pos[8]=~s/ID\=*//;
		$pos[8]=~s/EXPRESS\=*//;
		$pos[8]=~s/\"//g;
		$pos[8]=~s/\s+//g;
		my @feature=split /\;/,$pos[8];
		if ($list[5]<$list[7]) {
			my $seq=substr ($pre{$list[0]},$list[5]-1,$list[8]-$list[5]+1);
			print FASTA ">$list[0]\n$seq\n";
			if ($pos[6] eq "+") {
				my $pre_start=$pos[3]+$list[5]-1;
				my $pre_end=$pos[3]+$list[8]-1;
				#S091595 .       miRNA   194     211     .       -       .       ACC="MI1"; ID="ath-miR5021";
				print GFF "$pos[0]\tMITP\tmiRNA\t$pre_start\t$pre_end\t\.\t$pos[6]\t\.\tACC=\"$feature[0]\"; ID=\"$feature[1]\"; EXPRESS=\"$feature[2]\";\n";
			}else {
				my $pre_end=$pos[4]-$list[5]+1;
				my $pre_start=$pos[4]-$list[8]+1;
				print GFF "$pos[0]\tMITP\tmiRNA\t$pre_start\t$pre_end\t\.\t$pos[6]\t\.\tACC=\"$feature[0]\"; ID=\"$feature[1]\"; EXPRESS=\"$feature[2]\";\n";
			}
		}else {
			my $seq=substr ($pre{$list[0]},$list[7]-1,$list[6]-$list[7]+1);
			print FASTA ">$list[0]\n$seq\n";
			if ($pos[6] eq "+") {
				my $pre_start=$pos[3]+$list[7]-1;
				my $pre_end=$pos[3]+$list[6]-1;
				print GFF "$pos[0]\tMITP\tmiRNA\t$pre_start\t$pre_end\t\.\t$pos[6]\t\.\tACC=\"$feature[0]\"; ID=\"$feature[1]\"; EXPRESS=\"$feature[2]\";\n";
			}else {
				my $pre_end=$pos[4]-$list[7]+1;
				my $pre_start=$pos[4]-$list[6]+1;
				print GFF "$pos[0]\tMITP\tmiRNA\t$pre_start\t$pre_end\t\.\t$pos[6]\t\.\tACC=\"$feature[0]\"; ID=\"$feature[1]\"; EXPRESS=\"$feature[2]\";\n";
			}
		}
	}
	close FILTER;
	close FASTA;
	close GFF;

	$time=localtime;
	print LOG "The hairpin sequences and gff file have been extracted: $time.\n\n";

	#Filter redundant mature record from mature and hairpin fasta and gff file
	print LOG "Filter redundant mature record from mature and hairpin fasta and gff files.\n";

	%pre=();%pos=();
	my (%repeat,);
	open (IN,"<$opt{prefix}\_mature.gff")||die("Cannot read $opt{prefix}\_mature.gff.\n");
	while (<IN>) {
		chomp;
		#S091595 .       miRNA   194     211     .       -       .       ACC="MI1"; ID="ath-miR5021";
		my @list=split /\t/,$_;
		$list[8]=~s/ACC\=//;
		$list[8]=~s/ID\=//;
		$list[8]=~s/EXPRESS\=//;
		$list[8]=~s/\"//g;
		$list[8]=~s/\;//g;
		my @feature=split /\s+/,$list[8];
		if (defined $id{$feature[0]}) {
			$repeat{$list[0]}{$list[6]}{"$list[3]\t$list[4]"}[0]=$list[3];
			$repeat{$list[0]}{$list[6]}{"$list[3]\t$list[4]"}[1]=$feature[0];
			$repeat{$list[0]}{$list[6]}{"$list[3]\t$list[4]"}[2]=$feature[2];
			$repeat{$list[0]}{$list[6]}{"$list[3]\t$list[4]"}[3]=$_;
		}
	}
	close IN;

	open (GFF,">$opt{prefix}\_mature_unique.gff")||die("Cannot write to $opt{prefix}\_mature_unique.gff file.\n");
	my $seq_unique=0;
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
					if (($end-$start+1)/($end1-$start1+1)>=$opt{matoverrate}/100 or ($end-$start+1)/($end2-$start2+1)>=$opt{matoverrate}/100) {
						if ($repeat{$chr}{$strand}{$record[$i]}[2]>=$repeat{$chr}{$strand}{$record[$i+1]}[2]) {
							delete $repeat{$chr}{$strand}{$record[$i+1]};
							$record[$i+1]=$record[$i];
						}else {
							delete $repeat{$chr}{$strand}{$record[$i]};
						}
					}
				}
			}
			@record=sort {$repeat{$chr}{$strand}{$a}[0]<=>$repeat{$chr}{$strand}{$b}[0]} keys %{$repeat{$chr}{$strand}};
			for (my $i=0;$i<@record;$i++) {
				print GFF "$repeat{$chr}{$strand}{$record[$i]}[3]\n";
				$id{$repeat{$chr}{$strand}{$record[$i]}[1]}=1;
				$seq_unique++;
			}
		}
	}
	close GFF;
	%repeat=();

	open (STAT,">>$opt{prefix}.statistics")||die("Cannot write to $opt{prefix}.statistics.\n");
	print STAT "Extracted_seq_unique\t$seq_unique\n\n";
	close STAT;

	$time=localtime;
	print LOG "The redundant mature record in mature gff file has been filtered: $time.\n";

	open (FASTA,">$opt{prefix}\_mature_unique.fa")||die("Cannot write to $opt{prefix}\_mature_unique.fa file.\n");
	open (IN,"<$opt{prefix}\_mature.fa")||die("Cannot read $opt{prefix}\_mature.fa.\n");
	my $flag=0;
	while (<IN>) {
		chomp;
		if (/^>(\S+)/) {
			if (defined $id{$1} and $id{$1}>0) {
				print FASTA ">$1\n";
				$flag=1;
			}else {
				$flag=0;
			}
			next;
		}
		print FASTA "$_\n" if ($flag==1);
	}
	close IN;
	close FASTA;

	$time=localtime;
	print LOG "The redundant mature sequence in mature fasta file has been filtered: $time.\n";

	open (GFF,">$opt{prefix}\_hairpin_unique.gff")||die("Cannot write to $opt{prefix}\_hairpin_unique.gff file.\n");
	open (IN,"<$opt{prefix}\_hairpin.gff")||die("Cannot read $opt{prefix}\_hairpin.gff.\n");
	while (<IN>) {
		chomp;
		#S091595 .       miRNA   194     211     .       -       .       ACC="MI1"; ID="ath-miR5021";
		my @list=split /\t/,$_;
		$list[8]=~s/ACC\=//;
		$list[8]=~s/ID\=//;
		$list[8]=~s/EXPRESS\=//;
		$list[8]=~s/\"//g;
		$list[8]=~s/\;//g;
		my @feature=split /\s+/,$list[8];
		if (defined $id{$feature[0]} and $id{$feature[0]}>0) {
			print GFF "$_\n";
		}
	}
	close IN;
	close GFF;

	$time=localtime;
	print LOG "The redundant hairpin record in hairpin gff file has been filtered: $time.\n";

	open (FASTA,">$opt{prefix}\_hairpin_unique.fa")||die("Cannot write to $opt{prefix}\_hairpin_unique.fa file.\n");
	open (IN,"<$opt{prefix}\_hairpin.fa")||die("Cannot read $opt{prefix}\_hairpin.fa.\n");
	$flag=0;
	while (<IN>) {
		chomp;
		if (/^>(\S+)/) {
			if (defined $id{$1} and $id{$1}>0) {
				print FASTA ">$1\n";
				$flag=1;
			}else {
				$flag=0;
			}
			next;
		}
		print FASTA "$_\n" if ($flag==1);
	}
	close IN;
	close FASTA;

	$time=localtime;
	print LOG "The redundant hairpin sequence in hairpin fasta file has been filtered: $time.\n";
	%id=();

	#Do RNAfold for unique hairpin sequence
	print LOG "Do RNAfold for unique hairpin sequences.\n";
	system "less $opt{prefix}\_hairpin_unique.fa | RNAfold --noconv --noPS > $opt{prefix}\_hairpin_unique.RNAfold";

	$time=localtime;
	print LOG "The RNAfold for unique hairpin sequences has been done: $time.\n\n";
}

if ($opt{step} eq "cluster" or $opt{step} eq "filter" or $opt{step} eq "extract_seq" or $opt{step} eq "candidate" or $opt{step} eq "candidate_filter") {

	#Obtain hairpin sequence and structure attribute values
	print LOG "Obtain hairpin sequence and structure attribute values.\n";
	open (OUT,">$opt{prefix}\_hairpin_unique.attribute")||die("Cannot write to $opt{prefix}\_hairpin_unique.attribute file.\n");
	print OUT "#miRNA\tHairpin_length\tPair_base\tPair_percent(%)\tA_content(%)\tU_content(%)\tG_content(%)\tC_content(%)\tAU_content(%)\tGC_content(%)\tMFE\tAMFE\tMFEI\n";
	open (IN,"<$opt{prefix}\_hairpin_unique.RNAfold")||die("Cannot read $opt{prefix}\_hairpin_unique.RNAfold.\n");
	while (!eof(IN)) {
		my ($key,$seq,$struct,);
		$key = <IN>; $seq = <IN>; $struct = <IN>;
		chomp ($key,$seq,$struct);
		$key=~s/^>(\S+).*/$1/;
		my $length=length ($seq);
		my $pair=0;
		my $mfe=undef;
		if ( $struct=~/^(\S+)\s+\(([^\)]+)/ ) {
			my $align=$1;my $mfe=$2;
			my @align=split //,$align;
			for (my $i=0;$i<@align;$i++) {
				if ($align[$i] eq "(") {
					$pair++;
				}elsif ($align[$i] eq ")") {
					$pair++;
				}
			}
			my ($a_content,$u_content,$g_content,$c_content,$au_content,$gc_content)=&base_content($seq);
			my $pair_percent=sprintf "%.2f",($pair)/($length)*100;
			my $amfe=sprintf "%.2f",$mfe/$length*100;
			my $mfei=();
			if ($gc_content>0) {
				$mfei=sprintf "%.2f",$amfe/$gc_content;
			}else {
				$mfei=sprintf "%.2f",$amfe/4;
			}
			print OUT "$key\t$length\t$pair\t$pair_percent\t$a_content\t$u_content\t$g_content\t$c_content\t$au_content\t$gc_content\t$mfe\t$amfe\t$mfei\n";
		}
	}
	close IN;
	close OUT;

	$time=localtime;
	print LOG "The hairpin sequence and structure attribute values have been obtained: $time.\n\n";

	%id=();
	#Filter candidate miRNA according to the miRNA sequence and structure attribute values
	print LOG "Filter candidate miRNA according to the miRNA sequence and structure attribute values.\n";
	open (OUT,">$opt{prefix}\_hairpin_unique.filter")||die("Cannot write to $opt{prefix}\_hairpin_unique.filter file.\n");
	print OUT "#miRNA\tHairpin_length\tPair_base\tPair_percent(%)\tA_content(%)\tU_content(%)\tG_content(%)\tC_content(%)\tAU_content(%)\tGC_content(%)\tMFE\tAMFE\tMFEI\n";
	open (IN,"<$opt{prefix}\_hairpin_unique.attribute")||die("Cannot read $opt{prefix}\_hairpin_unique.attribute.\n");
	my $attribute_filter=0;
	while (<IN>) {
		chomp;
		#MI152895        68      46      67.65   10.29   11.76   32.35   45.59   22.06   77.94   -45.00  -66.18  -0.85
		next if (/^\#/);
		my @list=split /\t/,$_;		
		next if ($list[12]>$opt{mfei});
		next if ($list[10]>$opt{mfe});
		next if ($list[1]<$opt{hairlen});
		next if (defined $opt{pairbase} and $list[2]<$opt{pairbase});
		next if (defined $opt{pairpercent} and $list[3]<$opt{pairpercent});
		next if (defined $opt{amfe} and $list[11]>$opt{amfe});
		next if (defined $opt{minapercent} and $list[4]<$opt{minapercent});
		next if (defined $opt{maxapercent} and $list[4]>$opt{maxapercent});
		next if (defined $opt{minupercent} and $list[4]<$opt{minupercent});
		next if (defined $opt{maxupercent} and $list[4]>$opt{maxupercent});
		next if (defined $opt{mingpercent} and $list[4]<$opt{mingpercent});
		next if (defined $opt{maxgpercent} and $list[4]>$opt{maxgpercent});
		next if (defined $opt{mincpercent} and $list[4]<$opt{mincpercent});
		next if (defined $opt{maxcpercent} and $list[4]>$opt{maxcpercent});
		next if (defined $opt{minaupercent} and $list[4]<$opt{minaupercent});
		next if (defined $opt{maxaupercent} and $list[4]>$opt{maxaupercent});
		next if (defined $opt{mingcpercent} and $list[4]<$opt{mingcpercent});
		next if (defined $opt{maxgcpercent} and $list[4]>$opt{maxgcpercent});
		print OUT "$_\n";
		$id{$list[0]}=$_;
		$attribute_filter++;
	}
	close IN;
	close OUT;

	open (STAT,">>$opt{prefix}.statistics")||die("Cannot write to $opt{prefix}.statistics.\n");
	print STAT "Candidate_attribute_filter\t$attribute_filter\n";
	close STAT;

	if (defined $opt{selfalign}) {

		system "perl $programdir/fish.pl -tb list -tf fa -b $opt{prefix}\_hairpin_unique.filter -f $opt{prefix}\_hairpin_unique.fa -o $opt{prefix}\_hairpin_unique.filter.fa";
		system "perl $programdir/fish.pl -tb list -tf fa -b $opt{prefix}\_hairpin_unique.filter -f $opt{prefix}\_mature_unique.fa -o $opt{prefix}\_mature_unique.filter.fa";
		
		if ($opt{alignsoft} eq "blast") {
			system "formatdb -i $opt{prefix}\_hairpin_unique.filter.fa -p F";
			system "formatdb -i $opt{prefix}\_mature_unique.filter.fa -p F";
			system "blastall -p blastn -d $opt{prefix}\_hairpin_unique.filter.fa -i $opt{prefix}\_hairpin_unique.filter.fa -v 1000 -b 1000 -W 7 -o $opt{prefix}\_hairpin_unique.filter.blast";
			system "blastall -p blastn -d $opt{prefix}\_mature_unique.filter.fa -i $opt{prefix}\_mature_unique.filter.fa -v 1000 -b 1000 -W 7 -o $opt{prefix}\_mature_unique.filter.blast";
			system "perl $programdir/EblastN.pl -i $opt{prefix}\_hairpin_unique.filter.blast -o $opt{prefix}\_hairpin_unique.filter.eblastn";
			system "perl $programdir/EblastN.pl -i $opt{prefix}\_mature_unique.filter.blast -o $opt{prefix}\_mature_unique.filter.eblastn";

			my (%redundant1,%redundant2,%mirna,);
			open (IN,"<$opt{prefix}\_mature_unique.filter.eblastn")||die("fail to open $opt{prefix}\_mature_unique.filter.eblastn.\n");
			while (<IN>) {
				chomp;
				next if (/^Query/);
				#MI34    18      1       18      1       18      19      36.2    5e-06   18/18   100     MI21405 ptc-miR169f
				my @list=split /\t/,$_;
				next if ($list[0] eq $list[11]);
				next if ($list[4] > $list[5]);
				my $ident=$list[10];
				my $arate1=($list[3]-$list[2]+1)/$list[1];
				my $arate2=($list[5]-$list[4]+1)/$list[6];	
				if ($ident>=$opt{selfidentity} and ($arate1>=$opt{selfoverrate} or $arate2>=$opt{selfoverrate})) {
					my $flag=0;
					foreach my $key (keys %redundant1) {
						my @id=split /\,/,$redundant1{$key};
						for (my $i=0;$i<@id;$i++) {
							if ($list[0] eq $id[$i] or $list[11] eq $id[$i]) {
								$flag=1;
								$redundant1{$key}.=$list[0]."," if ($redundant1{$key}!~/$list[0],/);
								$redundant1{$key}.=$list[11]."," if ($redundant1{$key}!~/$list[11],/);
							}
							last if ($flag==1);
						}
						last if ($flag==1);
					}
					$redundant1{$list[0]}.=$list[0]."," if ($flag==0);
					$redundant1{$list[0]}.=$list[11]."," if ($flag==0);
				}
			}
			close IN;

			open (IN,"<$opt{prefix}\_hairpin_unique.filter.eblastn")||die("fail to open $opt{prefix}\_hairpin_unique.filter.eblastn.\n");
			while (<IN>) {
				chomp;
				next if (/^Query/);
				#MI34    18      1       18      1       18      19      36.2    5e-06   18/18   100     MI21405 ptc-miR169f
				my @list=split /\t/,$_;
				next if ($list[0] eq $list[11]);
				next if ($list[4] > $list[5]);
				my $ident=$list[10];
				my $arate1=($list[3]-$list[2]+1)/$list[1];
				my $arate2=($list[5]-$list[4]+1)/$list[6];
				if ($ident>=$opt{selfidentity} and ($arate1>=$opt{selfoverrate} or $arate2>=$opt{selfoverrate})) {
					if (!defined $mirna{$list[0]} or !defined $mirna{$list[11]}) {
						foreach my $key1 (keys %redundant1) {
							if ($redundant1{$key1}=~/$list[0],/ and $redundant1{$key1}=~/$list[11],/) {
								my $flag=0;
								foreach my $key (keys %redundant2) {
									my @id=split /\,/,$redundant2{$key};
									for (my $i=0;$i<@id;$i++) {
										if ($list[0] eq $id[$i] or $list[11] eq $id[$i]) {
											$flag=1;
											$redundant2{$key}.=$list[0]."," if ($redundant2{$key}!~/$list[0],/);
											$redundant2{$key}.=$list[11]."," if ($redundant2{$key}!~/$list[11],/);
										}
										last if ($flag==1);
									}
									last if ($flag==1);
								}
								$redundant2{$list[0]}.=$list[0]."," if ($flag==0);
								$redundant2{$list[0]}.=$list[11]."," if ($flag==0);
								$mirna{$list[0]}++;
								$mirna{$list[11]}++;
							}
						}
					}
				}
			}
			close IN;

			open (OUT,">$opt{prefix}\_hairpin_unique.filter2")||die("fail to open $opt{prefix}\_hairpin_unique.filter2.\n");
			foreach my $key (keys %redundant2) {
				my @id=split /\,/,$redundant2{$key};
				for (my $i=0;$i<@id;$i++) {
					delete $id{$id[$i]} if ($id[$i] ne $key);
				}
			}
			my $selfalign_filter=0;
			foreach my $id (%id) {
				print OUT "$id{$id}\n";
				$selfalign_filter++;
			}
			close OUT;

			open (STAT,">>$opt{prefix}.statistics")||die("Cannot write to $opt{prefix}.statistics.\n");
			print STAT "Candidate_selfalign_filter\t$selfalign_filter\n\n";
			close STAT;

			%redundant1=();
			%redundant2=();
			%mirna=();
		}else {
			system "blat $opt{prefix}\_hairpin_unique.filter.fa $opt{prefix}\_hairpin_unique.filter.fa -tileSize=8 -oneOff=1 -minMatch=1 -minIdentity=80 -noHead $opt{prefix}\_hairpin_unique.filter.blat";
			system "blat $opt{prefix}\_mature_unique.filter.fa $opt{prefix}\_mature_unique.filter.fa -tileSize=8 -oneOff=1 -minMatch=1 -minIdentity=80 -noHead $opt{prefix}\_mature_unique.filter.blat";

			my (%redundant1,%redundant2,%mirna,);
			open (IN,"<$opt{prefix}\_mature_unique.filter.blat")||die("fail to open $opt{prefix}\_mature_unique.filter.blat.\n");
			while (<IN>) {
				chomp;
				next if (/^[\D]+/);
				my @list=split /\t/,$_;
				next if ($list[9] eq $list[13]);
				next if ($list[8] eq "-");
				my $ident=$list[0]/($list[0]+$list[1])*100;
				my $arate1=($list[12]-$list[11])/$list[10];
				my $arate2=($list[16]-$list[15])/$list[14];	
				if ($ident>=$opt{selfidentity} and ($arate1>=$opt{selfoverrate} or $arate2>=$opt{selfoverrate})) {
					my $flag=0;
					foreach my $key (keys %redundant1) {
						my @id=split /\,/,$redundant1{$key};
						for (my $i=0;$i<@id;$i++) {
							if ($list[9] eq $id[$i] or $list[13] eq $id[$i]) {
								$flag=1;
								$redundant1{$key}.=$list[9]."," if ($redundant1{$key}!~/$list[9],/);
								$redundant1{$key}.=$list[13]."," if ($redundant1{$key}!~/$list[13],/);
							}
							last if ($flag==1);
						}
						last if ($flag==1);
					}
					$redundant1{$list[9]}.=$list[9]."," if ($flag==0);
					$redundant1{$list[9]}.=$list[13]."," if ($flag==0);
				}
			}
			close IN;

			open (IN,"<$opt{prefix}\_hairpin_unique.filter.blat")||die("fail to open $opt{prefix}\_hairpin_unique.filter.blat.\n");
			while (<IN>) {
				chomp;
				next if (/^[\D]+/);
				my @list=split /\t/,$_;
				next if ($list[9] eq $list[13]);
				next if ($list[8] eq "-");
				my $ident=$list[0]/($list[0]+$list[1])*100;
				my $arate1=($list[12]-$list[11])/$list[10];
				my $arate2=($list[16]-$list[15])/$list[14];	
				if ($ident>=$opt{selfidentity} and ($arate1>=$opt{selfoverrate} or $arate2>=$opt{selfoverrate})) {
					if (!defined $mirna{$list[0]} or !defined $mirna{$list[11]}) {
						foreach my $key1 (keys %redundant1) {
							if ($redundant1{$key1}=~/$list[0],/ and $redundant1{$key1}=~/$list[11],/) {
								my $flag=0;
								foreach my $key (keys %redundant2) {
									my @id=split /\,/,$redundant2{$key};
									for (my $i=0;$i<@id;$i++) {
										if ($list[9] eq $id[$i] or $list[13] eq $id[$i]) {
											$flag=1;
											$redundant2{$key}.=$list[9]."," if ($redundant2{$key}!~/$list[9],/);
											$redundant2{$key}.=$list[13]."," if ($redundant2{$key}!~/$list[13],/);
										}
										last if ($flag==1);
									}
									last if ($flag==1);
								}
								$redundant2{$list[9]}.=$list[9]."," if ($flag==0);
								$redundant2{$list[9]}.=$list[13]."," if ($flag==0);
								$mirna{$list[0]}++;
								$mirna{$list[11]}++;
							}
						}
					}
				}
			}
			close IN;

			open (OUT,">$opt{prefix}\_hairpin_unique.filter2")||die("fail to open $opt{prefix}\_hairpin_unique.filter2.\n");
			foreach my $key (keys %redundant2) {
				my @id=split /\,/,$redundant2{$key};
				for (my $i=0;$i<@id;$i++) {
					delete $id{$id[$i]} if ($id[$i] ne $key);
				}
			}
			my $selfalign_filter=0;
			foreach my $id (%id) {
				print OUT "$id{$id}\n";
				$selfalign_filter++;
			}
			close OUT;

			open (STAT,">>$opt{prefix}.statistics")||die("Cannot write to $opt{prefix}.statistics.\n");
			print STAT "Candidate_selfalign_filter\t$selfalign_filter\n\n";
			close STAT;

			%redundant1=();
			%redundant2=();
			%mirna=();
		}
	}

	$time=localtime;
	print LOG "The candidate hairpin sequence was filtered according to the hairpin sequence and structure attribute values: $time.\n\n";

	#Obtain the mature and hairpin sequence and gff file for final candidate miRNA
	print LOG "Obtain the mature and hairpin sequence and gff file for final candidate miRNA.\n";

	if (-e "$opt{prefix}\_hairpin_unique.filter2") {
		system "perl $programdir/fish.pl -tb list -tf fa -b $opt{prefix}\_hairpin_unique.filter2 -f $opt{prefix}\_hairpin_unique.fa -o $opt{prefix}\_candidate_hairpin.fa";
		system "perl $programdir/fish.pl -tb list -tf fa -b $opt{prefix}\_hairpin_unique.filter2 -f $opt{prefix}\_mature_unique.fa -o $opt{prefix}\_candidate_mature.fa";
		system "perl $programdir/fish.pl -tb list -tf gff -b $opt{prefix}\_hairpin_unique.filter2 -f $opt{prefix}\_hairpin_unique.gff -o $opt{prefix}\_candidate_hairpin.gff";
		system "perl $programdir/fish.pl -tb list -tf gff -b $opt{prefix}\_hairpin_unique.filter2 -f $opt{prefix}\_mature_unique.gff -o $opt{prefix}\_candidate_mature.gff";
	}else {
		system "perl $programdir/fish.pl -tb list -tf fa -b $opt{prefix}\_hairpin_unique.filter -f $opt{prefix}\_hairpin_unique.fa -o $opt{prefix}\_candidate_hairpin.fa";
		system "perl $programdir/fish.pl -tb list -tf fa -b $opt{prefix}\_hairpin_unique.filter -f $opt{prefix}\_mature_unique.fa -o $opt{prefix}\_candidate_mature.fa";
		system "perl $programdir/fish.pl -tb list -tf gff -b $opt{prefix}\_hairpin_unique.filter -f $opt{prefix}\_hairpin_unique.gff -o $opt{prefix}\_candidate_hairpin.gff";
		system "perl $programdir/fish.pl -tb list -tf gff -b $opt{prefix}\_hairpin_unique.filter -f $opt{prefix}\_mature_unique.gff -o $opt{prefix}\_candidate_mature.gff";
	}
	$time=localtime;
	print LOG "The mature and hairpin sequence and gff file for final candidate miRNA have been obtained: $time.\n\n";

	if ($opt{figure} eq "yes") {

		#Obtain the second structure for hairpin sequence of candidate miRNA.
		print LOG "Obtain the second structure for hairpin sequence of candidate miRNA.\n";

		system "mkdir $opt{prefix}\_figure";

		chdir "$opt{output_dir}/$opt{prefix}\_figure";
		
		open (IN,"<$opt{output_dir}/$opt{prefix}\_hairpin_unique.RNAfold")||die("Cannot read $opt{output_dir}/$opt{prefix}\_hairpin_unique.RNAfold.\n");
		while (!eof(IN)) {
			my ($key,$seq,$struct,);
			$key = <IN>; $seq = <IN>; $struct = <IN>;
			chomp ($key,$seq,$struct);
			$key=~s/^>(\S+).*/$1/;
			if (defined $id{$key}) {
				open (OUT,">$key.seq")||die("Cannot write to $key.seq file.\n");
				print OUT ">$key\n$seq\n$struct\n";
				close OUT;

				system "less $key.seq | $programdir/b2ct > $key.ct";
				system "$programdir/sir_graph -p $key.ct >$key.ss";
				system "$programdir/epstopdf $key.ps";
				system "rm $key.seq $key.ct $key.ss $key.ps";
			}
		}
		close IN;

		chdir "$opt{output_dir}";

		$time=localtime;
		print LOG "The second structure for hairpin sequence of candidate miRNA has obtained: $time.\n\n";
	}
}

if ($opt{step} eq "cluster" or $opt{step} eq "filter" or $opt{step} eq "extract_seq" or $opt{step} eq "candidate" or $opt{step} eq "candidate_filter" or $opt{step} eq "expression") {

	if (defined $opt{expression}) {
		
		#Read statistics file
		print LOG "Read statistics file.\n";
		my $total_read=0;
		open (IN,"<$opt{prefix}.statistics")||die("Cannot read $opt{prefix}.statistics.\n");
		while (<IN>) {
			chomp;
			if (/Filtered_read\s+(\d+)/) {
				$total_read=$1;
				last;
			}			
		}
		close IN;

		$time=localtime;
		print LOG "Statistics file has been readed: $time.\n\n";

		print LOG "Read precursor gff file.\n";

		my %pos=();
		my %hash=();
		open (IN,"<$opt{prefix}\_pre.gff")||die("Cannot read $opt{prefix}\_pre.gff.\n");
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
		print LOG "The precursor gff file has been readed: $time.\n\n";

		print LOG "Read $opt{prefix}\_read.cluster_filter file.\n";

		open (IN,"<$opt{prefix}\_read.cluster_filter")||die("Cannot read $opt{prefix}\_read.cluster_filter.\n");
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
		print LOG "The $opt{prefix}\_read.cluster_filter file has been readed: $time.\n\n";

		print LOG "Read $opt{prefix}\_pre.filter file and calulate the expression for candidate miRNA.\n";

		open (OUT,">$opt{prefix}\_candidate.expression")||die("Cannot write to $opt{prefix}\_candidate.expression file.\n");
		print OUT "#miRNA\tchr\tstrand\thairpin_start\thairpin_end\thairpin_read\thairpin_express(read_number/total_read*1M)\tmature_start\tmature_end\tmature_read\tmature_express\tstar_start\tstar_end\tstar_read\tstar_express\n";
		my %pre=();
		open (IN,"<$opt{prefix}\_pre.filter")||die("fail to open $opt{prefix}\_pre.filter.\n");
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
					ECOM:
					if ($hend<$pstart) {
						$h++;
						if (defined $hash[$h]) {
							($hstart,$hend)=split /\t/,$hash[$h];
							goto ECOM;
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
							goto ECOM;
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
		print LOG "miRNA expresssion values have been obtained: $time.\n\n";
	}
}

if ($opt{step} eq "cluster" or $opt{step} eq "filter" or $opt{step} eq "extract_seq" or $opt{step} eq "candidate" or $opt{step} eq "candidate_filter" or $opt{step} eq "expression" or $opt{step} eq "class") {

	if (defined $opt{conserveseq}) {
		#Identify conserve miRNA from candidate miRNA
		print LOG "Identify conserve miRNA from candidate miRNA.\n";

		my ($key,);

		$key=undef;
		$string=undef;
		open (OUT,">$opt{conserveseq}.new")||die("fail to open $opt{conserveseq}.new.\n");
		open (IN,"<$opt{conserveseq}")||die("fail to open $opt{conserveseq}.\n");
		while (<IN>) {
			chomp;
			if (/^\#/) {
				next;
			}elsif (/^>(.*)/) {
				if (defined $key) {
					print OUT ">$key\n";
					$string=~tr/Uu/Tt/;
					print OUT "$string\n";
					$string=undef;
				}
				$key=$1;
				next;
			}
			$string.=$_;
		}
		print OUT ">$key\n";
		$string=~tr/Uu/Tt/;
		print OUT "$string\n";
		$key=undef;
		$string=undef;
		close IN;
		close OUT;

		%id=();
		my %sample=();
		open (IN,"<$opt{prefix}\_candidate_mature.fa")||die("fail to open $opt{prefix}\_candidate_mature.fa.\n");
		while (<IN>) {
			chomp;
			if (/^\#/) {
				next;
			}elsif (/^>(\S+)/) {
				$id{$1}=1;
			}
		}
		close IN;

		if ($opt{alignsoft} eq "blast") {

			system "formatdb -i $opt{conserveseq}.new -p F";
			system "blastall -p blastn -d $opt{conserveseq}.new -i $opt{prefix}\_candidate_mature.fa -v 1000 -b 1000 -W 7 -o $opt{prefix}\_candidate_mature.blast";
			system "perl $programdir/EblastN.pl -i $opt{prefix}\_candidate_mature.blast -o $opt{prefix}\_candidate_mature.eblastn";

			my %conserve=();
			open (IN,"<$opt{prefix}\_candidate_mature.eblastn")||die("Cannot read $opt{prefix}\_candidate_mature.eblastn.\n");
			while (<IN>) {
				next if (/^Query/);
				my @list=split(/\t/,$_);
				my ($overlap,$total)=split /\//,$list[9];
				my $mis=$list[6]-(abs($list[5]-$list[4])+1)+($total-$overlap);
				if ($list[3]-$list[2]+1>=$opt{minimum} and $list[3]-$list[2]+1<=$opt{maximum} and $list[10]>=$opt{identity} and $mis<=$opt{mismatch}) {
					$conserve{$list[0]}.=$list[11]."\t";
				}
			}
			close IN;

			open (PREDICT,">$opt{prefix}\_candidate.class")||die("Cannot write to $opt{prefix}\_candidate.class.\n");
			my @id=keys %conserve;
			my $conservenum=$#id+1;

			open (STAT,">>$opt{prefix}.statistics")||die("Cannot write to $opt{prefix}.statistics.\n");
			print STAT "Conserve_miRNA\t$conservenum\n";
			close STAT;

			print PREDICT "#Conserve miRNA $conservenum\n";
			print PREDICT "#ID\t#conserveID\n";
			foreach my $id (@id) {
				print PREDICT "$id\t$conserve{$id}\n";
				my @name=split /\t/,$conserve{$id};
				if ($name[0]=~/[A-Za-z]+\-([A-Za-z]+\-\d+)/) {
					$sample{$1}++;
				}
			}
			my $novelnum=(keys %id)-$conservenum;

			open (STAT,">>$opt{prefix}.statistics")||die("Cannot write to $opt{prefix}.statistics.\n");
			print STAT "Novel_miRNA\t$novelnum\n";
			close STAT;

			print PREDICT "\n#Novel miRNA $novelnum\n";
			print PREDICT "#ID\n";
			foreach my $id (keys %id) {
				if (!defined $conserve{$id}) {
					print PREDICT "$id\n";
				}
			}
		}else {

			system "blat $opt{conserveseq}.new $opt{prefix}\_candidate_mature.fa -tileSize=8 -oneOff=1 -minMatch=1 -minIdentity=80 -noHead $opt{prefix}\_candidate_mature.blat";
			
			my %conserve=();
			open (IN,"<$opt{prefix}\_candidate_mature.blat")||die("Cannot read $opt{prefix}\_candidate_mature.blat.\n");
			while (<IN>) {
				chomp;
				my @list=split /\t/,$_;
				my $identity=$list[0]/($list[0]+$list[1])*100;
				my $mis=$list[1]+$list[5]+$list[7];
				if ($list[12]-$list[11]>=$opt{minimum} and $list[12]-$list[11]<=$opt{maximum} and $identity>=$opt{identity} and $mis<=$opt{mismatch}) {
					$conserve{$list[9]}.=$list[13]."\t";
				}
			}
			close IN;

			open (PREDICT,">$opt{prefix}\_candidate.class")||die("Cannot write to $opt{prefix}\_candidate.class.\n");
			my @id=keys %conserve;
			my $conservenum=$#id+1;

			open (STAT,">>$opt{prefix}.statistics")||die("Cannot write to $opt{prefix}.statistics.\n");
			print STAT "Conserve_miRNA\t$conservenum\n";
			close STAT;

			print PREDICT "#Conserve miRNA $conservenum\n";
			print PREDICT "#ID\t#conserveID\n";
			foreach my $id (@id) {
				print PREDICT "$id\t$conserve{$id}\n";
				my @name=split /\t/,$conserve{$id};
				if ($name[0]=~/[A-Za-z]+\-([A-Za-z]+\-\d+)/) {
					$sample{$1}++;
				}
			}
			my $novelnum=(keys %id)-$conservenum;

			open (STAT,">>$opt{prefix}.statistics")||die("Cannot write to $opt{prefix}.statistics.\n");
			print STAT "Novel_miRNA\t$novelnum\n";
			close STAT;

			print PREDICT "\n#Novel miRNA $novelnum\n";
			print PREDICT "#ID\n";
			foreach my $id (keys %id) {
				if (!defined $conserve{$id}) {
					print PREDICT "$id\n";
				}
			}
		}

		$time=localtime;
		print LOG "The conserve miRNA in candidate miRNA has been identified: $time.\n";

		print LOG "Do miRNA compare for selected multiple species and this sample.\n";

		my %fam=();my %mirna=();

		open (IN,"<$opt{conserveseq}")||die("fail to open $opt{conserveseq}.\n");
		while (<IN>) {
			chomp;
			#>cel-let-7 MI0000001 Caenorhabditis elegans let-7 stem-loop
			#>cel-mir-1 MI0000003 Caenorhabditis elegans miR-1 stem-loop
			if (/^>([^\-]+)\-([^\-]+\-\d+)/) {
				if (defined $sample{$2}) {
					$mirna{$1}{$2}++;
					$fam{$2}++;
				}
			}
		}
		close IN;

		open (OUT,">$opt{prefix}\_multiple_species.compare")||die("Cannot write to $opt{prefix}\_multiple_species.compare file.\n");

		foreach my $fam (keys %fam) {
			$fam{$fam}+=$sample{$fam};
		}

		my @fam=sort {$fam{$b} <=> $fam{$a}} keys %fam;
		print OUT "Species\t";
		foreach my $fam (@fam) {
			print OUT "$fam\t";
		}
		print OUT "Total\n";

		#print the miRNA number of different families in sample
		my $total=0;
		print OUT "new\t";
		my $species_total=0;
		foreach my $fam (@fam) {
			$species_total+=$sample{$fam};
		}

		foreach my $fam (@fam) {
			printf OUT "$sample{$fam}(%.2f)\t",$sample{$fam}/$species_total*100;
		}
		print OUT "$species_total\n";
		$total+=$species_total;

		#print the miRNA number of different families in other species
		foreach my $species (keys %mirna) {
			print OUT "$species\t";
			$species_total=0;
			foreach my $fam (@fam) {
				if (defined $mirna{$species}{$fam}) {
					$species_total+=$mirna{$species}{$fam};
				}
			}
			foreach my $fam (@fam) {
				if (defined $mirna{$species}{$fam}) {
					printf OUT "$mirna{$species}{$fam}(%.2f)\t",$mirna{$species}{$fam}/$species_total*100;
				}else {
					printf OUT "0(0.00)\t";
				}
			}
			print OUT "$species_total\n";
			$total+=$species_total;
		}

		print OUT "Total\t";
		foreach my $fam (@fam) {
			printf OUT "$fam{$fam}(%.2f)\t",$fam{$fam}/$total*100;
		}
		print OUT "$total\n";
		close OUT;

		$time=localtime;
		print LOG "The miRNA comparison for selected multiple species has been finished: $time.\n\n";
	}
}

if ($opt{step} eq "cluster" or $opt{step} eq "filter" or $opt{step} eq "extract_seq" or $opt{step} eq "candidate" or $opt{step} eq "candidate_filter" or $opt{step} eq "expression" or $opt{step} eq "class" or $opt{step} eq "target") {
	my ($key,$string,$time,);

	if (defined $opt{target_seq}) {

		#Identify target for candidate miRNA
		print LOG "Identify target for candidate miRNA.\n";

		my %bp;
		$bp{"AT"} = 0;
		$bp{"TA"} = 0;
		$bp{"GC"} = 0;
		$bp{"CG"} = 0;
		$bp{"GT"} = 0.5;
		$bp{"TG"} = 0.5;
		$bp{"AC"} = 1;
		$bp{"CA"} = 1;
		$bp{"AG"} = 1;
		$bp{"GA"} = 1;
		$bp{"TC"} = 1;
		$bp{"CT"} = 1;
		$bp{"AA"} = 1;
		$bp{"TT"} = 1;
		$bp{"CC"} = 1;
		$bp{"GG"} = 1;

		$key=undef;
		my $info=undef;
		$string=undef;
		my %target=();
		my %tinfo=();
		open (OUT,">$opt{target_seq}.new")||die("fail to open $opt{target_seq}.new.\n");
		open (IN,"<$opt{target_seq}")||die("fail to open $opt{target_seq}.\n");
		while (<IN>) {
			chomp;
			if (/^\#/) {
				next;
			}elsif (/^>(\S+)\s*(.*)/) {
				if (defined $key) {
					print OUT ">$key $info\n";
					$string=~tr/Uu/Tt/;
					$string=uc $string;
					print OUT "$string\n";
					$target{$key}=$string;
					$tinfo{$key}=$info;
					$string=undef;
				}
				$key=$1;
				$info=$2;
				next;
			}
			$string.=$_;
		}
		print OUT ">$key $info\n";
		$string=~tr/Uu/Tt/;
		$string=uc $string;
		print OUT "$string\n";
		$target{$key}=$string;
		$tinfo{$key}=$info;
		$key=undef;
		$info=undef;
		$string=undef;
		close IN;
		close OUT;

		my %mat=();
		my @mat=();
		open (OUT,">$opt{prefix}\_candidate_mature.fa.new")||die("fail to open $opt{prefix}\_candidate_mature.fa.new.\n");
		open (IN,"<$opt{prefix}\_candidate_mature.fa")||die("fail to open $opt{prefix}\_candidate_mature.fa.\n");
		while (<IN>) {
			chomp;
			if (/^\#/) {
				next;
			}elsif (/^>(\S+)/) {
				if (defined $key) {
					print OUT ">$key\n";
					$string=~tr/Uu/Tt/;
					$string=uc $string;
					print OUT "$string\n";
					$mat{$key}=$string;
					push(@mat,$key);
					$string=undef;
				}
				$key=$1;
				next;
			}
			$string.=$_;
		}
		print OUT ">$key\n";
		$string=~tr/Uu/Tt/;
		$string=uc $string;
		print OUT "$string\n";
		$mat{$key}=$string;
		push (@mat,$key);
		$key=undef;
		$string=undef;
		close IN;
		close OUT;

		if ($opt{target_tool} eq "blast") {
			system "formatdb -i $opt{target_seq}.new -p F";
			system "blastall -p blastn -d $opt{target_seq}.new -i $opt{prefix}\_candidate_mature.fa.new -v 1000 -b 1000 -W 7 -o $opt{prefix}\_target.blast";
			system "perl $programdir/EblastN.pl -i $opt{prefix}\_target.blast -o $opt{prefix}\_target.eblastn";

			open (OUT,">$opt{prefix}\_miRNA.target.blast")||die("Cannot read $opt{prefix}\_miRNA.target.blast.\n");
			open (IN,"<$opt{prefix}\_target.eblastn")||die("Cannot read $opt{prefix}\_target.eblastn.\n");
			while (<IN>) {
				chomp;
				next if (/^Query/);
				my @list=split(/\t/,$_);
				next if (abs($list[5]-$list[4]) != $list[3]-$list[2]);
				if ($list[4]>$list[5]) {
					my $target=substr($target{$list[11]},($list[5]-($list[1]-$list[3]))-1,(($list[4]+($list[2]-1))-($list[5]-($list[1]-$list[3]))+1));
					my $mirna=$mat{$list[0]};
					$mirna=reverse $mirna;
					my @target=split //,$target;
					my @mirna=split //,$mirna;
					my $gu=0;
					my $mis=0;
					my $mis10_11=0;
					my $mis2_12=0;
					my $maxconmis=0;
					my $misstart=undef;
					my $homology_string;
					for(my $i=0;$i<@target;$i++){
						if ($i==0) {
							if ($bp{"$target[$i]$mirna[$i]"}==1) {
								$mis++;
								if (!defined $misstart) {
									$misstart=$i;
									$maxconmis=1 if ($maxconmis<1);
								}else{
									$maxconmis=$i-$misstart+1 if ($maxconmis<$i-$misstart+1);
								}
								$homology_string .=" ";
							}elsif ($bp{"$target[$i]$mirna[$i]"}==0.5) {
								$mis+=0.5;
								$gu++;
								$misstart=undef;
								$homology_string .=".";
							}else {
								$misstart=undef;
								$homology_string .=":";
							}
						}elsif ($i>=1 and $i<=8) {
							if ($bp{"$target[$i]$mirna[$i]"}==1) {
								$mis++;
								$mis2_12++;
								if (!defined $misstart) {
									$misstart=$i;
									$maxconmis=1 if ($maxconmis<1);
								}else{
									$maxconmis=$i-$misstart+1 if ($maxconmis<$i-$misstart+1);
								}
								$homology_string .=" ";
							}elsif ($bp{"$target[$i]$mirna[$i]"}==0.5) {
								$mis+=0.5;
								$mis2_12+=0.5;
								$gu++;
								$misstart=undef;
								$homology_string .=".";
							}else {
								$misstart=undef;
								$homology_string .=":";
							}
						}elsif ($i>=9 and $i<=10) {
							if ($bp{"$target[$i]$mirna[$i]"}==1) {
								$mis++;
								$mis2_12++;
								$mis10_11++;
								if (!defined $misstart) {
									$misstart=$i;
									$maxconmis=1 if ($maxconmis<1);
								}else{
									$maxconmis=$i-$misstart+1 if ($maxconmis<$i-$misstart+1);
								}
								$homology_string .=" ";
							}elsif ($bp{"$target[$i]$mirna[$i]"}==0.5) {
								$mis+=0.5;
								$mis2_12+=0.5;
								$mis10_11+=0.5;
								$gu++;
								$misstart=undef;
								$homology_string .=".";
							}else {
								$misstart=undef;
								$homology_string .=":";
							}
						}elsif ($i==11) {
							if ($bp{"$target[$i]$mirna[$i]"}==1) {
								$mis++;
								$mis2_12++;
								if (!defined $misstart) {
									$misstart=$i;
									$maxconmis=1 if ($maxconmis<1);
								}else{
									$maxconmis=$i-$misstart+1 if ($maxconmis<$i-$misstart+1);
								}
								$homology_string .=" ";
							}elsif ($bp{"$target[$i]$mirna[$i]"}==0.5) {
								$mis+=0.5;
								$mis2_12+=0.5;
								$gu++;
								$misstart=undef;
								$homology_string .=".";
							}else {
								$misstart=undef;
								$homology_string .=":";
							}
						}else {
							if ($bp{"$target[$i]$mirna[$i]"}==1) {
								$mis++;
								if (!defined $misstart) {
									$misstart=$i;
									$maxconmis=1 if ($maxconmis<1);
								}else{
									$maxconmis=$i-$misstart+1 if ($maxconmis<$i-$misstart+1);
								}
								$homology_string .=" ";
							}elsif ($bp{"$target[$i]$mirna[$i]"}==0.5) {
								$mis+=0.5;
								$gu++;
								$misstart=undef;
								$homology_string .=".";
							}else {
								$misstart=undef;
								$homology_string .=":";
							}
						}
					}
					if ($mis10_11==0 and $maxconmis<=2 and $mis2_12<=1 and $mis<=4 and $gu<=3 and $mis-$gu*0.5<=3) {
						$target=~s/T/U/g;
						$mirna=~s/T/U/g;
						printf OUT ">$list[0]\t$list[11]\t$list[12]\n\nTarget%9.0f $target %9.0f\n                $homology_string          \n miRNA%9.0f $mirna %9.0f\n\n",($list[5]-($list[1]-$list[3])),($list[4]+($list[2]-1)),$list[1],1;
					}
				}
			}
			close IN;
			close OUT;
		}elsif ($opt{target_tool} eq "blat") {
				
			system "blat $opt{target_seq}.new $opt{prefix}\_candidate_mature.fa.new -tileSize=8 -oneOff=1 -minMatch=1 -minIdentity=80 -noHead $opt{prefix}\_target.blat";
			
			open (OUT,">$opt{prefix}\_miRNA.target.blat")||die("Cannot read $opt{prefix}\_miRNA.target.blat.\n");
			open (IN,"<$opt{prefix}\_target.blat")||die("Cannot read $opt{prefix}\_target.blat.\n");
			while (<IN>) {
				chomp;
				my @list=split /\t/,$_;
				next if ($list[5]+$list[7]>0);
				if ($list[8] eq "-") {
					my $target=substr($target{$list[13]},($list[15]+1-($list[10]-$list[12])-1),(($list[16]+($list[11]+1-1))-($list[15]+1-($list[10]-$list[12]))+1));
					my $mirna=$mat{$list[9]};
					$mirna=reverse $mirna;
					my @target=split //,$target;
					my @mirna=split //,$mirna;
					my $gu=0;
					my $mis=0;
					my $mis10_11=0;
					my $mis2_12=0;
					my $maxconmis=0;
					my $misstart=undef;
					my $homology_string;
					for(my $i=0;$i<@target;$i++){
						if ($i==0) {
							if ($bp{"$target[$i]$mirna[$i]"}==1) {
								$mis++;
								if (!defined $misstart) {
									$misstart=$i;
									$maxconmis=1 if ($maxconmis<1);
								}else{
									$maxconmis=$i-$misstart+1 if ($maxconmis<$i-$misstart+1);
								}
								$homology_string .=" ";
							}elsif ($bp{"$target[$i]$mirna[$i]"}==0.5) {
								$mis+=0.5;
								$gu++;
								$misstart=undef;
								$homology_string .=".";
							}else {
								$misstart=undef;
								$homology_string .=":";
							}
						}elsif ($i>=1 and $i<=8) {
							if ($bp{"$target[$i]$mirna[$i]"}==1) {
								$mis++;
								$mis2_12++;
								if (!defined $misstart) {
									$misstart=$i;
									$maxconmis=1 if ($maxconmis<1);
								}else{
									$maxconmis=$i-$misstart+1 if ($maxconmis<$i-$misstart+1);
								}
								$homology_string .=" ";
							}elsif ($bp{"$target[$i]$mirna[$i]"}==0.5) {
								$mis+=0.5;
								$mis2_12+=0.5;
								$gu++;
								$misstart=undef;
								$homology_string .=".";
							}else {
								$misstart=undef;
								$homology_string .=":";
							}
						}elsif ($i>=9 and $i<=10) {
							if ($bp{"$target[$i]$mirna[$i]"}==1) {
								$mis++;
								$mis2_12++;
								$mis10_11++;
								if (!defined $misstart) {
									$misstart=$i;
									$maxconmis=1 if ($maxconmis<1);
								}else{
									$maxconmis=$i-$misstart+1 if ($maxconmis<$i-$misstart+1);
								}
								$homology_string .=" ";
							}elsif ($bp{"$target[$i]$mirna[$i]"}==0.5) {
								$mis+=0.5;
								$mis2_12+=0.5;
								$mis10_11+=0.5;
								$gu++;
								$misstart=undef;
								$homology_string .=".";
							}else {
								$misstart=undef;
								$homology_string .=":";
							}
						}elsif ($i==11) {
							if ($bp{"$target[$i]$mirna[$i]"}==1) {
								$mis++;
								$mis2_12++;
								if (!defined $misstart) {
									$misstart=$i;
									$maxconmis=1 if ($maxconmis<1);
								}else{
									$maxconmis=$i-$misstart+1 if ($maxconmis<$i-$misstart+1);
								}
								$homology_string .=" ";
							}elsif ($bp{"$target[$i]$mirna[$i]"}==0.5) {
								$mis+=0.5;
								$mis2_12+=0.5;
								$gu++;
								$misstart=undef;
								$homology_string .=".";
							}else {
								$misstart=undef;
								$homology_string .=":";
							}
						}else {
							if ($bp{"$target[$i]$mirna[$i]"}==1) {
								$mis++;
								if (!defined $misstart) {
									$misstart=$i;
									$maxconmis=1 if ($maxconmis<1);
								}else{
									$maxconmis=$i-$misstart+1 if ($maxconmis<$i-$misstart+1);
								}
								$homology_string .=" ";
							}elsif ($bp{"$target[$i]$mirna[$i]"}==0.5) {
								$mis+=0.5;
								$gu++;
								$misstart=undef;
								$homology_string .=".";
							}else {
								$misstart=undef;
								$homology_string .=":";
							}
						}
					}
					if ($mis10_11==0 and $maxconmis<=2 and $mis2_12<=1 and $mis<=4 and $gu<=3 and $mis-$gu*0.5<=3) {
						$target=~s/T/U/g;
						$mirna=~s/T/U/g;
						printf OUT ">$list[9]\t$list[13]\t$tinfo{$list[13]}\n\nTarget%9.0f $target %9.0f\n                $homology_string          \n miRNA%9.0f $mirna %9.0f\n\n",($list[15]+1-($list[10]-$list[12])),($list[16]+($list[11]+1-1)),$list[10],1;
					}
				}
			}
		}elsif ($opt{target_tool} eq "targetfinder") {

			open (OUT,">$opt{prefix}\_miRNA.target.targetfinder")||die("Cannot read $opt{prefix}\_miRNA.target.targetfinder.\n");
			foreach $key (@mat) {
				system "$programdir/targetfinder.pl -s $mat{$key} -d $opt{target_seq}.new -q $key > $key.out";
				open (RESULT,"<$key.out")||die("Cannot read $key.out.\n");
				while (<RESULT>) {
					chomp;
					if (/^(query=.*)/) {
						print OUT ">$_\n";
					}else {
						print OUT "$_\n";
					}
				}
				close RESULT;
				system "rm $key.out";
			}
			close OUT;
		}elsif ($opt{target_tool} eq "miranda") {
			open (OUT,">$opt{prefix}\_miRNA.target.miranda")||die("Cannot read $opt{prefix}\_miRNA.target.miranda.\n");
			print OUT "#Query\tTarget\tTot Score\tTot Energy\tMax Score\tMax Energy\tStrand\tLen1\tLen2\tPositions\n";
			open (IN,"<$opt{prefix}\_candidate_mature.fa.new")||die("fail to open $opt{prefix}\_candidate_mature.fa.new.\n");
			while (<IN>) {
				chomp;
				if (/^\#/) {
					next;
				}elsif (/^>(\S+)/) {
					if (defined $key) {
						open (RECORD,">$key.seq")||die("Cannot write to $key.seq.\n");
						print RECORD ">$key\n$string\n";
						close RECORD;

						system "miranda $key.seq $opt{target_seq}.new -out $key.out";
						open (RESULT,"<$key.out")||die("fail to open $key.out.\n");
						my $targetinfo;
						while (<RESULT>) {
							chomp;
							if (/Read Sequence:(.*)/) {
								$targetinfo=$1;
							}elsif (/>>(.*)/) {
								print OUT ">$key\t$targetinfo\n$1\n";
							}
						}
						close RESULT;
						system "rm $key.seq";
						system "rm $key.out";
						$string=undef;
					}
					$key=$1;
					next;
				}
				$string.=$_;
			}
			open (RECORD,">$key.seq")||die("Cannot write to $key.seq.\n");
			print RECORD ">$key\n$string\n";
			close RECORD;
			system "miranda $key.seq $opt{target_seq}.new -out $key.out";
			open (RESULT,"<$key.out")||die("Cannot read $key.out.\n");
			my $targetinfo;
			while (<RESULT>) {
				chomp;
				if (/Read Sequence:(.*)/) {
					$targetinfo=$1;
				}elsif (/>>(.*)/) {
					print OUT ">$key\t$targetinfo\n$1\n";
				}
			}
			close RESULT;
			system "rm $key.seq";
			system "rm $key.out";
			$string=undef;
			close IN;
			close OUT;
		}elsif ($opt{target_tool} eq "RNAhybrid") {
			system "RNAhybrid -t $opt{target_seq}.new -q $opt{prefix}\_candidate_mature.fa.new -s $opt{dataset} > $opt{prefix}\_miRNA.target.RNAhybrid";
		}

		$time=localtime;
		print LOG "Target for candidate miRNA has been Identified: $time.\n";
	}
}

$time=localtime;
print LOG "The miRNA identification has been finished: $time.\n\n";
close LOG;

sub base_content {
	my $seq=shift;
	my $A=($seq=~tr/Aa/Aa/);
	my $C=($seq=~tr/Cc/Cc/);
	my $G=($seq=~tr/Gg/Gg/);
	my $T=($seq=~tr/TtUu/TtUu/);
	my $N=($seq=~tr/Nn/Nn/);
	my $a_content=sprintf "%.2f",($A)/($G+$C+$A+$T+$N)*100;
	my $u_content=sprintf "%.2f",($T)/($G+$C+$A+$T+$N)*100;
	my $g_content=sprintf "%.2f",($G)/($G+$C+$A+$T+$N)*100;
	my $c_content=sprintf "%.2f",($C)/($G+$C+$A+$T+$N)*100;
	my $au_content=sprintf "%.2f",($A+$T)/($G+$C+$A+$T+$N)*100;
	my $gc_content=sprintf "%.2f",($G+$C)/($G+$C+$A+$T+$N)*100;
	return ($a_content,$u_content,$g_content,$c_content,$au_content,$gc_content);
}

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

sub low_complexity {
	my $str=shift;
	my $type=0;
	my $max=0;
	my $num;
	my $base;
	my @str=split //,uc($str);
	if (@str>0) {
		foreach my $key (@str) {
			if (!defined $base) {
				$base=$key;
				$num=1;
			}else {
				if ($key eq $base) {
					$num++;
				}else {
					$max=$num if ($max<$num);
					$base=$key;
					$num=1;
				}
			}
		}
		my $A=($str=~tr/A/A/);
		my $C=($str=~tr/C/C/);
		my $G=($str=~tr/G/G/);
		my $T=($str=~tr/T/T/);
		if ($max<=5 and ($A+$C)/(length($str))<0.9 and ($A+$G)/(length($str))<0.9 and ($A+$T)/(length($str))<0.9 and ($C+$G)/(length($str))<0.9 and ($C+$T)/(length($str))<0.9 and ($G+$T)/(length($str))<0.9) {
			$type=1;
		}
	}
	return $type;
}
