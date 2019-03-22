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
	Date: Jun 8, 2013
	Version: $version

	Usage: perl $program --read_cluster|-rc <*_read.cluster> { --filter_gff|-fg <*.gff> | --filter_fa|-ff <*.fa> & --genome|-g <*.fa> }
	
	The following options are necessary.
	--read_cluster|-rc	<*_read.cluster>	:read cluster result
	--filter_gff|-fg		<*.gff>		:the filter record file in gff2 format, all clusters overlap with these records will be removed (please do not including intron regions)
	--filter_fa|-ff		<*.fa>		:the filter record file in fasta format, all clusters overlap with these records will be removed (please do not including intron regions)
	--genome|-g	<*.fa>	:genome sequence file in fasta format

	The following options are optional.
	--output_dir|-o		<output_dir>	:default is sample under the current directory, you should change it if you have more than one sample to distinguish different samples
	--prefix|-p		<prefix=MI>	:prefix for file and miRNA ID, if you want to compare different samples, you should use different prefix, the prefix should be end with letters, default is MI
	--strand_specific|-ss			:if you assign this parameter, it means the data is strand specific, if not, as default, it means the data is not strand specific
	--mismatch|-m		<mismatch=3>	:maximum mismatch value for mapping result (including indels), default is 3
	--identity|-i		<identity=90>	:minimum identity percent(%) (in mapped regions), default is 90	
	--minimum|-min		<minimum=18>	:minimum miRNA length, default is 18
	--maximum|-max		<maximum=25>	:maximum miRNA length, default is 25
	--alignsoft|-as	<blast|blat>	:in default, blast was used.
	--filter_rate|-fr		<filter_rate=50>	:minimum filter overlap rate percent(%) between filter region and read cluster, default is 50
	--help|-h				:print the usage information

USAGE

#Gather input
&GetOptions(
	"read_cluster|rc=s"	=>\$opt{read_cluster},
	"filter_gff|fg=s"	=>\$opt{filter_gff},
	"filter_fa|ff=s"	=>\$opt{filter_fa},
	"genome|g=s"		=>\$opt{genome},
	"output_dir|o=s"	=>\$opt{output_dir},
	"prefix|p=s"		=>\$opt{prefix},
	"strand_specific|ss"=>\$opt{strand_specific},
	"mismatch|m=i"		=>\$opt{mismatch},
	"identity|i=i"		=>\$opt{identity},	
	"minimum|min=i"		=>\$opt{minimum},
	"maximum|max=i"		=>\$opt{maximum},
	"alignsoft|as=s"	=>\$opt{alignsoft},
	"filter_rate|fr=i"	=>\$opt{filter_rate},
	"help|h"	=>\$opt{help},
);

#Verify input
if (defined $opt{filter_gff}) {
	if (!defined $opt{read_cluster} or defined $opt{help}) {
		die "$usage\n";
	}
}elsif (defined $opt{filter_fa} or defined $opt{genome}) {
	if (!defined $opt{read_cluster} or !defined $opt{filter_fa} or !defined $opt{genome} or defined $opt{help}) {
		die "$usage\n";
	}
}else {
	die "$usage\n";
}

#Default parameters
$opt{prefix}	= "MI" unless defined $opt{prefix};
$opt{mismatch}	= 3 unless defined $opt{mismatch};
$opt{identity}	= 90 unless defined $opt{identity};
$opt{minimum}	= 18 unless defined $opt{minimum};
$opt{maximum}	= 25 unless defined $opt{maximum};
$opt{alignsoft}	= "blast" unless defined $opt{alignsoft};
$opt{filter_rate}	= 50 unless defined $opt{filter_rate};

# Absolute path
$opt{read_cluster}  = abs_path $opt{read_cluster} if defined $opt{read_cluster};
$opt{filter_gff}  = abs_path $opt{filter_gff} if defined $opt{filter_gff};
$opt{filter_fa}   = abs_path $opt{filter_fa} if defined $opt{filter_fa};
$opt{genome}      = abs_path $opt{genome} if defined $opt{genome};

#Output directory
$opt{output_dir} = $ENV{'PWD'}."/sample"  unless defined $opt{output_dir};
$opt{output_dir} = abs_path $opt{output_dir} if defined $opt{output_dir};
`mkdir -p $opt{output_dir}` unless -d $opt{output_dir};
chdir $opt{output_dir};

#Build commond line
my $cmdline = $program;
$cmdline .= " -rc $opt{read_cluster}"  if defined $opt{read_cluster};
$cmdline .= " -fg $opt{filter_gff}"  if defined $opt{filter_gff};
$cmdline .= " -ff $opt{filter_fa}"  if defined $opt{filter_fa};
$cmdline .= " -g $opt{genome}" if defined $opt{genome};
$cmdline .= " -o $opt{output_dir}"  if defined $opt{output_dir};
$cmdline .= " -p $opt{prefix}" if defined $opt{prefix};
$cmdline .= " -ss $opt{strand_specific}" if defined $opt{strand_specific};
$cmdline .= " -m $opt{mismatch}" if defined $opt{mismatch};
$cmdline .= " -i $opt{identity}" if defined $opt{identity};
$cmdline .= " -min $opt{minimum}" if defined $opt{minimum};
$cmdline .= " -max $opt{maximum}" if defined $opt{maximum};
$cmdline .= " -as $opt{alignsoft}" if defined $opt{alignsoft};
$cmdline .= " -fr $opt{filter_rate}"  if defined $opt{filter_rate};

print `date`, "\n>>> $program <<<\n\n";
print $cmdline,"\n\n";

my ($chr,$string,%gen,%len,%id,%hash,$time,);

if (defined $opt{filter_fa}) {
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
}

if (defined $opt{filter_gff}) {
	#Read $opt{prefix}_read.cluster file.\n";
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
	print "The $opt{read_cluster} file has been readed: $time.\n\n";

	#Filter clusters which belong to regions in filter file (gff2 format)
	print "Filter clusters which belong to regions in filter file.\n";

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

	my $out="$opt{prefix}\_".(basename($opt{filter_gff})).".read";
	
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
				COM:
				if ($hend<$fstart) {
					$h++;
					if (defined $hash[$h]) {
						($hstart,$hend)=split /\t/,$hash[$h];
						goto COM;
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
						goto COM;
					}else {
						print FILTER "$filter{$chr}{$strand}{$filter[$f]}[1] EXPRESS$connector\"$express\"\n";
						next;
					}
				}elsif ($hstart>$fend) {
					print FILTER "$filter{$chr}{$strand}{$filter[$f]}[1] EXPRESS$connector\"$express\"\n";
					if (defined $min_overlap_h) {
						$h=$min_overlap_h;
						($hstart,$hend)=split /\t/,$hash[$h];
						$min_overlap_h=undef;
					}
					next;
				}
			}
		}
	}
	close FILTER;

	$time=localtime;
	print "The clusters which belong to regions in filter file have been filtered: $time.\n\n";

	#Print filtered cluster result
	print "Print filtered cluster result.\n";
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
	print "The filtered clusters result have been printed: $time.\n\n";

}elsif (defined $opt{filter_fa}) {

	#Read $opt{prefix}_read.cluster file.\n";
	print "Read $opt{read_cluster} file.\n";

	open (OUT,">$opt{read_cluster}.fa")||die("Cannot write to $opt{read_cluster}.fa file.\n");
	open (IN,"<$opt{read_cluster}")||die("Cannot read $opt{read_cluster}.\n");
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
	print "The $opt{read_cluster} file has been readed: $time.\n\n";

	#Filter clusters which belong to regions in filter file (gff2 format)
	print "Filter clusters which belong to regions in filter file.\n";

	if ($opt{alignsoft} eq "blast") {
		system "formatdb -i $opt{filter_fa} -p F";
		system "blastall -p blastn -d $opt{filter_fa} -i $opt{read_cluster}.fa -v 1000 -b 1000 -W 7 -o $opt{read_cluster}.blast";
		system "perl $programdir/EblastN.pl -i $opt{read_cluster}.blast -o $opt{read_cluster}.eblastn";

		my %filter=();
		open (IN,"<$opt{read_cluster}.eblastn")||die("Cannot read $opt{read_cluster}.eblastn.\n");
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

		my $out="$opt{prefix}\_".(basename($opt{filter_fa})).".read";

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
		print "The clusters which belong to regions in filter file have been filtered: $time.\n\n";

		#Print filtered cluster result
		print "Print filtered cluster result.\n";
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
		print "The filtered clusters result have been printed: $time.\n\n";

	}elsif ($opt{alignsoft} eq "blat") {
		system "blat $opt{filter_fa} $opt{read_cluster}.fa -tileSize=8 -oneOff=1 -minMatch=1 -minIdentity=80 -noHead $opt{read_cluster}.blat";
		
		my %filter=();
		open (IN,"<$opt{read_cluster}.blat")||die("Cannot read $opt{read_cluster}.blat.\n");
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

		my $out="$opt{prefix}\_".(basename($opt{filter_fa})).".read";

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
		print "The clusters which belong to regions in filter file have been filtered: $time.\n\n";

		#Print filtered cluster result
		print "Print filtered cluster result.\n";
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
		print "The filtered clusters result have been printed: $time.\n\n";
	}
}else {
	print "Does not assign filter file.\n";
	system "cp $opt{read_cluster} $opt{read_cluster}_filter";
	$time=localtime;
	print "Just do a copy named $opt{read_cluster}_filter for $opt{read_cluster} file: $time.\n\n";

	open (STAT,">>$opt{prefix}.statistics")||die("Cannot write to $opt{prefix}.statistics.\n");
	print STAT "Filtered_cluster_number\t0\n";
	print STAT "Filtered_cluster_read\t0\n\n";
	close STAT;
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
	my $bin=unpack("B*",$flag);
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
