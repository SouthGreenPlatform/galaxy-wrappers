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

	Usage: perl candidate.pl --mat_fa|-mf <*_mature.fa> --mat_gff|-mg <*_mature.gff> --pre_fa|-pf <*_pre.fa> --pre_gff|-pg <*_pre.gff>
	
	The following options are necessary.
	--mat_fa|-mf	<*_mature.fa>	:candidate mature sequence file in fasta format
	--mat_gff|-mg	<*_mature.gff>	:candidate mature gff file in gff2 format
	--pre_fa|-pf	<*_pre.fa>	:candidate precursor sequence file in fasta format
	--pre_gff|-pg	<*_pre.gff>	:candidate precursor gff file in gff2 format

	The following options are optional.
	--output_dir|-o		<output_dir>	:default is sample under the current directory, you should change it if you have more than one sample to distinguish different samples
	--prefix|-p		<prefix=MI>	:prefix for file and miRNA ID, if you want to compare different samples, you should use different prefix, the prefix should be end with letters, default is MI
	--minbasepair|-mbp	<minbasepair=16>	:minimum base-pairs between mature and star of miRNA comparison, default is 16
	--maxbasebulge|-mbb	<maxbasebulge=3>	:maximum base bulge in mature and star miRNA comparison, default is 3
	--maxunpairbase|-mub	<maxunpairbase=6>	:maximum unpair base number in mature or star miRNA region, default is 6
	--matoverrate|-mor	<matoverrate=60>	:minimum redundant mature sequence overlap rate percent(%), default is 60
	--help|-h				:print the usage information

USAGE

#Gather input
&GetOptions(
	"mat_fa|mf=s"		=>\$opt{mat_fa},
	"mat_gff|mg=s"		=>\$opt{mat_gff},
	"pre_fa|pf=s"		=>\$opt{pre_fa},
	"pre_gff|pg=s"		=>\$opt{pre_gff},
	"output_dir|o=s"	=>\$opt{output_dir},
	"prefix|p=s"		=>\$opt{prefix},
	"minbasepair|mbp=i"	=>\$opt{minbasepair},
	"maxbasebulge|mbb=i"=>\$opt{maxbasebulge},
	"maxunpairbase|mub=i"=>\$opt{maxunpairbase},
	"matoverrate|mor=i"	=>\$opt{matoverrate},
	"help|h"	=>\$opt{help},
);

#Verify input
if (!defined $opt{mat_fa} or !defined $opt{mat_gff} or !defined $opt{pre_fa} or !defined $opt{pre_gff} or defined $opt{help}) {
	die "$usage\n";
}

#Default parameters
$opt{prefix}	= "MI" unless defined $opt{prefix};
$opt{minbasepair}	= 16 unless defined $opt{minbasepair};
$opt{maxbasebulge}	= 3 unless defined $opt{maxbasebulge};
$opt{maxunpairbase}	= 6 unless defined $opt{maxunpairbase};
$opt{matoverrate}	= 60 unless defined $opt{matoverrate};

# Absolute path
$opt{mat_fa}      = abs_path $opt{mat_fa} if defined $opt{mat_fa};
$opt{mat_gff}     = abs_path $opt{mat_gff} if defined $opt{mat_gff};
$opt{pre_fa}      = abs_path $opt{pre_fa} if defined $opt{pre_fa};
$opt{pre_gff}     = abs_path $opt{pre_gff} if defined $opt{pre_gff};

#Output directory
$opt{output_dir} = $ENV{'PWD'}."/sample"  unless defined $opt{output_dir};
$opt{output_dir} = abs_path $opt{output_dir} if defined $opt{output_dir};
`mkdir -p $opt{output_dir}` unless -d $opt{output_dir};
chdir $opt{output_dir};

#Build commond line
my $cmdline = $0;
$cmdline .= " -mf $opt{mat_fa}" if defined $opt{mat_fa};
$cmdline .= " -mg $opt{mat_gff}" if defined $opt{mat_gff};
$cmdline .= " -pf $opt{pre_fa}" if defined $opt{pre_fa};
$cmdline .= " -pg $opt{pre_gff}" if defined $opt{pre_gff};
$cmdline .= " -o $opt{output_dir}"  if defined $opt{output_dir};
$cmdline .= " -p $opt{prefix}" if defined $opt{prefix};
$cmdline .= " -mbp $opt{minbasepair}" if defined $opt{minbasepair};
$cmdline .= " -mbb $opt{maxbasebulge}" if defined $opt{maxbasebulge};
$cmdline .= " -mub $opt{maxunpairbase}" if defined $opt{maxunpairbase};
$cmdline .= " -mor $opt{matoverrate}" if defined $opt{matoverrate};

print `date`, "\n>>> candidate.pl <<<\n\n";
print $cmdline,"\n\n";

my ($time,);

#Do RNAfold for obtained precursor sequence
print "Do RNAfold for obtained precursor sequence.\n";
system "less $opt{pre_fa} | RNAfold --noconv --noPS > $opt{prefix}\_pre.RNAfold";

$time=localtime;
print "The RNAfold for obtained precursor sequence has been done: $time.\n\n";

#Filter precursor based on RNAfold result
print "Filter precursor based on RNAfold result.\n";
my ($key,%mat,);
open (MAT,"<$opt{mat_fa}")||die("Cannot read $opt{mat_fa}.\n");
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
print "The RNAfold result of precursor has been filtered: $time.\n\n";

#Extract hairpin sequence and create hairpin gff file
print "Extract hairpin sequence and create hairpin gff file for filtered precursor.\n";
my (%pre,%pos,%id,);
open (PRE,"<$opt{pre_fa}")||die("Cannot read $opt{pre_fa}.\n");
while (<PRE>) {
	chomp;
	if (/^>(\S+)/) {
		$key=$1;
	}else{
		$pre{$key}.=$_;
	}
}
close PRE;

open (IN,"<$opt{pre_gff}")||die("Cannot read $opt{pre_gff}.\n");
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
			print GFF "$pos[0]\tMIP\tmiRNA\t$pre_start\t$pre_end\t\.\t$pos[6]\t\.\tACC=\"$feature[0]\"; ID=\"$feature[1]\"; EXPRESS=\"$feature[2]\";\n";
		}else {
			my $pre_end=$pos[4]-$list[5]+1;
			my $pre_start=$pos[4]-$list[8]+1;
			print GFF "$pos[0]\tMIP\tmiRNA\t$pre_start\t$pre_end\t\.\t$pos[6]\t\.\tACC=\"$feature[0]\"; ID=\"$feature[1]\"; EXPRESS=\"$feature[2]\";\n";
		}
	}else {
		my $seq=substr ($pre{$list[0]},$list[7]-1,$list[6]-$list[7]+1);
		print FASTA ">$list[0]\n$seq\n";
		if ($pos[6] eq "+") {
			my $pre_start=$pos[3]+$list[7]-1;
			my $pre_end=$pos[3]+$list[6]-1;
			print GFF "$pos[0]\tMIP\tmiRNA\t$pre_start\t$pre_end\t\.\t$pos[6]\t\.\tACC=\"$feature[0]\"; ID=\"$feature[1]\"; EXPRESS=\"$feature[2]\";\n";
		}else {
			my $pre_end=$pos[4]-$list[7]+1;
			my $pre_start=$pos[4]-$list[6]+1;
			print GFF "$pos[0]\tMIP\tmiRNA\t$pre_start\t$pre_end\t\.\t$pos[6]\t\.\tACC=\"$feature[0]\"; ID=\"$feature[1]\"; EXPRESS=\"$feature[2]\";\n";
		}
	}
}
close FILTER;
close FASTA;
close GFF;

$time=localtime;
print "The hairpin sequences and gff file have been extracted: $time.\n\n";

#Filter redundant mature record from mature and hairpin fasta and gff file
print "Filter redundant mature record from mature and hairpin fasta and gff files.\n";

%pre=();%pos=();
my (%repeat,);
open (IN,"<$opt{mat_gff}")||die("Cannot read $opt{mat_gff}.\n");
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
print "The redundant mature record in mature gff file has been filtered: $time.\n";

open (FASTA,">$opt{prefix}\_mature_unique.fa")||die("Cannot write to $opt{prefix}\_mature_unique.fa file.\n");
open (IN,"<$opt{mat_fa}")||die("Cannot read $opt{mat_fa}.\n");
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
print "The redundant mature sequence in mature fasta file has been filtered: $time.\n";

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
print "The redundant hairpin record in hairpin gff file has been filtered: $time.\n";

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
print "The redundant hairpin sequence in hairpin fasta file has been filtered: $time.\n";
%id=();

#Do RNAfold for unique hairpin sequence
print "Do RNAfold for unique hairpin sequences.\n";
system "less $opt{prefix}\_hairpin_unique.fa | RNAfold --noconv --noPS > $opt{prefix}\_hairpin_unique.RNAfold";

$time=localtime;
print "The RNAfold for unique hairpin sequences has been done: $time.\n\n";
