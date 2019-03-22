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

my $usage=<<USAGE; #******* Instruction of this program *********#

	Author: Wanfei Liu & Chengqi Xin
	Email: <liuwf\@big.ac.cn> & <xinchq\@big.ac.cn>
	Date: Jun 8, 2013
	Version: $version

	Usage: perl $program --target_seq|-ts <*.fa> --candidate_mat|-cm <*_candidate_mature.fa>
	
	The following options are necessary.
	--target_seq|-ts	<*.fa>	:The target sequence file. If you assign this parameter, the program will do target prediction for candidate miRNA.
	--candidate_mat|-cm	<*_candidate_mature.fa>	:The candidate mature miRNA sequence file in fasta format

	The following options are optional.
	--output_dir|-o		<output_dir>	:default is sample under the current directory, you should change it if you have more than one sample to distinguish different samples
	--prefix|-p		<prefix=MI>	:prefix for file and miRNA ID, if you want to compare different samples, you should use different prefix, the prefix should be end with letters, default is MI
	--alignsoft|-as		<blast|blat>	:in default, blast was used.
	--target_tool|-tt	<*>	:The program for target prediction (blast, blat, targetfinder, miranda and RNAhybrid). For plant, you can identify target genes by our rule set according to blast or blat alignment (in default) or by targetfinder; for animal, you can use miranda or RNAhybrid.
	--dataset|-ds	<dataset=3utr_fly|3utr_worm|3utr_human>	:three data set name in RNAhybird program for target prediction, in default is 3utr_human. This parameter is only need when you use RNAhybrid to predict target
	--help|-h				:print the usage information

USAGE

#Gather input
&GetOptions(
	"target_seq|ts=s"	=>\$opt{target_seq},
	"candidate_mat|cm=s"	=>\$opt{candidate_mat},
	"output_dir|o=s"	=>\$opt{output_dir},
	"prefix|p=s"		=>\$opt{prefix},
	"alignsoft|as=s"	=>\$opt{alignsoft},
	"target_tool|tt=s"	=>\$opt{target_tool},
	"dataset|ds=s"	=>\$opt{dataset},
	"help|h"	=>\$opt{help},
);

#Verify input
if (!defined $opt{target_seq} or !defined $opt{candidate_mat} or defined $opt{help}) {
	die "$usage\n";
}

#Default parameters
$opt{prefix}	= "MI" unless defined $opt{prefix};
$opt{alignsoft}	= "blast" unless defined $opt{alignsoft};
$opt{target_tool}	= $opt{alignsoft} unless defined $opt{target_tool};
$opt{dataset}	= '3utr_human' unless defined $opt{dataset};

# Absolute path
$opt{target_seq}   = abs_path $opt{target_seq} if defined $opt{target_seq};
$opt{candidate_mat}= abs_path $opt{candidate_mat} if defined $opt{candidate_mat};

#Output directory
$opt{output_dir} = $ENV{'PWD'}."/sample"  unless defined $opt{output_dir};
$opt{output_dir} = abs_path $opt{output_dir} if defined $opt{output_dir};
`mkdir -p $opt{output_dir}` unless -d $opt{output_dir};
chdir $opt{output_dir};

#Build commond line
my $cmdline = $program;
$cmdline .= " -ts $opt{target_seq}" if defined $opt{target_seq};
$cmdline .= " -cm $opt{candidate_mat}" if defined $opt{candidate_mat};
$cmdline .= " -o $opt{output_dir}"  if defined $opt{output_dir};
$cmdline .= " -p $opt{prefix}" if defined $opt{prefix};
$cmdline .= " -as $opt{alignsoft}" if defined $opt{alignsoft};
$cmdline .= " -tt $opt{target_tool}" if defined $opt{target_tool};
$cmdline .= " -ds $opt{dataset}" if defined $opt{dataset};

print `date`, "\n>>> $program <<<\n\n";
print $cmdline,"\n\n";

my ($time,);

#Identify target for candidate miRNA
print "Identify target for candidate miRNA.\n";

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

my $key=undef;
my $info=undef;
my $string=undef;
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
open (OUT,">$opt{candidate_mat}.new")||die("fail to open $opt{candidate_mat}.new.\n");
open (IN,"<$opt{candidate_mat}")||die("fail to open $opt{candidate_mat}.\n");
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
			push (@mat,$key);
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
	system "blastall -p blastn -d $opt{target_seq}.new -i $opt{candidate_mat}.new -v 1000 -b 1000 -W 7 -o $opt{prefix}\_target.blast";
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
		
	system "blat $opt{target_seq}.new $opt{candidate_mat}.new -tileSize=8 -oneOff=1 -minMatch=1 -minIdentity=80 -noHead $opt{prefix}\_target.blat";
	
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
	open (IN,"<$opt{candidate_mat}.new")||die("fail to open $opt{candidate_mat}.new.\n");
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
	
	system "RNAhybrid -t $opt{target_seq}.new -q $opt{candidate_mat}.new -s $opt{dataset} > $opt{prefix}\_miRNA.target.RNAhybrid";
}

$time=localtime;
print "Target for candidate miRNA has been Identified: $time.\n";
