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

	Usage: perl $program --conserveseq|-cs <*.fa> --candidate_mat|-cm <*_candidate_mature.fa>
	
	The following options are necessary.
	--conserveseq|-cs	<*.fa>	:The mature sequence of conserve miRNA. This file format must like the mature sequence file in miRBase.
	--candidate_mat|-cm	<*_candidate_mature.fa>	:The candidate mature miRNA sequence file in fasta format

	The following options are optional.
	--output_dir|-o		<output_dir>	:default is sample under the current directory, you should change it if you have more than one sample to distinguish different samples
	--prefix|-p		<prefix=MI>	:prefix for file and miRNA ID, if you want to compare different samples, you should use different prefix, the prefix should be end with letters, default is MI
	--mismatch|-m		<mismatch=3>	:maximum mismatch value for mapping result (including indels), default is 3
	--identity|-i		<identity=90>	:minimum identity percent(%) (in mapped regions), default is 90
	--minimum|-min		<minimum=18>	:minimum miRNA length, default is 18
	--maximum|-max		<maximum=25>	:maximum miRNA length, default is 25
	--alignsoft|-as	<blast|blat>	:in default, blast was used.
	--triletterabbr|-tla	<triletterabbr=new>	:three letter abbreviation for studied species, it used to compare miRNA conservation in multiple species, default is new
	--help|-h				:print the usage information

USAGE

#Gather input
&GetOptions(
	"conserveseq|cs=s"	=>\$opt{conserveseq},
	"candidate_mat|cm=s"	=>\$opt{candidate_mat},
	"output_dir|o=s"	=>\$opt{output_dir},
	"prefix|p=s"		=>\$opt{prefix},
	"mismatch|m=i"	=>\$opt{mismatch},
	"identity|i=i"	=>\$opt{identity},
	"minimum|min=s"	=>\$opt{minimum},
	"maximum|max=s"	=>\$opt{maximum},
	"alignsoft|as=s"	=>\$opt{alignsoft},
	"triletterabbr|tla=s"	=>\$opt{triletterabbr},
	"help|h"	=>\$opt{help},
);

#Verify input
if (!defined $opt{conserveseq} or !defined $opt{candidate_mat} or defined $opt{help}) {
	die "$usage\n";
}

#Default parameters
$opt{prefix}	= "MI" unless defined $opt{prefix};
$opt{mismatch}	= 3 unless defined $opt{mismatch};
$opt{identity}	= 90 unless defined $opt{identity};
$opt{minimum}	= 18 unless defined $opt{minimum};
$opt{maximum}	= 25 unless defined $opt{maximum};
$opt{alignsoft}	= "blast" unless defined $opt{alignsoft};
$opt{triletterabbr}	= "new" unless defined $opt{triletterabbr};

# Absolute path
$opt{conserveseq}  = abs_path $opt{conserveseq} if defined $opt{conserveseq};
$opt{candidate_mat}= abs_path $opt{candidate_mat} if defined $opt{candidate_mat};

#Output directory
$opt{output_dir} = $ENV{'PWD'}."/sample"  unless defined $opt{output_dir};
$opt{output_dir} = abs_path $opt{output_dir} if defined $opt{output_dir};
`mkdir -p $opt{output_dir}` unless -d $opt{output_dir};
chdir $opt{output_dir};

#Build commond line
my $cmdline = $program;
$cmdline .= " -cs $opt{conserveseq}" if defined $opt{conserveseq};
$cmdline .= " -cm $opt{candidate_mat}" if defined $opt{candidate_mat};
$cmdline .= " -o $opt{output_dir}"  if defined $opt{output_dir};
$cmdline .= " -p $opt{prefix}" if defined $opt{prefix};
$cmdline .= " -m $opt{mismatch}" if defined $opt{mismatch};
$cmdline .= " -i $opt{identity}" if defined $opt{identity};
$cmdline .= " -min $opt{minimum}" if defined $opt{minimum};
$cmdline .= " -max $opt{maximum}" if defined $opt{maximum};
$cmdline .= " -e $opt{alignsoft}" if defined $opt{alignsoft};
$cmdline .= " -tla $opt{triletterabbr}" if defined $opt{triletterabbr};

print `date`, "\n>>> $program <<<\n\n";
print $cmdline,"\n\n";

my ($key,$string,$time,);

#Identify conserve miRNA from candidate miRNA
print "Identify conserve miRNA from candidate miRNA.\n";

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

my %id=();
my %sample=();
open (IN,"<$opt{candidate_mat}")||die("fail to open $opt{candidate_mat}.\n");
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
	system "blastall -p blastn -d $opt{conserveseq}.new -i $opt{candidate_mat} -v 1000 -b 1000 -W 7 -o $opt{prefix}\_candidate_mature.blast";
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

	system "blat $opt{conserveseq}.new $opt{candidate_mat} -tileSize=8 -oneOff=1 -minMatch=1 -minIdentity=80 -noHead $opt{prefix}\_candidate_mature.blat";
	
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
print "The conserve miRNA in candidate miRNA has been identified: $time.\n";

print "Do miRNA compare for selected multiple species and this sample.\n";

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
print "The miRNA comparison for selected multiple species has been finished: $time.\n\n";

