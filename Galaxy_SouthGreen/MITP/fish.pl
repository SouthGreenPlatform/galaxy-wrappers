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

	This program can fishing in one file according to bait in another file among fasta file, gff file and list file.
	
	Usage: perl fish.pl --tbait|-tb <fa|gff|list> --tfish|-tf <fa|gff|list> --bait|-b <*> --fish|-f <*> --output|-o <*>
	
	The following options are necessary.
	--tbait|-tb	<fa|gff|list>	:type of bait file (fa, gff or list)
	--tfish|-tf	<fa|gff|list>	:type of fish file (fa, gff or list)
	--bait|-b		:bait file
	--fish|-f		:fish file
	--output|-o		:output file

	The following options are optional.
	--cbait|-cb		:bait colum in bait file, default is 1 for fa and list file, 9 for gff file
	--cfish|-cf		:fish colum in fish file, default is 1 for fa and list file, 9 for gff file
	--contrary|-c		:fish those not in the bait file
	--help|-h		:print the usage information

	Comment: This program can deal with all file formats used in MIP program (fa, gff and list). list file is just like excel tables which have several colums and rows. The program gets the baits (key that is universal in both bait and fish file) from the bait file at first, and then searches the baits in the fish file. If -contrary is specifed, retrive those items which are not exist in the bait file.

USAGE

#Gather input
&GetOptions(
	"tbait|tb=s"	=>\$opt{tbait},
	"tfish|tf=s"	=>\$opt{tfish},
	"bait|b=s"	=>\$opt{bait},
	"fish|f=s"	=>\$opt{fish},
	"output|o=s"	=>\$opt{output},
	"cbait|cb=i"	=>\$opt{cbait},
	"cfish|cf=i"	=>\$opt{cfish},
	"contrary|c"	=>\$opt{contrary},
	"help|h"	=>\$opt{help},
);

#Verify input
if (!defined $opt{tbait} or !defined $opt{tfish} or !defined $opt{bait} or !defined $opt{fish} or !defined $opt{output} or defined $opt{help}) {
	die "$usage\n";
}

# Absolute path
$opt{bait}        = abs_path $opt{bait};
$opt{fish}        = abs_path $opt{fish};
$opt{output}      = abs_path $opt{output};

#Build commond line
my $cmdline = $0;
$cmdline .= " -tb $opt{tbait}" if defined $opt{tbait};
$cmdline .= " -tf $opt{tfish}" if defined $opt{tfish};
$cmdline .= " -b $opt{bait}" if defined $opt{bait};
$cmdline .= " -f $opt{fish}" if defined $opt{fish};
$cmdline .= " -o $opt{output}" if defined $opt{output};
$cmdline .= " -cb $opt{cbait}" if defined $opt{cbait};
$cmdline .= " -cf $opt{cfish}" if defined $opt{cfish};
$cmdline .= " -c $opt{contrary}" if defined $opt{contrary};

print `date`, "\n>>> fish.pl <<<\n\n";
print $cmdline,"\n\n";

print "Read the bait file.\n";

my %bait;
if ($opt{tbait} eq 'gff') {
	read_gff($opt{bait},\%bait);
}elsif ($opt{tbait} eq 'fa') {
	read_fa($opt{bait},\%bait);
}elsif ($opt{tbait} eq 'list') {
	read_list($opt{bait},\%bait);
}else {
	print "The bait file format can not be recognized.\n";
	die "$usage\n";
}

my $time=localtime;
print "The bait file has been readed: $time.\n\n";

print "Process the fish file.\n";

if ($opt{tfish} eq 'gff') {
	out_gff($opt{fish},$opt{output},\%bait);
}elsif($opt{tfish} eq 'fa'){
	out_fa($opt{fish},$opt{output},\%bait);
}elsif($opt{tfish} eq 'list'){
	out_list($opt{fish},$opt{output},\%bait);
}else {
	print "The fish file format can not be recognized.\n";
	die "$usage\n";
}

$time=localtime;
print "The fish file has been processed: $time.\n\n";

sub read_gff{
	my $file=shift;
	my $bait_hp=shift;
	my $bait_colum=($opt{cbait}) ? $opt{cbait} : 9;
	open(IN,$file)||die("fail to open $file\n");
	while (<IN>) {
		chomp;
		next if (/^\#/);
		my @temp=split (/\t/,$_);
		if ($temp[$bait_colum-1]=~/\"([^\"]+)\"/) {
			$$bait_hp{$1}=1;
		}
	}
	close(IN);
}

sub read_list{
	my $file=shift;
	my $bait_hp=shift;
	my $bait_colum=($opt{cbait}) ? $opt{cbait} : 1;
	open(IN,$file)||die("fail to open $file\n");
	while (<IN>) {
		chomp;
		next if (/^\#/);
		my @temp=split(/\t/,$_);
		$$bait_hp{$temp[$bait_colum-1]}=1;
	}
	close(IN);
}

sub read_fa{
	my $file=shift;
	my $bait_hp=shift;
	my $bait_colum=($opt{cbait}) ? $opt{cbait} : 1;
	open(IN,$file)||die("fail to open $file\n");
	$/=">";<IN>;$/="\n";
	while (<IN>) {
		my $title=$_;
		chomp $title;
		$/=">";
		my $seq=<IN>;
		chomp $seq;
		$/="\n";
		my @temp=split(/\s+/,$title);
		$$bait_hp{$temp[$bait_colum-1]}=1;
	}
	close(IN);
}


sub out_gff{
	my $file=shift;
	my $outfile=shift;
	my $fish_hp=shift;	
	my $fish_colum=($opt{cfish}) ? $opt{cfish} : 9;

	open (OUT,">$outfile")||die("fail to open $outfile.\n");
	open(IN,$file)||die("fail to open $file\n");
	while (<IN>) {
		chomp;
		next if (/^\#/);
		my @temp=split(/\t/,$_);
		if ($temp[$fish_colum-1]=~/\"([^\"]+)\"/) {
			if (defined $$fish_hp{$1} && !defined $opt{contrary}) {
				print OUT "$_\n";
			}elsif (!defined $$fish_hp{$1} && defined $opt{contrary}) {
				print OUT "$_\n";
			}
		}
	}
	close(IN);
	close OUT;
}

sub out_list{	
	my $file=shift;
	my $outfile=shift;
	my $fish_hp=shift;
	my $fish_colum=($opt{cfish}) ? $opt{cfish} : 1;

	open (OUT,">$outfile")||die("fail to open $outfile.\n");
	open(IN,$file)||die("fail to open $file\n");
	while (<IN>) {
		chomp;
		next if (/^\#/);
		my @temp=split(/\t/,$_);
		if (defined $$fish_hp{$temp[$fish_colum-1]} && !defined $opt{contrary}) {
			print OUT "$_\n";
		}elsif (!defined $$fish_hp{$temp[$fish_colum-1]} && defined $opt{contrary}) {
			print OUT "$_\n";
		}
	}
	close(IN);
	close OUT;
}

sub out_fa{
	my $file=shift;
	my $outfile=shift;
	my $fish_hp=shift;
	my $fish_colum=($opt{cfish}) ? $opt{cfish} : 1;

	open (OUT,">$outfile")||die("fail to open $outfile.\n");
	open(IN,$file)||die("fail to open $file\n");
	$/=">";<IN>;$/="\n";
	while (<IN>) {
		my $title=$_;
		chomp $title;
		$/=">";
		my $seq=<IN>;
		chomp $seq;
		$/="\n";
		my @temp=split(/\s+/,$title);
		if (defined $$fish_hp{$temp[$fish_colum-1]} && !defined $opt{contrary}) {
			print OUT ">$title\n$seq";
		}elsif (!defined $$fish_hp{$temp[$fish_colum-1]} && defined $opt{contrary}) {
			print OUT ">$title\n$seq";
		}
	}
	close(IN);
	close OUT;
}
