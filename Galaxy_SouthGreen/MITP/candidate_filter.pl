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

	Usage: perl $program --hair_rnafold|-hr <*_hairpin_unique.RNAfold> --mat_fa|-mf <*_mature.fa> --mat_gff|-mg <*_mature.gff> --hair_fa|-hf <*_hairpin.fa> --hair_gff|-hg <*_hairpin.gff>
	
	The following options are necessary.
	--hair_rnafold|-hr	<*_hairpin_unique.RNAfold>	:the RNAfold result for *hairpin_unique.fa
	--mat_fa|-mf	<*_mature.fa>	:candidate mature sequence file in fasta format
	--mat_gff|-mg	<*_mature.gff>	:candidate mature gff file in gff2 format
	--hair_fa|-hf	<*_hairpin.fa>	:candidate hairpin sequence file in fasta format
	--hair_gff|-hg	<*_hairpin.gff>	:candidate hairpin gff file in gff2 format

	The following options are optional.
	--output_dir|-o		<output_dir>	:default is sample under the current directory, you should change it if you have more than one sample to distinguish different samples
	--prefix|-p		<prefix=MI>	:prefix for file and miRNA ID, if you want to compare different samples, you should use different prefix, the prefix should be end with letters, default is MI
	--alignsoft|-as	<blast|blat>	:in default, blast was used.
	--mfei		<mfei=-0.85>	:the mfei means the minimal folding free energy index (MFEI). In default, keep miRNA which mfei value is equal or lower than -0.85
	--mfe		<mfe=-25>	:the mfe means the minimal folding free energy (MFE). In default, keep miRNA which mfe value is equal or lower than -25
	--hairlen	<len=50>	:in default, keep miRNA which hairpin length is equal or larger than 50
	
	we also provide other filter parameters if you want to further filter candidate miRNA. In default, these parameters are not used

	--pairbase			:minimum pair base number
	--pairpercent			:minimum pair base percent(%)
	--amfe				:the amfe means adjusted mfe (AMFE) represented the mfe of 100 nucleotides. You can set the maximum amfe value 
	--minapercent & --maxapercent	:minimum and maximum A base percent(%)
	--minupercent & --maxapercent	:minimum and maximum U base percent(%)
	--mingpercent & --maxapercent	:minimum and maximum G base percent(%)
	--mincpercent & --maxapercent	:minimum and maximum C base percent(%)
	--minaupercent & --maxapercent	:minimum and maximum AU base percent(%)
	--mingcpercent & --maxapercent	:minimum and maximum GC base percent(%)

	we can do self alignment for candidate miRNA to filter candidate miRNA with same mature and hairpin sequence, but located in different genome region.
	--selfalign|-sa			:if you assign this parameter, it means the program will filter candidate miRNA according to self alignment for mature and hairpin sequence, as default, does not do anything.
	--selfidentity|-si		<selfidentity=100>	:minimum identity percent(%) for filter self alignment result, default is 100
	--selfoverrate|-sor		<overrate=80>	:minimum overlap rate percent(%) for read cluster, default is 80

	we can draw second structure for candidate miRNA.
	--figure	<yes|no>	:default is yes, the program will produce second structure figure in pdf format for all candidate miRNA in a subdirectory named $opt{prefix}\_figure at output directory.

	--help|-h				:print the usage information

USAGE

#Gather input
&GetOptions(
	"hair_rnafold|hr=s"		=>\$opt{hair_rnafold},
	"mat_fa|mf=s"		=>\$opt{mat_fa},
	"mat_gff|mg=s"		=>\$opt{mat_gff},
	"hair_fa|hf=s"		=>\$opt{hair_fa},
	"hair_gff|hg=s"		=>\$opt{hair_gff},
	"output_dir|o=s"	=>\$opt{output_dir},
	"prefix|p=s"		=>\$opt{prefix},
	"alignsoft|as=s"	=>\$opt{alignsoft},
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
	"help|h"	=>\$opt{help},
);

#Verify input
if (!defined $opt{hair_rnafold} or !defined $opt{mat_fa} or !defined $opt{mat_gff} or !defined $opt{hair_fa} or !defined $opt{hair_gff} or defined $opt{help}) {
	die "$usage\n";
}

#Default parameters
$opt{prefix}	= "MI" unless defined $opt{prefix};
$opt{alignsoft}	= "blast" unless defined $opt{alignsoft};
$opt{mfei}	= -0.85 unless defined $opt{mfei};
$opt{mfe}	= -25 unless defined $opt{mfe};
$opt{hairlen}	= 50 unless defined $opt{hairlen};
$opt{selfidentity}	= 100 unless defined $opt{selfidentity};
$opt{selfoverrate}	= 80 unless defined $opt{selfoverrate};
$opt{figure}	= "yes" unless defined $opt{figure};

# Absolute path
$opt{hair_rnafold}= abs_path $opt{hair_rnafold} if defined $opt{hair_rnafold};
$opt{mat_fa}      = abs_path $opt{mat_fa} if defined $opt{mat_fa};
$opt{mat_gff}     = abs_path $opt{mat_gff} if defined $opt{mat_gff};
$opt{hair_fa}     = abs_path $opt{hair_fa} if defined $opt{hair_fa};
$opt{hair_gff}    = abs_path $opt{hair_gff} if defined $opt{hair_gff};

#Output directory
$opt{output_dir} = $ENV{'PWD'}."/sample"  unless defined $opt{output_dir};
$opt{output_dir} = abs_path $opt{output_dir} if defined $opt{output_dir};
`mkdir -p $opt{output_dir}` unless -d $opt{output_dir};
chdir $opt{output_dir};

#Build commond line
my $cmdline = $program;
$cmdline .= " -hr $opt{hair_rnafold}" if defined $opt{hair_rnafold};
$cmdline .= " -mf $opt{mat_fa}" if defined $opt{mat_fa};
$cmdline .= " -mg $opt{mat_gff}" if defined $opt{mat_gff};
$cmdline .= " -hf $opt{hair_fa}" if defined $opt{hair_fa};
$cmdline .= " -hg $opt{hair_gff}" if defined $opt{hair_gff};
$cmdline .= " -o $opt{output_dir}"  if defined $opt{output_dir};
$cmdline .= " -p $opt{prefix}" if defined $opt{prefix};
$cmdline .= " -as $opt{alignsoft}" if defined $opt{alignsoft};
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

print `date`, "\n>>> $program <<<\n\n";
print $cmdline,"\n\n";

my ($time,);

#Obtain hairpin sequence and structure attribute values
print "Obtain hairpin sequence and structure attribute values.\n";
open (OUT,">$opt{prefix}\_hairpin_unique.attribute")||die("Cannot write to $opt{prefix}\_hairpin_unique.attribute file.\n");
print OUT "#miRNA\tHairpin_length\tPair_base\tPair_percent(%)\tA_content(%)\tU_content(%)\tG_content(%)\tC_content(%)\tAU_content(%)\tGC_content(%)\tMFE\tAMFE\tMFEI\n";
open (IN,"<$opt{hair_rnafold}")||die("Cannot read $opt{hair_rnafold}.\n");
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
print "The hairpin sequence and structure attribute values have been obtained: $time.\n\n";

my %id=();
#Filter candidate miRNA according to the miRNA sequence and structure attribute values
print "Filter candidate miRNA according to the miRNA sequence and structure attribute values.\n";
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

	system "perl $programdir/fish.pl -tb list -tf fa -b $opt{prefix}\_hairpin_unique.filter -f $opt{hair_fa} -o $opt{prefix}\_hairpin_unique.filter.fa";
	system "perl $programdir/fish.pl -tb list -tf fa -b $opt{prefix}\_hairpin_unique.filter -f $opt{mat_fa} -o $opt{prefix}\_mature_unique.filter.fa";
	
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
print "The candidate hairpin sequence was filtered according to the hairpin sequence and structure attribute values: $time.\n\n";

#Obtain the mature and hairpin sequence and gff file for final candidate miRNA
print "Obtain the mature and hairpin sequence and gff file for final candidate miRNA.\n";

if (-e "$opt{prefix}\_hairpin_unique.filter2") {
	system "perl $programdir/fish.pl -tb list -tf fa -b $opt{prefix}\_hairpin_unique.filter2 -f $opt{hair_fa} -o $opt{prefix}\_candidate_hairpin.fa";
	system "perl $programdir/fish.pl -tb list -tf fa -b $opt{prefix}\_hairpin_unique.filter2 -f $opt{mat_fa} -o $opt{prefix}\_candidate_mature.fa";
	system "perl $programdir/fish.pl -tb list -tf gff -b $opt{prefix}\_hairpin_unique.filter2 -f $opt{hair_gff} -o $opt{prefix}\_candidate_hairpin.gff";
	system "perl $programdir/fish.pl -tb list -tf gff -b $opt{prefix}\_hairpin_unique.filter2 -f $opt{mat_gff} -o $opt{prefix}\_candidate_mature.gff";
}else {
	system "perl $programdir/fish.pl -tb list -tf fa -b $opt{prefix}\_hairpin_unique.filter -f $opt{hair_fa} -o $opt{prefix}\_candidate_hairpin.fa";
	system "perl $programdir/fish.pl -tb list -tf fa -b $opt{prefix}\_hairpin_unique.filter -f $opt{mat_fa} -o $opt{prefix}\_candidate_mature.fa";
	system "perl $programdir/fish.pl -tb list -tf gff -b $opt{prefix}\_hairpin_unique.filter -f $opt{hair_gff} -o $opt{prefix}\_candidate_hairpin.gff";
	system "perl $programdir/fish.pl -tb list -tf gff -b $opt{prefix}\_hairpin_unique.filter -f $opt{mat_gff} -o $opt{prefix}\_candidate_mature.gff";
}
$time=localtime;
print "The mature and hairpin sequence and gff file for final candidate miRNA have been obtained: $time.\n\n";

if ($opt{figure} eq "yes") {

	#Obtain the second structure for hairpin sequence of candidate miRNA.
	print "Obtain the second structure for hairpin sequence of candidate miRNA.\n";

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
	print "The second structure for hairpin sequence of candidate miRNA has obtained: $time.\n\n";
}

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
