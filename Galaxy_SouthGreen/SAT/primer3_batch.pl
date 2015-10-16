#!/usr/bin/perl

use strict;
use Getopt::Long;
use Cwd 'abs_path';
my ($current_path, $nameprg) = $0 =~ /(.+)\/(.+)/;

my $usage = qq~Usage:perl $0 <args> [<opts>]
where <args> are:
    -i, --input         <FASTA file>
    -l  --list_targets  <list of targets>
    -o, --output        <name for output file>

where <opts> are:
    -t, --tm            <minTm,maxTm,optimalTm>
    -s, --size          <minsize,maxsize,optimalsize>
    -g, --gc            <minGC,maxGC>
~;
$usage .= "\n";

my ($fasta_file,$list_of_targets,$outfile);
my $Tm = "50,63,55";
my $size = "18,27,20";
my $gc = "20,80";

my $flanking_size = 2000;

GetOptions(
	"input=s"         => \$fasta_file,
	"list_targets=s"  => \$list_of_targets,
	"output=s"        => \$outfile,
	"tm=s"            => \$Tm,
	"size=s"          => \$size,
	"gc=s"            => \$gc
);

my ($minTm,$maxTm,$optTm) = split(",",$Tm);
my ($minsize,$maxsize,$optsize) = split(",",$size);
my ($mingc,$maxgc) = split(",",$gc);

die $usage
  if ( !$fasta_file || !$list_of_targets || !$outfile);


open(OUT,">$outfile");

my %sequences;
my $id = "";
my $seq = "";
open(FASTA,$fasta_file);
while(<FASTA>)
{
	if (/>([^\s]+)/)
	{
		$id = $1;
		$seq = "";
	}
	else
	{
		my $line = $_;
		$line=~s/\n//gs;
		$line=~s/\r//gs;
		$seq.=$line;
		$sequences{$id} = $seq;
	}
}
close(FASTA);


my $nb = 0;
open(TARGET,$list_of_targets);
while(<TARGET>)
{
	$nb++;
	my @infos = split("\t",$_);
	my $id = $infos[0];
	my $sequence = $sequences{$id};
	
	my $start = $infos[5];
	my $end = $infos[6];
	my $size = $end - $start + 1;
	my $target = $start . "_" . $end;
	
	
	# get subsequence around the target to accelerate the primer search
	my $begin = 0;
	my $length = $flanking_size * 2;
	my $new_target_position = $start;
	if ($start > $flanking_size)
	{
		$begin = $start - $flanking_size;
		$new_target_position = $flanking_size;
	}
	if (length($sequence) - $begin < ($flanking_size * 2))
	{
		$length = length($sequence) - $begin;
	}
	my $subseq = substr($sequence,$begin,$length);
	
	if ($id && $target)
	{
		my $pwd = `pwd`;
		$pwd=~s/\n//g;
	
		open(PRIMER3_INPUT,">$pwd/$id.$target.primer3_input");
		print PRIMER3_INPUT "SEQUENCE_ID=$id.$target\n";
		print PRIMER3_INPUT "SEQUENCE_TEMPLATE=$subseq\n";
		print PRIMER3_INPUT "SEQUENCE_TARGET=$new_target_position,$size\n";
		print PRIMER3_INPUT "PRIMER_TASK=pick_detection_primers\n";
		print PRIMER3_INPUT "PRIMER_PICK_LEFT_PRIMER=1\n";
		print PRIMER3_INPUT "PRIMER_PICK_INTERNAL_OLIGO=1\n";
		print PRIMER3_INPUT "PRIMER_PICK_RIGHT_PRIMER=1\n";
		
		print PRIMER3_INPUT "PRIMER_OPT_SIZE=$optsize\n";
		print PRIMER3_INPUT "PRIMER_MIN_SIZE=$minsize\n";
		print PRIMER3_INPUT "PRIMER_MAX_SIZE=$maxsize\n";
		
		print PRIMER3_INPUT "PRIMER_MIN_TM=$minTm\n";
		print PRIMER3_INPUT "PRIMER_MAX_TM=$maxTm\n";
		print PRIMER3_INPUT "PRIMER_OPT_TM=$optTm\n";
		
		print PRIMER3_INPUT "PRIMER_MIN_GC=$mingc\n";
		print PRIMER3_INPUT "PRIMER_MAX_GC=$maxgc\n";
		
		print PRIMER3_INPUT "PRIMER_MAX_NS_ACCEPTED=1\n";
		print PRIMER3_INPUT "PRIMER_PRODUCT_SIZE_RANGE=150-250 100-300 301-400 401-500 501-600 601-700 701-850 851-1000\n";
		print PRIMER3_INPUT "P3_FILE_FLAG=1\n";
		print PRIMER3_INPUT "SEQUENCE_INTERNAL_EXCLUDED_REGION=$new_target_position,$size\n";
		print PRIMER3_INPUT "PRIMER_EXPLAIN_FLAG=1\n";
		print PRIMER3_INPUT "=\n";
		close(PRIMER3_INPUT);
		
		#my $command = "$current_path/primer3-2.2.3/src/primer3_core < $pwd/$id.$target.primer3_input >>$pwd/$id.$target.primer3_output";
		my $command = "primer3_core < $pwd/$id.$target.primer3_input >>$pwd/$id.$target.primer3_output";
		system($command);
		
		
		open(PRIMER_OUT,"$pwd/$id.$target.primer3_output");
		while(<PRIMER_OUT>)
		{
			if (!/SEQUENCE_TEMPLATE/)
			{
				print OUT $_;	
			}
		}
		close(PRIMER_OUT);
		
		unlink("$pwd/$id.$target.primer3_output");
		unlink("$pwd/$id.$target.primer3_input");
		unlink("$pwd/$id.rev");
		unlink("$pwd/$id.for");
		unlink("$pwd/$id.int");
	}
}
close(TARGET);
close(OUT);

