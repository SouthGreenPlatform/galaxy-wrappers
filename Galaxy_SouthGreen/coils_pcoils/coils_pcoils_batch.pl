#!/usr/bin/perl

use strict;
use Switch;
use Getopt::Long;
use Bio::SeqIO;

my $usage = qq~Usage:$0 <args> [<opts>]
where <args> are:
    -f, --fasta         <Fasta input>
    -o, --out           <output>
    -w, --window        <window length>
    -t, --threshold     <threshold: Residues below this score are predicted to be in a coiled coil. Default: 0.025> 
    -p, --path          <path of the executable>   
~;
$usage .= "\n";

my ($fasta,$out,$window_length,$path,$threshold);
$threshold = 0.025;

GetOptions(
	"fasta=s"    => \$fasta,
	"out=s"      => \$out,
	"window=s"   => \$window_length,
        "threshold=s"   => \$threshold,
        "path=s"     => \$path
);


die $usage
  if ( !$fasta ||  !$out || !$window_length || !$path);




my $in  = Bio::SeqIO->new(-file => $fasta , '-format' => 'Fasta');
open(my $OUT,">$out");
while ( my $seq = $in->next_seq() ) 
{
	my $id = $seq -> id();
	$id =~s/\|/_/g;
	my $out = Bio::SeqIO->new(-file => ">$id.fa" , '-format' => 'Fasta');
	$out->write_seq($seq);
	system("$path/PCOILS_v1.0.1/run_coils -win $window_length < $id.fa >$id.out");

	my @positions;
	my $go = 0;
	my $start;
	my @intervals;
	my $position;
	open(my $O,"$id.out");
	<$O>;<$O>;
	while(<$O>)
	{
		if (/\s+(\d+)\s+\w\s+\w\s+\d+\.*\d*\s+(\d+\.*\d*)/)
		{
			$position = $1;
			my $pscore = $2;
			if ($pscore >= $threshold && !$start)
			{
				$start = $position;
			}
			if ($pscore < $threshold && $start)
			{
				my $end = $position-1;
				push(@intervals,"$start-$end");
				$start = 0;
			}
		}
	}
	close($O);
	if ($start && !@intervals)
	{
		my $end = $position-1;
		push(@intervals,"$start-$end");
	}
	if (scalar @intervals)
	{
		print $OUT "$id";
		foreach my $interval(@intervals)
		{
			print $OUT "	" . $interval;
		}
		print $OUT "\n";
	}
	unlink("$id.fa");
	#unlink("$id.out");
}
close($OUT);
