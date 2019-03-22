#!/usr/local/bin/perl -w
#This progrem transfer the BLAST result include the Query sites and Sbject sites to the excel file.
$vision="Vision: 1.0";
use Getopt::Std;
getopts "o:i:e:a:l:s:";

print "*************\n*$vision*\n*************\n";
if ((!defined $opt_i)|| (!defined $opt_o) ) {
	die "************************************************************************
	Usage:transferall.pl -i filename -o outfile.
	  -h : help and usage.
	  -v : vision.
	  -e : expect value(default 10)
	  -a : identity% (default 0)
	  -l : alignment length (default 0)
	  -s : score (default 0)
************************************************************************\n";
}
$Expect = (defined $opt_e) ? $opt_e : 10;
$Length = (defined $opt_l) ? $opt_l : 0;
$Identity = (defined $opt_a) ? $opt_a : 0;
$Score = (defined $opt_s) ? $opt_s : 0;
open(Ofile,">$opt_o");print "Running.....\n";
open (F,$opt_i) || die"can't open $opt_i:$!\n";
$i=$q=$s=$a=0;
print Ofile "Query name\tLetter\tQueryX\tQueryY\tSbjctX\tSbjctY\tLength\tScore\tE value\tOverlap/total\tIdentity\tSbject Name\tAnontation\n";
while (<F>)
{
	if (/Query= (\S+)/)
	{
		if($i==1)
		{
			#print Ofile "$query^$letter^Query:$qbegin----$qend^Sbject:$sbegin----$send^$name^$annotation^$length^$score^$expect^$identity^$over\n";
			if($score >$Score && $expect<=$Expect && $over>=$Identity && $identity2 >= $Length)
			{
				$ovalap_total = "$identity1/$identity2";
				print Ofile "$query\t$letter\t$qbegin\t$qend\t$sbegin\t$send\t$length\t$score\t$expect\t$ovalap_total\t$over\t$name\t$annotation\n";
			}
			$i=$q=$s=0;
			$query=$letter=$qbegin=$qend=$sbegin=$send=$name=$annotation=$length=$score=$expect=$identity=$over=0;
		}
		$query=$1;
	}
	elsif (/\((\S+)\s+letters\)/)
	{
		$letter=$1;$letter=~s/,//;
	}
	elsif (/^>(\S*)\s+(.*)/)
	{
		if($i==1)
		{
			if($score >$Score && $expect<=$Expect && $over>=$Identity && $identity2 >= $Length)
			{
				$ovalap_total = "$identity1/$identity2";
				print Ofile "$query\t$letter\t$qbegin\t$qend\t$sbegin\t$send\t$length\t$score\t$expect\t$ovalap_total\t$over\t$name\t$annotation\n";
			}
			$i=$q=$s=0;
			$qbegin=$qend=$sbegin=$send=$name=$annotation=$length=$score=$expect=$identity=$over=0;
		}
		$name=$1;
		$annotation=$2;
		$a=1;
	}
	elsif (/Length = (\d+)/)
	{
		$length=$1;
		$a=0;
	}
	elsif ($a==1)
	{
		chomp;
		$annotation.=$_;
		$annotation=~s/\s+/ /g;
	} #This sentence could get the very long annotation that is longer than one line;
	elsif (/Score = (.+) bits.+Expect =\s+(\S+)\s*/ or /Score = (.+) bits.+Expect\(\d+\) =\s+(\S+)\s*/)
	{
		if($i==1)
		{
			if($score >$Score && $expect<=$Expect && $over>=$Identity && $identity2 >= $Length)
			{
				$ovalap_total = "$identity1/$identity2";
				print Ofile "$query\t$letter\t$qbegin\t$qend\t$sbegin\t$send\t$length\t$score\t$expect\t$ovalap_total\t$over\t$name\t$annotation\n";
			}
			$i=$q=$s=0;
			$qbegin=$qend=$sbegin=$send=$score=$expect=$identity=$over=0;
		}
	$score=$1;$expect=$2;$expect=~s/^e/1e/;
	}
	elsif (/Identities = (\d+)\/(\d+)\s+\((.{0,4})%\)/)
	{
		$identity1=$1;
		$identity2=$2;
		$over=$3;
	}
	elsif (/Query\:\s+(\d+)\s+\S+\s+(\d+)/)
	{
		if($q==0)
		{
			$qbegin=$1;
		}
		$qend=$2;
		$q=1;
	}
	elsif (/Sbjct\:\s+(\d+)\s+\S+\s+(\d+)/)
	{
		if($s==0)
		{
			$sbegin=$1;
		}
		$send=$2;
		$s=$i=1;
	}
}
if(($score >$Score && $expect<=$Expect && $over>=$Identity && $identity2 >= $Length)&&($i==1))
{
	$ovalap_total = "$identity1/$identity2";
	print Ofile "$query\t$letter\t$qbegin\t$qend\t$sbegin\t$send\t$length\t$score\t$expect\t$ovalap_total\t$over\t$name\t$annotation\n";
}
close(F);
close(Ofile);

print "Done!\n";
