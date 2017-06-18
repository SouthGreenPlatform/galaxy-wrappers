#!/usr/bin/perl

use strict;
use Switch;
use Getopt::Long;
use Bio::SeqIO;

my $HAPLOPHYLE_EXE = "java -Xmx2048m -jar /usr/local/bioinfo/sniplay/Haplophyle/NetworkCreator_fat.jar";
my $NEATO_EXE = "neato";
my $CONVERT_EXE = "convert";
my $RSCRIPT_EXE = "/usr/local/bioinfo/R/default/bin/Rscript";

my $usage = qq~Usage:$0 <args> [<opts>]
where <args> are:
    -i, --input         <input>
<opts>
    -r, --rscript       <R scripts in zipfile>
~;
$usage .= "\n";

my ($infile,$r_zip);


GetOptions(
	"input=s"    => \$infile,
	"rscript=s"  => \$r_zip
);


die $usage
  if ( !$infile);


if (-e $r_zip)
{
	system("unzip $r_zip");
}
  
my %inputs;
my $gene;
my $data = "";
open(my $I,$infile);
while(<$I>) 
{
	if (/===(.*)===/)
	{
		$gene = $1;
	}
	else
	{
		$inputs{$gene}.= $_;
	}
}
close($I);


foreach my $gene(keys(%inputs))
{
	open(F,">$gene.distinct_haplo.fas");
	print F $inputs{$gene};
	close(F);
	
	my $input   = "$gene.distinct_haplo.fas";
    my $outfile = "$gene.network.dot";
    my $out_png = "$gene.network.png";

	my $command = "$HAPLOPHYLE_EXE -in $input -out $outfile >> $gene.haplophyle.log 2>&1";
	system($command);

	

	
	open(OUTFILE,"$outfile");
	open(OUTFILE2,">$outfile.vis");
	print OUTFILE2 "var nodes = [\n";
	open(INFILE,"$input");
	while(<INFILE>)
	{
		if (/>haplo(\d+)\|(\d+)/)
		{
			print OUTFILE2 "{id: $1, label: 'Haplo$1',shape:'image',radius:'$2'}\n";
		}
	}
	close(INFILE);
	print OUTFILE2 "];\n";
	print OUTFILE2 "var edges = [\n";
	my $n = 0;
	while(<OUTFILE>)
	{
		my $line = $_;
		if (/^haplo(\d+) -- haplo(\d+)/)
		{
			my $from = $1;
			my $to = $2;
			
			print OUTFILE2 "{from: $from, to: $to,length:'7',style: 'line',color:'black',},\n";
		}	

	}
	close(OUTFILE);
	print OUTFILE2 "];\n";
	close(OUTFILE2);
}

system("zip networks.zip *network.dot");
#system("zip pies.zip *pie.eps");

