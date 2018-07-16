#!/usr/bin/perl

use strict;
use Getopt::Long;

use Cwd;
my $dir = getcwd;


# -m $modules -p $probes -i $pdf -s $minsize -n $threshold 
my $usage = qq~Usage:$0 <args> [<opts>]
where <args> are:
    -c, --counts        <counts>
    -m, --modules       <modules>
    -f, --figures_pdf   <figures pdf>
    -p, --path_wgcna    <path for wgcna rscript>
<opts> are:
    -s, --size_min      <minimum size of cluster. Default: 30>
    -o, --overlap       <Threshold for topological overlap for connecting genes in the network. Default: 0.1>
~;
$usage .= "\n";

my ($counts,$traits,$modules,$pdf,$minsize,$threshold,$path_wgcna);
$minsize=30;
$threshold=0.1;

GetOptions(
        "counts=s"     => \$counts,
        "modules=s"    => \$modules,
        "figures_pdf=s"=> \$pdf,
        "size_min=s"   => \$minsize,
        "overlap=s"    => \$threshold,
	"path_wgcna=s" => \$path_wgcna
);


die $usage
  if ( !$counts || !$modules || !$pdf || !$path_wgcna);

system("perl $path_wgcna/createInputsWGCNA.pl -i $counts -c $counts.counts -t $counts.traits"); 
system("Rscript $path_wgcna/wgcna.R -c $counts.counts -t $counts.traits -m $modules.modules -p $modules.probes -i $pdf -s $minsize -n $threshold");

my %hash_genes;
open(P,"$modules.probes");
while(<P>){
	my $line = $_;
	$line =~s/\n//g;$line =~s/\r//g;
	$line =~s/\"//g;
	my ($num,$gene) = split(/\t/,$line);
	$hash_genes{$num} = $gene;
}
close(P);

my %hash_cluster;
open(M,"$modules.modules");
while(<M>){
	my $line = $_;
        $line =~s/\n//g;$line =~s/\r//g;
        $line =~s/\"//g;
        my ($num,$cluster) = split(/\t/,$line);
	my $gene = $hash_genes{$num};
	$hash_cluster{$cluster}.=$gene.",";
}
close(M);

open(MOD,">$modules");
foreach my $cl(keys(%hash_cluster)){
	my @genes = split(",",$hash_cluster{$cl});
	foreach my $gene(@genes){
		print MOD "$gene	$cl\n";
	}
}
close(MOD);
