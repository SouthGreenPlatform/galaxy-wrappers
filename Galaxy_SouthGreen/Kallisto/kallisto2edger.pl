#!/usr/bin/perl

use strict;

use strict;
use Getopt::Long;


my $usage = qq~Usage:$0 <args> [<opts>]
where <args> are:
    -f, --first         <first input files for condition1>
    -s, --second        <second input files for condition2>
    -n, --name1         <name for condition1>
    -c, --condition2    <name for condition2>
    -o, --output        <output file>
~;
$usage .= "\n";

my ($first,$second,$out,$name1,$name2);
my $name1 = "cond1";
my $name2 = "cond2";

GetOptions(
        "first=s"     => \$first,
	"second=s"    => \$second,
        "out=s"       => \$out,
        "name1=s"     => \$name1,
        "condition2=s"=> \$name2
);

die $usage
  if ( !$first || !$second || !$out || !$name1 || !$name2);



my $num = 0;
my %hash;
my @files1 = split(",",$first);
my @genes;
my @files_cond1;
foreach my $file1(@files1){
	$num++;
	open(F,$file1);
	push(@files_cond1,$file1);
	<F>;
	while(<F>){
		my @i = split("\t",$_);
		my $gene = $i[0];
		my $count = $i[3];
		my $rounded = sprintf("%.0f", $count);
		if ($num==1){push(@genes,$gene);}
		$hash{$gene}{"cond1"}{$file1} = $rounded;	
	}
	close(F);
}
my @files2 = split(",",$second);
my @files_cond2;
foreach my $file2(@files2){
	$num++;
        open(F,$file2);
        push(@files_cond2,$file2);
        <F>;
        while(<F>){
                my @i = split("\t",$_);
                my $gene = $i[0];
                my $count = $i[3];
                my $rounded = sprintf("%.0f", $count);
                $hash{$gene}{"cond2"}{$file2} = $rounded;
        }
        close(F);
}
open(O,">$out");
print O "#";
for (my $k=1;$k<= scalar @files_cond1;$k++){
	print O "	G1:$name1";
}
for (my $k=1;$k<= scalar @files_cond2;$k++){
        print O "	G2:$name2";
}
print O "\n";
print O "#Feature";
for (my $k=1;$k<= scalar @files_cond1;$k++){
        print O "	$name1.$k";
}
for (my $k=1;$k<= scalar @files_cond2;$k++){
        print O "	$name2.$k";
}
#for (my $k=1;$k<=$num;$k++){
#	print O "	$k";
#}
print O "\n";
foreach my $gene(@genes){
	print O "$gene";
	# condition1
	foreach my $file(@files_cond1){
		my $count = $hash{$gene}{"cond1"}{$file};
		print O "\t$count";
	}
	# condition2
	foreach my $file(@files_cond2){
                my $count = $hash{$gene}{"cond2"}{$file};
                print O "\t$count";
        }	
	print O "\n";
}
close(O);


