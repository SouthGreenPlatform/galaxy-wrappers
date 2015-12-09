#!/usr/local/bin/perl
use strict;

my $ACE = shift;
my $outfasta = shift;
my $outqual = shift;

open (OUTFASTA, ">$outfasta");
open (OUTQUAL, ">$outqual");

my $length;
my $nbreRead;

open ACE, "$ACE" or die "cannot open $ACE\n";

my $contigName;
	my $switch = 0;

while (<ACE>) {
	my $ligne = $_;
	
	if(/^\s*\n/){
		$switch =0;
	}
	
	if ($switch == 1){
		print OUTFASTA $ligne;
	}
	if ($switch == 2)
	{	
		print OUTQUAL $ligne;
	}
	
	if (/CO (.+) (\d+) (\d+) \d+ .+\n/) {
		$switch = 1;
		$contigName = $1;
		$length = $2;
		$nbreRead =$3;
		print OUTFASTA ">$contigName\n";
	}
	
	if ( /^BQ\n/ ) {
		$switch = 2;
		print OUTQUAL ">$contigName\n";
		$contigName = "";	
	}
	
	if ( /^AF/ ) {
		$switch = 0;
	}
	
		
	

	

}

close OUTFASTA;
close OUTQUAL;
close ACE;