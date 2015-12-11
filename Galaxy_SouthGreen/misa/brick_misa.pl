#!usr/local/bin/perl

use strict;
use File::Copy;


my $file=shift;

my $result_file=shift;
my $stat_file=shift;

#print "\tStart Misa\n";
system ("/home/galaxydev/galaxy/tools/SouthGreen/misa/misa.pl $file $stat_file ");

#print "\tReport Misa\n";
system ("/home/galaxydev/galaxy/tools/SouthGreen/misa/p3_in.pl $file.misa ");


#print "\tStart Primer3\n";
system ("primer3_core < $file.p3in > $file.p3out ");

#print "\tReport Primer3\n";
system ("/home/galaxydev/galaxy/tools/SouthGreen/misa/p3_out.pl $file.p3out $file.misa $result_file ");

