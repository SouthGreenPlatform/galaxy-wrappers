#!/usr/bin/perl

# ------------------------------------------------------------------
# Copyright (C) 2004 INRA <eugene@ossau.toulouse.inra.fr>
#
# This program is open source; you can redistribute it and/or modify
# it under the terms of the Artistic License (see LICENSE file).
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#
# You should have received a copy of Artistic License along with
# this program; if not, please see http://www.opensource.org
#
# $Id: sensor_plugin_builder.pl,v 1.4 2004-09-16 11:55:07 cros Exp $
# ------------------------------------------------------------------
# File:     sensor_plugin_builder.pl
# Contents: Just a little script to build new plugin :  
#             - create two files Sensor."name".h and Sensor."name".cc
#             - modify the makefile  
# ------------------------------------------------------------------


use strict;
use warnings;
use IO::Handle;


# input
my $sensorName = "";
#my $nbParam    = "";
my $makefile;

# output
my $sensorH;
my $sensorCC;
my $newMakefile;

($#ARGV == 0) || die ("Usage: $0 <Sensors plugins Makefile>\n");

print STDOUT "\n=================== Build new plugin ====================\n";
print STDOUT "| Input :\n";
print STDOUT "    - New sensor name (Eg: Est, EuStop)...";
$sensorName = <STDIN>;
$sensorName = ucfirst($sensorName);
chomp $sensorName;
while($sensorName eq "") {
    print STDOUT "\t=> Warning: Invalid name..........";
    $sensorName = <STDIN>;
    $sensorName = ucfirst($sensorName);
    chomp $sensorName;
}
$sensorH  = "Sensor.".$sensorName.".h";
$sensorCC = "Sensor.".$sensorName.".cc";

#print STDOUT "    - How many parameters do you need.....";
#$nbParam = <STDIN>;
#chomp $nbParam;
#while($nbParam =~ /[^0-9]/ || $nbParam eq "") {
#    print STDOUT "\t=> Warning: Integer required......";
#    $nbParam = <STDIN>;
#    chomp $nbParam;
#}

print STDOUT "---------------------------------------------------------\n";
print STDOUT "| Output :\n";


###################################################
#                Sensor."name".h                  #
###################################################
print STDOUT "    - $sensorH..........";
open  (SEN_H,">$sensorH") || die "ERROR: Could not open $sensorH";
print SEN_H  "#ifndef  SENSOR_".uc($sensorName)."_H_INCLUDED
#define  SENSOR_".uc($sensorName)."_H_INCLUDED\n
#include \"../../EuGene/Sensor.h\"\n
/*************************************************************
 **                      Sensor".$sensorName."
 *************************************************************/
class Sensor".$sensorName." : public Sensor
{
  private:
    // Methods and Data structures:
    //   Eg: std::vector<int>  ...;
    //       void ReadFile (char[FILENAME_MAX+1], int);

  public:
    Sensor".$sensorName."          (int);
    virtual ~Sensor".$sensorName." ();
    virtual void Init       (DNASeq *);
    virtual void GiveInfo   (DNASeq *X, int, DATA *);
    virtual void Plot       (DNASeq *X);
    virtual void PostAnalyse(Prediction *);
};
\nextern \"C\" Sensor".$sensorName."* builder0( int n ) { return new Sensor".$sensorName."(n);}
\n#endif";
close SEN_H;
print STDOUT "created\n";


###################################################
#                Sensor."name".cc                 #
###################################################
print STDOUT   "    - $sensorCC.........";
open  (SEN_CC, ">$sensorCC") || die "ERROR: Could not open $sensorCC";
print SEN_CC   "#include \"Sensor.".$sensorName.".h\"\n
/*************************************************************
 **                       Sensor".$sensorName."
 *************************************************************/

extern Parameters PAR;

// ----------------------
// Default constructor.
// ----------------------
Sensor".$sensorName." :: Sensor".$sensorName." (int n) : Sensor(n)
{
  // Save parameters to limit the map access number
  // ... = PAR.getD(\"".$sensorName."....parameter name...\");
}

// ----------------------
//  Default destructor.
// ----------------------
Sensor".$sensorName." :: ~Sensor".$sensorName." ()
{
  // Clear the data structures
}

// ----------------------
//  Init.
// ----------------------
void Sensor".$sensorName." :: Init (DNASeq *X)
{
  char tempname[FILENAME_MAX+1];

  // Type initialisation
  // Eg: type = Type_Start/Stop/Splice/Content/Multiple/Unknown;

  // Clear the data structures 

  //fprintf(stderr, \"Reading file.................\");
  //fflush(stderr);
  //strcpy(tempname,PAR.getC(\"fstname\"));
  //strcat(tempname,\"....extension name...\");
  //ReadFile(tempname, X->SeqLen);

  //if (PAR.getI(\"Output.graph\")) Plot(X);
}

// ------------------------
//  GiveInfo.
// ------------------------
void Sensor".$sensorName." :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  // Sensor signal : x <=> tStart/tStop/Start/Stop/Acc/Don
  //  d->sig[DATA::...x...].weight[Signal::Forward]   += log(...);
  //  d->sig[DATA::...x...].weight[Signal::ForwardNo] += log(...);
  //  d->sig[DATA::...x...].weight[Signal::Reverse]   += log(...);
  //  d->sig[DATA::...x...].weight[Signal::ReverseNo] += log(...);

  // Sensor content :
  //  for(int i=0; i<6; i++)
  //    d->ContentScore[i] += log(...); //Exons
  //  d->ContentScore[6] += log(...);   //IntronF
  //  d->ContentScore[7] += log(...);   //IntronR
  //  d->ContentScore[8] += log(...);   //InterG
  //  d->ContentScore[9] += log(...);   //UTR5'F
  //  d->ContentScore[10]+= log(...);   //UTR5'R
  //  d->ContentScore[11]+= log(...);   //UTR3'F
  //  d->ContentScore[12]+= log(...);   //UTR3'R
}

// ----------------------------
//  Plot Sensor information.
// ----------------------------
void Sensor".$sensorName." :: Plot(DNASeq *X)
{
}

// ------------------
//  Post analyse.
// ------------------
void Sensor".$sensorName." :: PostAnalyse(Prediction *pred)
{
}
";
close SEN_CC;
print STDOUT   "created\n";


###################################################
#               Makefile modification             #
###################################################
$newMakefile = "Makefile_"."$sensorName";
print STDOUT "    - $newMakefile..........";
flush STDOUT;
open(MAKE,"$ARGV[0]") || die "\n  ERROR: Could not open the sensors plugins Makefile $ARGV[0]\n";
open  (NEW_MAKE,">$newMakefile") || die "ERROR: Could not open file $newMakefile";
while (my $line=<MAKE>) {
  chomp $line;
  my $tmp = "$sensorName/Sensor.$sensorName";

  if ($line =~ /all : /) {
   print NEW_MAKE $line." $tmp.o\n"; }
  else {
   if ($line =~ /clean :/) {
    print NEW_MAKE "$tmp.o : $tmp.cc $tmp.h ../EuGene/Sensor.h\n";
    print NEW_MAKE "\t\$(CXX) \$(CFLAGS) -c $tmp.cc -o \$@\n";
    print NEW_MAKE "\t\$(CXX) \$(CFLAGS) \$(DLL)  \$@ -o $tmp\$(DLLEXT)\n";
    print NEW_MAKE "\tmv $tmp\$(DLLEXT) \$(DIR_PLUGINS)\n\n";
    print NEW_MAKE $line."\n"; }
   else { print NEW_MAKE $line."\n"; }
  }
}
print STDOUT "created\n";
print STDOUT "---------------------------------------------------------\n";
print STDOUT "| After the new sensor implementation, don't forget  to |\n";
print STDOUT "| modify the parameter file (parameters, use, prioty).  |\n";
print STDOUT "=========================================================\n\n";
