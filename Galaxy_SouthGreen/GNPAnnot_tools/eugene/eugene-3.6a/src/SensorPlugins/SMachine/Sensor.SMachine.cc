// VERSION BIDOUILLEE - THOMAS
// ------------------------------------------------------------------
// Copyright (C) 2004 INRA <eugene@ossau.toulouse.inra.fr>
//
// This program is open source; you can redistribute it and/or modify
// it under the terms of the Artistic License (see LICENSE file).
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
//
// You should have received a copy of Artistic License along with
// this program; if not, please see http://www.opensource.org
//
// $Id: Sensor.SMachine.cc,v 1.21 2009-01-12 14:42:43 sallet Exp $
// ------------------------------------------------------------------
// File:     Sensor.SMachine.cc
// Contents: Sensor SMachine
// ------------------------------------------------------------------

#include "Sensor.SMachine.h"
#include"../../System.h"

extern Parameters PAR;

#define NORM(x,n) Max(3.0,(x)+(n))/(n)

/*************************************************************
 **                        SensorSpliceMachine              **
 *************************************************************/

// ----------------------
//  Default constructor.
// ----------------------
SensorSMachine :: SensorSMachine (int n, DNASeq *X) : Sensor(n)
{
  char *seqname;
  char tempname[FILENAME_MAX+1];

  type = Type_Acc|Type_Don|Type_Start;

  isScaled = PAR.getI("SMachine.isScaled",GetNumber()) ;

  fprintf(stderr, "Probing SpliceMachine (splice sites)..........");  
  fflush(stderr);

  seqname = PAR.getC("fstname");
  strcpy(tempname,seqname);
  inputFormat_ = to_string(PAR.getC("SMachine.format", GetNumber(),1));

  if ( inputFormat_ == "GFF3" ) // load from GFF3 file
  {
    strcat(tempname,".spliceM.gff3");
    ReadMachineGff3(tempname, X->SeqLen);
    CheckSplices(X,vPosAccF, vPosDonF, vPosAccR, vPosDonR);
    fprintf(stderr,"  done\n");
    fprintf(stderr, "Probing SpliceMachine (starts)................");  
    fflush(stderr);
	
  }
  else // load from native format file
  {
    strcat(tempname,".spliceMAD");
    if (!ProbeFile(NULL,tempname)) SpliceMachine();
 
    ReadSMachineSplices(tempname, X->SeqLen);
    fprintf(stderr,"  done\n");
    CheckSplices(X,vPosAccF, vPosDonF, vPosAccR, vPosDonR);

    fprintf(stderr, "Probing SpliceMachine (starts)................");  
    fflush(stderr);

    strcpy(tempname,seqname);
    strcat(tempname,".spliceMSt");
    if (!ProbeFile(NULL,tempname)) SpliceMachine();
    ReadSMachineStarts(tempname, X->SeqLen);
  }
  fprintf(stderr,"  done\n");
  CheckStart(X,vPosF, vPosR);
}

// ----------------------
//  Default destructor.
// ----------------------
SensorSMachine :: ~SensorSMachine ()
{
  vPosAccF.clear();  vPosAccR.clear();
  vPosDonF.clear();  vPosDonR.clear();
  vValAccF.clear();  vValAccR.clear();
  vValDonF.clear();  vValDonR.clear();

  vPosF.clear();     vValF.clear();
  vPosR.clear();     vValR.clear();
}

// -------------------------------------
//  Scaling modes for the 2 signal edges
// -------------------------------------
inline double ScaleIt(double w, double B, double P, int Scaled)
{
  
  switch (Scaled) {
  case 0:
    return B*w-P;
    break;

  case 1:
  case 2:
    return B*log(w)-P;
    break;
    
  default:
    fprintf(stderr,"Error: incorrect value for parameter IsScaled\n");
    exit(1);
  }
}

inline double ScaleItNo(double w, double B,double P, int Scaled)
{
  switch (Scaled) {
    
  case 0:
  case 2:
    return 0.0;
    break;
    
  case 1:
    return log(1.0-pow(w,B)*exp(-P));
    break;
    
  default:
    fprintf(stderr,"Error: incorrect value for parameter IsScaled\n");
    exit(1);
  }
}
// ----------------------
//  Init SMachine.
// ----------------------
void SensorSMachine :: Init (DNASeq *X)
{
  accB = PAR.getD("SMachine.accB*",GetNumber());
  accP = PAR.getD("SMachine.accP*",GetNumber());
  donB = PAR.getD("SMachine.donB*",GetNumber());
  donP = PAR.getD("SMachine.donP*",GetNumber());

  transSpliceB = PAR.getD("SMachine.tSpliceB*",GetNumber());

  startP = PAR.getD("SMachine.startP*",GetNumber());
  startB = PAR.getD("SMachine.startB*",GetNumber());

  indexR = indexF = 0;
  iAccF = iDonF = iAccR = iDonR = 0;
  PositionGiveInfo = -1;

  if (PAR.getI("Output.graph")) Plot(X);
}
// -----------------------------
//  Read SpliceMachine Splice Site files
// -----------------------------
void SensorSMachine :: ReadSMachineSplices(char *name, int SeqLen)
{
  FILE *fp;
  char buf[FILENAME_MAX];
  char type[13];  // max len for "acceptor_rev"
  size_t len = 0;
  int i,pos,end =0;
  double force;

  if (!(fp = fopen(name, "r"))) {
    fprintf(stderr, "Cannot open splice sites file %s\n", name);
    exit(2);
  }
  // skip first line
  fgets(buf, FILENAME_MAX-1, fp);
  
  len =1;
  while (!end) {
    i = fscanf(fp,"%d %s %lf\n", &pos, type, &force);
    len ++;
    if (i == EOF) { end = 1; break;}
    if (i < 3) { end =2; break;}

    if (strcmp(type,"acceptor") == 0) {
      vPosAccF.push_back(pos);
      vValAccF.push_back(force);
    } else if (strcmp(type,"donor") == 0) {
      vPosDonF.push_back(pos-1);
      vValDonF.push_back(force);
    } else if (strcmp(type,"acceptor_rev") == 0) {
      vPosAccR.push_back(pos-1);
      vValAccR.push_back(force);
    } else if (strcmp(type,"donor_rev") == 0) {
      vPosDonR.push_back(pos);
      vValDonR.push_back(force);
    } else end = 2;
  }
  fclose(fp);

  if (end ==2) {
    fprintf(stderr, "Error in SpliceMachine splice site file %s, line %d\n", name, len);
    exit(2);
  }
}
// --------------------------
//  Read start file.
// --------------------------
void SensorSMachine :: ReadSMachineStarts(char *name, int Len)
{
  
  FILE *fp;
  char type[10];  // max len for "start_rev"
  int len = 0;
  int i,pos,end =0;
  double force;

  if (!(fp = fopen(name, "r"))) {
    fprintf(stderr, "Cannot open start file %s\n", name);
    exit(2);
  }

  len =1;
  while (!end) {
    i = fscanf(fp,"%d %s %lf\n", &pos, type, &force);
    len ++;
    if (i == EOF) {end = 1; break;}
    if (i < 3) {end =2; break;}

    if (strcmp(type,"start") == 0) {
      vPosF.push_back(pos-1);
      vValF.push_back(force);
    } else if (strcmp(type,"start_rev") == 0) {
      vPosR.push_back(pos);
      vValR.push_back(force);
    } else end = 2;
  }
  fclose(fp);
  
  if (end ==2) {
    fprintf(stderr, "Error in SpliceMachine start file %s, line %d\n", name, len);
    exit(2);
  }
}

// ------------------------
//  GiveInfo signal SMachine. 
// ------------------------
void SensorSMachine :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  bool update = false;
  if ( (PositionGiveInfo == -1) || (pos != PositionGiveInfo+1) ) update = true; // update indexes on vectors
  PositionGiveInfo = pos;

  // Accepteur Forward
  if(!vPosAccF.empty()) {
    if (update) 
      iAccF = lower_bound(vPosAccF.begin(), vPosAccF.end(), pos)-vPosAccF.begin();
    if((iAccF<(int)vPosAccF.size()) && (vPosAccF[iAccF] == pos)) {
      d->sig[DATA::Acc].weight[Signal::Forward] += ScaleIt(vValAccF[iAccF],accB,accP,isScaled);
      /*fprintf(stderr, "% 5.5i AccF\t%lf\n", pos, X->IsAcc(pos-2, 1));*/
      d->sig[DATA::Acc].weight[Signal::Forward] += log(X->IsAcc(pos-2, 1));
      d->sig[DATA::tStart].weight[Signal::Forward] += ScaleIt(vValAccF[iAccF],transSpliceB,0,0);
      d->sig[DATA::Acc].weight[Signal::ForwardNo] += ScaleItNo(vValAccF[iAccF],accB,accP,isScaled);
      iAccF++;
    }
  }
  
  // Accepteur Reverse
  if (!vPosAccR.empty()) {
    if (update)
      iAccR = lower_bound(vPosAccR.begin(), vPosAccR.end(), pos)-vPosAccR.begin();

    if((iAccR<(int)vPosAccR.size()) && (vPosAccR[iAccR] == pos)) {
      d->sig[DATA::Acc].weight[Signal::Reverse] += ScaleIt(vValAccR[iAccR], accB, accP, isScaled);
      /*fprintf(stderr, "% 5.5i AccR\t%lf\n", pos, X->IsAcc(pos+1, -1));*/
      d->sig[DATA::Acc].weight[Signal::Reverse] += log(X->IsAcc(pos+1, -1));
      d->sig[DATA::tStart].weight[Signal::Reverse] += ScaleIt(vValAccR[iAccR],transSpliceB,0,0);
      d->sig[DATA::Acc].weight[Signal::ReverseNo] += ScaleItNo(vValAccR[iAccR], accB, accP, isScaled);
      iAccR++;
    }
  }

  // Donneur Forward
  if(!vPosDonF.empty()) {
    if (update) 
      iDonF = lower_bound(vPosDonF.begin(), vPosDonF.end(), pos)-vPosDonF.begin();

    if ((iDonF<(int)vPosDonF.size()) && (vPosDonF[iDonF] == pos)) {
      d->sig[DATA::Don].weight[Signal::Forward] += ScaleIt(vValDonF[iDonF], donB, donP,isScaled);
      /*fprintf(stderr, "% 5.5i DonF\t%lf\n", pos, X->IsDon(pos, 1));*/
      d->sig[DATA::Don].weight[Signal::Forward] += log(X->IsDon(pos, 1));
      d->sig[DATA::Don].weight[Signal::ForwardNo] += ScaleItNo(vValDonF[iDonF], donB, donP,isScaled);
      iDonF++;
    }
  }
  
  // Donneur Reverse
  if(!vPosDonR.empty()) {
    if (update) 
      iDonR = lower_bound(vPosDonR.begin(), vPosDonR.end(), pos)-vPosDonR.begin();

    if((iDonR<(int)vPosDonR.size()) && (vPosDonR[iDonR] == pos)) {
      d->sig[DATA::Don].weight[Signal::Reverse] += ScaleIt(vValDonR[iDonR], donB, donP,isScaled);
      /*fprintf(stderr, "% 5.5i DonR\t%lf\n", pos, X->IsDon(pos-1, -1));*/
      d->sig[DATA::Don].weight[Signal::Reverse] += log(X->IsDon(pos-1, -1));
      d->sig[DATA::Don].weight[Signal::ReverseNo] += ScaleItNo(vValDonR[iDonR], donB, donP,isScaled);
      iDonR++;
    }
  }

  // Start Forward
  if(!vPosF.empty()) {
    if (update) 
      indexF = lower_bound(vPosF.begin(), vPosF.end(), pos)-vPosF.begin();
    
    if((indexF<(int)vPosF.size()) && (vPosF[indexF] == pos)) {
      d->sig[DATA::Start].weight[Signal::Forward] += ScaleIt(vValF[indexF], startB,startP,isScaled);
      /*fprintf(stderr, "% 5.5i StartF\t%lf\n", pos, X->IsStart(pos, 1));*/
      d->sig[DATA::Start].weight[Signal::Forward] += log(X->IsEStart(pos, 1));
      d->sig[DATA::Start].weight[Signal::ForwardNo] += ScaleItNo(vValF[indexF], startB,startP,isScaled);
      indexF++;
    }
  }
  
  // Start Reverse
  if (!vPosR.empty()) {
    if (update) 
      indexR = lower_bound(vPosR.begin(), vPosR.end(), pos)-vPosR.begin();

    if((indexR<(int)vPosR.size()) && (vPosR[indexR] == pos)) {
      d->sig[DATA::Start].weight[Signal::Reverse] += ScaleIt(vValR[indexR], startB,startP,isScaled);
	/*for (int z=-3;z<3;z++) {*/
	      /*fprintf(stderr, "% 5.5i %+i StartR\t%lf\n", pos, z, X->IsStart(pos+z, -1));*/
	/*}*/
      /*fprintf(stderr, "% 5.5i StartR\t%lf\n", pos, X->IsStart(pos-1, -1));*/
      d->sig[DATA::Start].weight[Signal::Reverse] += log(X->IsEStart(pos-1, -1));
      d->sig[DATA::Start].weight[Signal::ReverseNo] += ScaleItNo(vValR[indexR], startB,startP,isScaled);
      indexR++;
    }
  }
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorSMachine :: Plot(DNASeq *X)
{
  double f;

  for (int i =0; i < (int)vPosAccF.size(); i++) {
	  f = ScaleIt(vValAccF[i], accB, accP, isScaled);
    PlotAcc(vPosAccF[i], 1, NORM(log(f),20.0));
  }
  
  for (int i =0; i < (int)vPosDonF.size(); i++) {
    f = ScaleIt(vValDonF[i], donB, donP, isScaled);
    PlotDon(vPosDonF[i], 1, NORM(log(f),20.0));
  }
  
  for (int i =0; i < (int)vPosAccR.size(); i++) {
    f = ScaleIt(vValAccR[i], accB, accP, isScaled);
    PlotAcc(vPosAccR[i], -1, NORM(log(f),20.0));
  }

  for (int i =0; i < (int)vPosDonR.size(); i++) {
    f = ScaleIt(vValDonR[i], donB, donP, isScaled);
    PlotDon(vPosDonR[i], -1, NORM(log(f),20.0));
  }

  for (int i =0; i < (int)vPosF.size(); i++) {
    f = ScaleIt(vValF[i], startB,startP,isScaled);
    PlotStart(vPosF[i], (vPosF[i]%3)+1, NORM(log(f),10.0));
  }
  
  for (int i =0; i < (int)vPosR.size(); i++) {
    f = ScaleIt(vValR[i], startB,startP,isScaled);
    PlotStart(vPosR[i], -((X->SeqLen-vPosR[i])%3)-1, NORM(log(f),10.0));
  }
}

// ------------------
//  Post analyse
// ------------------
void SensorSMachine :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}
// -----------------------
//  SpliceMachine launcher
// -----------------------
void SensorSMachine :: SpliceMachine()
{
  char *SMcmd=PAR.getC("SMachine.cmd");
  char s1[FILENAME_MAX+1];

  strcpy(s1,SMcmd);
  strcat(s1, " ");
  strcat(s1, PAR.getC("fstname"));

  system(s1);

  return;
}


// --------------------------
//  Read SMachine GFF3 file 
// --------------------------

void SensorSMachine :: ReadMachineGff3(char name[FILENAME_MAX+1], int SeqLen)
{
  
  GeneFeatureSet * geneFeatureSet = new GeneFeatureSet (name);

  vector< GeneFeature *>::iterator it = geneFeatureSet->getIterator();
  int nbElement=geneFeatureSet->getNbFeature();

  int i=0;
  while ( i<nbElement )
  {
    GeneFeature * tmpFeature = *it;
    string idSo=tmpFeature->getType();
    
    //Get SO code if feature correspond to the name or the synonym.
    if ( idSo.find("SO:") == string::npos )
    {
      string tmp=GeneFeatureSet::soTerms_->getIdFromName(idSo);
      idSo=tmp;
    }
     // Forward
    if ( tmpFeature->getLocus()->getStrand() == '+' ) 
    {
      if (idSo == "SO:0000163") //donor
      {
	vPosDonF.push_back( tmpFeature->getLocus()->getStart()-1 );
	vValDonF.push_back( tmpFeature->getScore() );
      }
      if (idSo == "SO:0000164")  //acceptor
      {
	vPosAccF.push_back(tmpFeature->getLocus()->getEnd() );
	vValAccF.push_back( tmpFeature->getScore() );
      }
      if (idSo == "SO:0000318")  //start
      {
	vPosF.push_back(tmpFeature->getLocus()->getStart() -1 );
	vValF.push_back( tmpFeature->getScore() );
      }

    }
    
    // Reverse
    if ( tmpFeature->getLocus()->getStrand() == '-' ) {
      if(idSo == "SO:0000163") {
	vPosDonR.push_back( tmpFeature->getLocus()->getEnd() );
	vValDonR.push_back( tmpFeature->getScore() );
      }
      if (idSo == "SO:0000164")  //acceptor
      {
	vPosAccR.push_back( tmpFeature->getLocus()->getStart()-1 );
	vValAccR.push_back( tmpFeature->getScore() );
      }
      if (idSo == "SO:0000318")  //start
      {
	vPosR.push_back( tmpFeature->getLocus()->getEnd() );
	vValR.push_back( tmpFeature->getScore() );
      }
    }
    it++;
    i++;
  }
  delete geneFeatureSet;

}

