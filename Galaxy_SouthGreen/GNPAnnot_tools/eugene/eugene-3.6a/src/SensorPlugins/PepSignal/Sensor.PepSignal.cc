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
// $Id: Sensor.PepSignal.cc,v 1.11 2009-01-12 14:38:27 sallet Exp $
// ------------------------------------------------------------------
// File:     Sensor.PepSignal.cc
// Contents: Sensor Peptide Signal
// ------------------------------------------------------------------

#include "Sensor.PepSignal.h"
#define NORM(x,n) (((n)+(Max(-(n),x)))/(n))

/*************************************************************
 **                       SensorPepSignal                   **
 *************************************************************/

extern Parameters PAR;

// ----------------------
// Default constructor.
// ----------------------
SensorPepSignal :: SensorPepSignal (int n, DNASeq *X) : Sensor(n)
{
  char *seqname;
  char tempname[FILENAME_MAX+1];
  seqname = PAR.getC("fstname");
  strcpy(tempname,seqname);
  strcat(tempname,".psignal");
  fprintf(stderr, "Probing PepSignal (starts)....");  
  fflush(stderr);

  inputFormat_ = to_string(PAR.getC("PepSignal.format", GetNumber(),1));

  if ( inputFormat_ == "GFF3" )
  {
    strcat(tempname,".gff3");
    ReadPepSignalGff3(tempname, X->SeqLen);
  }
  else
  {
     ReadPepSignalStarts(tempname, X->SeqLen);
  }
  fprintf(stderr,"done\n");
  fflush(stderr);
  CheckStart(X,vPosF, vPosR);
}

void SensorPepSignal :: ReadPepSignalGff3 (char name[FILENAME_MAX+1], int Len)
{
  
  GeneFeatureSet * geneFeatureSet = new GeneFeatureSet (name);
  vector <GeneFeature *>::iterator it = geneFeatureSet->getIterator();
  int nbElement=geneFeatureSet->getNbFeature();
  //geneFeatureSet->printFeature();
  int i=0;
  while ( i<nbElement )
  {
    //(*it)->second();
    GeneFeature * tmpFeature = *it;
    string idSo=tmpFeature->getType();
    if ( idSo.find("SO:") == string::npos )
    {
      string tmp=GeneFeatureSet::soTerms_->getIdFromName(idSo);
      idSo=tmp;
    }
     // Forward
    
    if ( tmpFeature->getLocus()->getStrand() == '+' ) 
    {
      vPosF.push_back( tmpFeature->getLocus()->getStart() -1);
      vValF.push_back( tmpFeature->getScore() );
    }
    if ( tmpFeature->getLocus()->getStrand() == '-' ) 
    {
      vPosR.push_back( tmpFeature->getLocus()->getStart()  );
      vValR.push_back( tmpFeature->getScore() );
    }
    i++;
    it++;
  }
  delete geneFeatureSet;
}

// --------------------------
//  Read start file.
// --------------------------
void SensorPepSignal :: ReadPepSignalStarts(char *name, int Len)
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
    i = fscanf(fp,"%d %s %lf %*s\n", &pos, type, &force);
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
    fprintf(stderr, "Error in  PepSignal start file %s, line %d\n", name, len);
    exit(2);
  }
}




// ----------------------
//  Default destructor.
// ----------------------
SensorPepSignal :: ~SensorPepSignal ()
{
  vPosF.clear();
  vValF.clear();
  vPosR.clear();
  vValR.clear();
}

// ----------------------
//  Init start.
// ----------------------
void SensorPepSignal :: Init (DNASeq *X)
{
  startP = PAR.getD("PepSignal.startP*",GetNumber());
  startB = PAR.getD("PepSignal.startB*",GetNumber());

  indexR = indexF = 0;
  PositionGiveInfo = -1;
  
  if (PAR.getI("Output.graph")) Plot(X);
}

// ------------------------
//  GiveInfo signal start.
// ------------------------
void SensorPepSignal :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  bool update = false;
  double f;

  if ( (PositionGiveInfo == -1) || (pos != PositionGiveInfo+1) ) update = true; // update indexes on vectors
  PositionGiveInfo = pos;
  
  // Start Forward
  if(!vPosF.empty()) {
    if (update) 
      indexF = lower_bound(vPosF.begin(), vPosF.end(), pos)-vPosF.begin();
    
    if((indexF<(int)vPosF.size()) && (vPosF[indexF] == pos)) {
      f = pow(vValF[indexF], startB)*(exp(-startP));
      d->sig[DATA::Start].weight[Signal::Forward] += log(f);
      d->sig[DATA::Start].weight[Signal::ForwardNo] += log(1.0-f);
      indexF++;
    }
  }
  
  // Start Reverse
  if (!vPosR.empty()) {
    if (update) 
      indexR = lower_bound(vPosR.begin(), vPosR.end(), pos)-vPosR.begin();
    
    if((indexR<(int)vPosR.size()) && (vPosR[indexR] == pos)) {
      f = pow(vValR[indexR], startB)*(exp(-startP));
      d->sig[DATA::Start].weight[Signal::Reverse] += log(f);
      d->sig[DATA::Start].weight[Signal::ReverseNo] += log(1.0-f);
      indexR++;
    }
  }


}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorPepSignal :: Plot(DNASeq *X)
{
  double f;

  for (int i =0; i < (int)vPosF.size(); i++) {
    PlotBarF(vPosF[i], (vPosF[i]%3)+1, 0.9, 0.2, 2);
  }

  for (int i =0; i < (int)vPosR.size(); i++) {
    PlotBarF(vPosR[i], -((X->SeqLen-vPosR[i])%3)-1, 0.9, 0.2, 2);
  }
}

// ------------------
//  Post analyse
// ------------------
void SensorPepSignal :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}


