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
// $Id: Sensor.NStart.cc,v 1.27 2009-01-12 14:36:23 sallet Exp $
// ------------------------------------------------------------------
// File:     Sensor.NStart.cc
// Contents: Sensor NStart
// ------------------------------------------------------------------

#include "Sensor.NStart.h"

extern Parameters PAR;

#define NORM(x,n) (((n)+(Max(-(n),x)))/(n))

/*************************************************************
 **                       SensorNStart                      **
 *************************************************************/

// ----------------------
// Default constructor.
// ----------------------
SensorNStart :: SensorNStart (int n, DNASeq *X) : Sensor(n)
{
  char tempname[FILENAME_MAX+1];

  type = Type_Start;
  
  fprintf(stderr, "Reading start file (NetStart).................");
  fflush(stderr);
  
  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".starts");

  inputFormat_ = to_string(PAR.getC("NStart.format", GetNumber(),1));
  if ( inputFormat_ == "GFF3" )
  {
    strcat(tempname,".gff3");
    ReadNStartGff3(tempname, X->SeqLen);
    fprintf(stderr,"forward, reverse done\n");
    fflush(stderr);
  }
  else
  {
    ReadNStartF(tempname, X->SeqLen);
    fprintf(stderr,"forward,");
    fflush(stderr);
  
    strcpy(tempname,PAR.getC("fstname"));
    strcat(tempname,".startsR");
    ReadNStartR(tempname, X->SeqLen);
    fprintf(stderr," reverse done\n");
    fflush(stderr);
  }
  CheckStart(X,vPosF, vPosR);

  // vectors for reverse are put in the increasing order
  reverse(vPosR.begin(), vPosR.end()); 
  reverse(vValR.begin(), vValR.end()); 
}

// ----------------------
//  Default destructor.
// ----------------------
SensorNStart :: ~SensorNStart ()
{
  vPosF.clear();
  vValF.clear();
  vPosR.clear();
  vValR.clear();
}

// ----------------------
//  Init start.
// ----------------------
void SensorNStart :: Init (DNASeq *X)
{
  startP = PAR.getD("NStart.startP*",GetNumber());
  startB = PAR.getD("NStart.startB*",GetNumber());

  indexR = indexF = 0;
  PositionGiveInfo = -1;
  
  if (PAR.getI("Output.graph")) Plot(X);
}

// --------------------------
//  Read start forward file.
// --------------------------
void SensorNStart :: ReadNStartF (char name[FILENAME_MAX+1], int Len)
{
  FILE *fp;
  int i,j = -1;
  double force;

  if (!(fp = fopen(name, "r"))) {
    fprintf(stderr, "cannot open start file %s\n", name);
    exit(2);
  }
  
  while (1) {
    i = fscanf(fp,"%d %lf %*s\n", &j, &force);
    if (i == EOF) break;
    if (i < 2) {
      fprintf(stderr, "Error in start file %s, position %d\n", name, j);
      exit(2);
    }
    vPosF.push_back( j-1 );
    vValF.push_back( force );
  }
  if (j == -1) fprintf(stderr,"WARNING: empty NetStart file !\n");
  fclose(fp);
}

// --------------------------
//  Read start reverse file.
// --------------------------
void SensorNStart :: ReadNStartR (char name[FILENAME_MAX+1], int Len)
{
  FILE *fp;
  int i,j = -1;
  double force;

  if (!(fp = fopen(name, "r"))) {
    fprintf(stderr, "cannot open start file %s\n", name);
    exit(2);
  }

  while (1) {
    i = fscanf(fp,"%d %lf %*s\n", &j, &force);
    if (i == EOF) break;
    if (i < 2) {
      fprintf(stderr, "Error in start file %s, position %d\n", name, j);
      exit(2);
    }
    j = Len-j+2;
    vPosR.push_back( j-1 );
    vValR.push_back( force );
  }
  if (j == -1) fprintf(stderr,"WARNING: empty NetStart file !\n");
  fclose(fp);
}


void SensorNStart :: ReadNStartGff3 (char name[FILENAME_MAX+1], int Len)
{
  
  GeneFeatureSet * geneFeatureSet = new GeneFeatureSet (name);
  vector< GeneFeature *>::iterator it = geneFeatureSet->getIterator();
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
    if ( idSo != "SO:0000318")
    {
	fprintf(stderr,"WARNING: NetStart plugin doesn't accept feature %s !\n",tmpFeature->getType().c_str());
 	i++;
    	it++;
	continue;
    }
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



// ------------------------
//  GiveInfo signal start.
// ------------------------
void SensorNStart :: GiveInfo (DNASeq *X, int pos, DATA *d)
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
void SensorNStart :: Plot(DNASeq *X)
{
  double f;

  for (int i =0; i < (int)vPosF.size(); i++) {
    f = pow(vValF[i], startB)*(exp(-startP));
    PlotStart(vPosF[i], (vPosF[i]%3)+1, NORM(log(f),4.0));
  }

  for (int i =0; i < (int)vPosR.size(); i++) {
    f = pow(vValR[i], startB)*(exp(-startP));
    PlotStart(vPosR[i], -((X->SeqLen-vPosR[i])%3)-1, NORM(log(f),4.0));
  }
}

// ------------------
//  Post analyse
// ------------------
void SensorNStart :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}
