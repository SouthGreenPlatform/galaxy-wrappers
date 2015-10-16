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
// $Id: Sensor.StartWAM.cc,v 1.11 2007-08-22 15:04:52 cnoirot Exp $
// ------------------------------------------------------------------
// File:     Sensor.StartWAM.cc
// Contents: Sensor StartWAM
// Definition of a start codon detection sensor based on a Weight Array Model
// ------------------------------------------------------------------

#include <ctype.h>
#include "Sensor.StartWAM.h"

extern Parameters PAR;

/*******************************************************
 **                  SensorStartWAM                   **
 *******************************************************/

// ----------------------
//  Default constructor.
// ----------------------
SensorStartWAM :: SensorStartWAM (int n, DNASeq *X) : Sensor(n)
{
  type = Type_Start;

  if (!IsInitialized) {
    char modelfilename[FILENAME_MAX+1];
    MarkovianOrder= PAR.getI("StartWAM.MarkovianOrder");
    NbNtBeforeATG = PAR.getI("StartWAM.NbNtBeforeATG");
    NbNtAfterATG = PAR.getI("StartWAM.NbNtAfterATG");
    MotifLength= NbNtBeforeATG + 3 + NbNtAfterATG;
    strcpy(modelfilename,PAR.getC("eugene_dir"));
    strcat(modelfilename,MODELS_DIR);
    strcat(modelfilename,"/");
    strcat(modelfilename,PAR.getC("StartWAM.modelfilename"));
    PlotScoreIncrease= 7.0;

    WAModel = new WAM(MarkovianOrder, MotifLength,"ACGT", modelfilename);
    IsInitialized = true;
  }
}

// ----------------------
//  Default destructor.
// ----------------------
SensorStartWAM :: ~SensorStartWAM ()
{
  delete WAModel;
}

// ---------------------
//  Init splice.
// ----------------------
void SensorStartWAM :: Init (DNASeq *X)
{
  ScaleCoef = PAR.getD("StartWAM.ScaleCoef*");
  ScalePenalty= PAR.getD("StartWAM.ScalePenalty*");

  if (PAR.getI("Output.graph")) Plot(X);
}

// ----------------
// Scaling
// ----------------
double SensorStartWAM :: ScaleWAMScore (double WAMScore) 
{
  return ( ScaleCoef  * WAMScore + ScalePenalty ) ;
}

// -----------------------
//  GiveInfo signal start.
// -----------------------
void SensorStartWAM :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  int i, j;
  double score;
  char* MotifExtnd = new char[MotifLength+MarkovianOrder+2]; 
  // Motif Extended = Motif + amount context 
  // (for markovian estimation of the first motif letter)
  MotifExtnd[MotifLength+MarkovianOrder+1] ='\0';

  ////////// START Forward (need enough context) //////////////
  if ( (X->IsEStart(pos,1)) &&
       (pos-NbNtBeforeATG-MarkovianOrder > 0) && 
       (pos+2+NbNtAfterATG < X->SeqLen) ) {
    score=0.0;
    j=0;
    for (i= pos-NbNtBeforeATG-MarkovianOrder; i<= pos+2+NbNtAfterATG; i++) {
      MotifExtnd[j]= toupper((*X)[i]);
      j++;
    }

    d->sig[DATA::Start].weight[Signal::Forward] += 
      ScaleWAMScore (WAModel->ScoreTheMotif(MotifExtnd));
  }

  ////////// START Reverse (need enough context)   //////////////

  if ( (X->IsEStart(pos-1,-1)) &&
       (pos-1+NbNtBeforeATG+MarkovianOrder < X->SeqLen) &&
       (pos-3-NbNtAfterATG > 0) ) {
    score=0.0;
    j=0;
    for (i= pos-1+NbNtBeforeATG+MarkovianOrder; i >= pos-3-NbNtAfterATG; i--) {
      MotifExtnd[j] = toupper((*X)(i));
      j++;
    }

    d->sig[DATA::Start].weight[Signal::Reverse] += 
      ScaleWAMScore (WAModel->ScoreTheMotif(MotifExtnd));
  }

  delete [] MotifExtnd;
}
// ----------------------------
//  Normalize Score for Plot
// ----------------------------
inline double SensorStartWAM :: NormalizePlot(double x, double n) 
{
  return Min(x,n)/n;
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorStartWAM :: Plot(DNASeq *X)
{  
  int pos;
  double plotweight;
  DATA data;

  for(pos=0; pos < X-> SeqLen ; pos++){

    data.sig[DATA::Start].weight[Signal::Forward] = 0.0 ;
    data.sig[DATA::Start].weight[Signal::Reverse] = 0.0 ;

    GiveInfo (X, pos, &data);

    if (data.sig[DATA::Start].weight[Signal::Forward] != 0) {
      plotweight= data.sig[DATA::Start].weight[Signal::Forward] + PlotScoreIncrease;
      if (data.sig[DATA::Start].weight[Signal::Forward] > 0)
      	PlotStart(pos,(pos%3)+1,NormalizePlot(plotweight,10.0));
    }
    if (data.sig[DATA::Start].weight[Signal::Reverse] != 0) {
      plotweight= data.sig[DATA::Start].weight[Signal::Reverse] + PlotScoreIncrease;
      if (data.sig[DATA::Start].weight[Signal::Reverse] > 0)
	PlotStart(pos,-((X->SeqLen-pos)%3)-1,NormalizePlot(plotweight,10.0));
    }
  }
}

// ------------------
//  Post analyse
// ------------------
void SensorStartWAM :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}
