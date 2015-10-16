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
// $Id: Sensor.SpliceWAM.cc,v 1.11 2007-08-22 15:03:53 cnoirot Exp $
// ------------------------------------------------------------------
// File:     Sensor.SpliceWAM.cc
// Contents: Sensor SpliceWAM
// A simple splice site detection sensor based on a Weight Array Model
// ------------------------------------------------------------------

#include <ctype.h>
#include "Sensor.SpliceWAM.h"

extern Parameters PAR;

//#define NORM(x,n) (((n)+(Max(-(n),x)))/(n))
#define NORM(x,n) ( (Min(x,n))/(n))
#define plotscoreincrease 10

/*************************************************************
 **                  SensorSpliceWAM                   **
 *************************************************************/

// ----------------------
//  Default constructor.
// ----------------------
SensorSpliceWAM :: SensorSpliceWAM (int n, DNASeq *X) : Sensor(n)
{
  type = Type_Acc|Type_Don;

  if (!IsInitialized) {
    char donmodelfilename[FILENAME_MAX+1];
    char accmodelfilename[FILENAME_MAX+1];
    MarkovianOrder= PAR.getI("SpliceWAM.MarkovianOrder");
    NbNtBeforeGT = PAR.getI("SpliceWAM.NbNtBeforeGT");
    NbNtAfterGT = PAR.getI("SpliceWAM.NbNtAfterGT");
    DonorSiteLength= NbNtBeforeGT + 2 + NbNtAfterGT;
    NbNtBeforeAG = PAR.getI("SpliceWAM.NbNtBeforeAG");
    NbNtAfterAG = PAR.getI("SpliceWAM.NbNtAfterAG");
    AcceptorSiteLength= NbNtBeforeAG + 2 + NbNtAfterAG;
    
    strcpy(donmodelfilename,PAR.getC("eugene_dir"));
    strcat(donmodelfilename,MODELS_DIR);
    strcat(donmodelfilename,"/");
    strcat(donmodelfilename,PAR.getC("SpliceWAM.donmodelfilename"));
    strcpy(accmodelfilename,PAR.getC("eugene_dir"));
    strcat(accmodelfilename,MODELS_DIR);
    strcat(accmodelfilename,"/");
    strcat(accmodelfilename,PAR.getC("SpliceWAM.accmodelfilename"));
    
    DonWAModel= new WAM(MarkovianOrder, DonorSiteLength,"ACGT", donmodelfilename);
    AccWAModel= new WAM(MarkovianOrder, AcceptorSiteLength, "ACGT", accmodelfilename);
    IsInitialized = true;
  }
}

// ----------------------
//  Default destructor.
// ----------------------
SensorSpliceWAM :: ~SensorSpliceWAM ()
{
  delete [] DonWAModel;
  delete [] AccWAModel;
}

// ---------------------
//  Init splice.
// ----------------------
void SensorSpliceWAM :: Init (DNASeq *X)
{
  AccScaleCoef = PAR.getD("SpliceWAM.AccScaleCoef*");
  AccScalePenalty= PAR.getD("SpliceWAM.AccScalePenalty*");
  DonScaleCoef = PAR.getD("SpliceWAM.DonScaleCoef*");
  DonScalePenalty= PAR.getD("SpliceWAM.DonScalePenalty*");

  if (PAR.getI("Output.graph")) Plot(X);
}

// -----------------------
//  GiveInfo signal splice.
// -----------------------
void SensorSpliceWAM :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  int i, j;
  double score;
  char* DonSite = new char[DonorSiteLength+MarkovianOrder+2];
  DonSite[DonorSiteLength+MarkovianOrder+1]= '\0';
  char* AccSite = new char[AcceptorSiteLength+MarkovianOrder+2];
  AccSite[AcceptorSiteLength+MarkovianOrder+1]= '\0';

  ////////// ACCEPTOR Forward (need enough context) //////////////
  if ( (((*X)[pos-2] == 'a') && ((*X)[pos-1] == 'g')) &&
       (pos-2-NbNtBeforeAG-MarkovianOrder > 0) && 
       (pos-1+NbNtAfterAG < X->SeqLen) ) {
    score=0.0;
    j=0;
    for (i= pos-2-NbNtBeforeAG-MarkovianOrder; i<= pos-1+NbNtAfterAG; i++) {
      AccSite[j] = toupper((*X)[i]);
      j++;
    }

    d->sig[DATA::Acc].weight[Signal::Forward] +=  (AccScaleCoef * AccWAModel->ScoreTheMotif(AccSite)) + AccScalePenalty;
  }

  ////////// ACCEPTOR Reverse (need enough context)   //////////////
  if ( (((*X)(pos) == 'g') && ((*X)(pos+1) == 'a')) &&
       (pos+1+NbNtBeforeAG+MarkovianOrder < X->SeqLen) &&
       (pos-NbNtAfterAG > 0) ) {
    score=0.0;
    j=0;
    for (i= pos+1+NbNtBeforeAG+MarkovianOrder; i >= pos-NbNtAfterAG; i--) {
      AccSite[j] = toupper((*X)(i));
      j++;
      }

    d->sig[DATA::Acc].weight[Signal::Reverse] += (AccScaleCoef * AccWAModel->ScoreTheMotif(AccSite)) + AccScalePenalty;
  }

  ////////// DONOR Forward (need enough context) //////////////
  if ( (((*X)[pos] == 'g') && ((*X)[pos+1] == 't')) &&
       (pos-NbNtBeforeGT-MarkovianOrder > 0) && 
       (pos+1+NbNtAfterGT < X->SeqLen) ) {
    score=0.0;
    j=0;
    for (i= pos-NbNtBeforeGT-MarkovianOrder; i<= pos+1+NbNtAfterGT; i++) {
      DonSite[j] = toupper((*X)[i]);
      j++;
    }

    d->sig[DATA::Don].weight[Signal::Forward] += (DonScaleCoef * DonWAModel->ScoreTheMotif(DonSite)) + DonScalePenalty;
  }

  ////////// DONOR Reverse (need enough context)   //////////////
  if ( (((*X)(pos-2) == 't') && ((*X)(pos-1) == 'g')) &&
       (pos-1+NbNtBeforeGT+MarkovianOrder < X->SeqLen) &&
       (pos-2-NbNtAfterGT > 0) ) {
    score=0.0;
    j=0;
    for (i= pos-1+NbNtBeforeGT+MarkovianOrder; i >= pos-2-NbNtAfterGT; i--) {
      DonSite[j] = toupper((*X)(i));
      j++;
    }
    d->sig[DATA::Don].weight[Signal::Reverse] += (DonScaleCoef * DonWAModel->ScoreTheMotif(DonSite)) + DonScalePenalty;
  }
  delete [] DonSite; delete [] AccSite;
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorSpliceWAM :: Plot(DNASeq *X)
{  
  int i,pos;
  DATA data;

  for(pos=0; pos < X-> SeqLen ; pos++){

    for (i= DATA::Acc; i <= DATA::Don; i++)
      data.sig[i].weight[Signal::Forward] = data.sig[i].weight[Signal::Reverse] = 0.0;

    GiveInfo (X, pos, &data);

    if (data.sig[DATA::Acc].weight[Signal::Forward] != 0 ) {
      data.sig[DATA::Acc].weight[Signal::Forward] += plotscoreincrease;
      if (data.sig[DATA::Acc].weight[Signal::Forward] > 0 )
	PlotAcc(pos,1,NORM(data.sig[DATA::Acc].weight[Signal::Forward],20.0));
    }
    if (data.sig[DATA::Don].weight[Signal::Forward] != 0 ) {
      data.sig[DATA::Don].weight[Signal::Forward] += plotscoreincrease;
      if (data.sig[DATA::Don].weight[Signal::Forward] > 0 )
	PlotDon(pos,1,NORM(data.sig[DATA::Don].weight[Signal::Forward],20.0));
    }
    if (data.sig[DATA::Acc].weight[Signal::Reverse] != 0 ) {
      data.sig[DATA::Acc].weight[Signal::Reverse] += plotscoreincrease;
      if (data.sig[DATA::Acc].weight[Signal::Reverse] > 0 )
	PlotAcc(pos,-1,NORM(data.sig[DATA::Acc].weight[Signal::Reverse],20.0));
    }
    if (data.sig[DATA::Don].weight[Signal::Reverse] != 0 ) {
      data.sig[DATA::Don].weight[Signal::Reverse] += plotscoreincrease;
      if (data.sig[DATA::Don].weight[Signal::Reverse] > 0 )
	PlotDon(pos,-1,NORM(data.sig[DATA::Don].weight[Signal::Reverse],20.0));
    }
/* to print all the sites: (differents colors if negatives)    
    if (data.sig[DATA::Acc].weight[Signal::Forward] != 0 )
      PlotBarF(pos,4,0.5,NORM((fabs)(data.sig[DATA::Acc].weight[Signal::Forward]),20.0), (data.sig[DATA::Acc].weight[Signal::Forward] >0) ? 4 : 7 );
    if (data.sig[DATA::Don].weight[Signal::Forward] != 0 )
      PlotBarF(pos,4,0.5,NORM((fabs)(data.sig[DATA::Don].weight[Signal::Forward]),20.0), (data.sig[DATA::Don].weight[Signal::Forward] >0) ? 5 : 6 );
    if (data.sig[DATA::Acc].weight[Signal::Reverse] != 0 )
      PlotBarF(pos,-4,0.5,NORM((fabs)(data.sig[DATA::Acc].weight[Signal::Reverse]),20.0), (data.sig[DATA::Acc].weight[Signal::Reverse] >0) ? 4 : 7 );
    if (data.sig[DATA::Don].weight[Signal::Reverse] != 0 )
      PlotBarF(pos,-4,0.5,NORM((fabs)(data.sig[DATA::Don].weight[Signal::Reverse]),20.0),(data.sig[DATA::Don].weight[Signal::Reverse] >0) ? 5 : 6 );
*/
  }
}

// ------------------
//  Post analyse
// ------------------
void SensorSpliceWAM :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}
