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
// $Id: Sensor.PatConst.cc,v 1.8 2005-08-30 12:40:01 bardou Exp $
// ------------------------------------------------------------------
// File:     Sensor.PatConst.cc
// Contents: Sensor PatConst (Constant Pattern)
// ------------------------------------------------------------------

#include <ctype.h>

#include "Sensor.PatConst.h"

extern Parameters PAR;

/*************************************************************
 **                       SensorPatConst
 *************************************************************/

// ----------------------
// Default constructor.
// ----------------------
SensorPatConst :: SensorPatConst (int n) : Sensor(n)
{
  newStatePos = PAR.getI("PatConst.newStatePos", GetNumber());
  patType     = PAR.getC("PatConst.type",        GetNumber());
  pattern     = PAR.getC("PatConst.pat",         GetNumber());
  patLen      = strlen(pattern);
  for(int i=0; i<patLen; i++)
    pattern[i] = tolower(pattern[i]);
  
  // Type initialisation

  if(!strcmp(patType, "start")) {
    sigTypeIndex = DATA::Start;   type = Type_Start;  }
  else if(!strcmp(patType, "insertion")) {
    sigTypeIndex = DATA::Ins;     type = Type_FS; }
  else if(!strcmp(patType, "deletion"))  {
    sigTypeIndex = DATA::Del;     type = Type_FS; }
  else if(!strcmp(patType, "transstart")) {
    sigTypeIndex = DATA::tStart;  type = Type_TStart; }
  else if(!strcmp(patType, "transstop")) {
    sigTypeIndex = DATA::tStop;   type = Type_TStop; }
  else if(!strcmp(patType, "donor")) {
    sigTypeIndex = DATA::Don;     type = Type_Don; }
  else if(!strcmp(patType, "acceptor")) {
    sigTypeIndex = DATA::Acc;     type = Type_Acc; }
  else if(!strcmp(patType, "stop")) {
    sigTypeIndex = DATA::Stop;     type = Type_Stop;}
  else {
    fprintf(stderr, "Error \"%s\" undefined type in the parameter file"
	    " for PatConst.pat[%d]\n", patType, GetNumber());
    fprintf(stderr, "Type must be : start, insertion, deletion, "
	    "transstart, transstop, stop,\n               "
	    "acceptor or donor.\n");
    exit(2);
  }
}

// ----------------------
//  Default destructor.
// ----------------------
SensorPatConst :: ~SensorPatConst ()
{
}
// ----------------------
//  Init.
// ----------------------
void SensorPatConst :: Init (DNASeq *X)
{
  patP   = PAR.getD("PatConst.patP*",  GetNumber());
  patPNo = PAR.getD("PatConst.patPNo*",GetNumber());

  if (PAR.getI("Output.graph")) Plot(X);
}

// ------------------------
//  GiveInfo.
// ------------------------
void SensorPatConst :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  int  i;
  bool isSigF = true;
  bool isSigR = true;
 
  for(i=0; i<patLen; i++)
    if((*X)[pos+i-newStatePos+1] != pattern[i]) {
      isSigF=false;
      break;
    }
  for(i=patLen-1; i>-1; i--)
    if((*X)(pos-i+newStatePos-2) != pattern[i]) {
      isSigR=false;
      break;
    }
  
  if(isSigF) {
    d->sig[sigTypeIndex].weight[Signal::Forward]   += patP;
    d->sig[sigTypeIndex].weight[Signal::ForwardNo] += patPNo;
  }
  if(isSigR) {
    d->sig[sigTypeIndex].weight[Signal::Reverse]   += patP;
    d->sig[sigTypeIndex].weight[Signal::ReverseNo] += patPNo;
  }
}

// ----------------------------
//  Plot Sensor information.
// ----------------------------
void SensorPatConst :: Plot(DNASeq *X)
{
  int  i,l;
  bool isSigF = true;
  bool isSigR = true;
 
  // Shame on us.Boyer-Moore or KMP needed

  for (l=0; l <= X->SeqLen;l++) {
    isSigF = true;
    isSigR = true;
    
    for (i=0; i<patLen; i++)
      if ((*X)[l+i-newStatePos+1] != pattern[i]) {
	isSigF=false;
	break;
      }
    
    for (i=patLen-1; i>-1; i--)
      if ((*X)(l-i+newStatePos-2) != pattern[i]) {
	isSigR=false;
	break;
      }
    
    if (isSigF) {
      if (sigTypeIndex == DATA::Start)
	PlotStart(l,(l%3)+1,0.5);
      else if (sigTypeIndex == DATA::Don)
	PlotDon(l,1,0.5);
      else if(sigTypeIndex == DATA::Acc)
	PlotAcc(l,1,0.5);
      else if(sigTypeIndex == DATA::Stop)
	PlotStop(l,(l%3)+1,1);
    }

    if(isSigR) {
      if (sigTypeIndex == DATA::Start)
	PlotStart(l,-((X->SeqLen-l)%3)-1,0.5);
      else if (sigTypeIndex == DATA::Don)
	PlotDon(l,-1,0.5);
      else if(sigTypeIndex == DATA::Acc)
	PlotAcc(l,-1,0.5);
      else if(sigTypeIndex == DATA::Stop)
	PlotStop(l,-((X->SeqLen-l)%3)-1,1);
    }
  }
}

// ------------------
//  Post analyse.
// ------------------
void SensorPatConst :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}
