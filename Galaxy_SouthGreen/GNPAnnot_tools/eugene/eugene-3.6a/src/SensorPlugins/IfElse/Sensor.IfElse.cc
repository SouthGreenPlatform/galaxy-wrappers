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
// $Id: Sensor.IfElse.cc,v 1.18 2005-08-30 12:39:57 bardou Exp $
// ------------------------------------------------------------------
// File:     Sensor.IfElse.cc
// Contents: Sensor IfElse
// ------------------------------------------------------------------

#include<iostream>

#include "Sensor.IfElse.h"

#include "../../MSensor.h"

extern Parameters PAR;
extern MasterSensor* MS;

/*************************************************************
 **                      SensorIfElse                       **
 *************************************************************/
// ----------------------
//  Default constructor.
// ----------------------
SensorIfElse :: SensorIfElse (int n, DNASeq *X) : Sensor(n)
{
  std::string name;

  // create instance of sensors
  name = ((std::string) "Sensor.") + PAR.getC("IfElse.SensorIf",n);
  sensorIf   = MS->MakeSensor( name.c_str(), n+1, X);
  name = ((std::string) "Sensor.") + PAR.getC("IfElse.SensorElse",n);
  sensorElse = MS->MakeSensor( name.c_str(), n+1, X);

  type = sensorIf->type | sensorElse->type;
}


// ----------------------
//  Init from new seq.
// ----------------------
void SensorIfElse :: Init (DNASeq *X)
{
  sensorIf->Init(X);
  sensorElse->Init(X);
}


// ----------------------
// ----------------------
void MergeInfo(DATA *ifData, DATA *elseData, DATA *d)
{
  int i,j;

  // on fusionne les signaux dans ifData
  for (j=0; j < DATA::LastSigType; j++) {
    for (i=0; i<= Signal::Reverse;  i++)
      if (!ifData->sig[j].IsSet(i)) {
	ifData->sig[j].weight[i] = elseData->sig[j].weight[i];
	ifData->sig[j].weight[i+2] = elseData->sig[j].weight[i+2];
      }
  }

  // on repercute dans d
  for (j=0; j < DATA::LastSigType; j++) {
    for (i=0; i<= Signal::Reverse;  i++)
      if (ifData->sig[j].IsSet(i)) {
	d->sig[j].weight[i] += ifData->sig[j].weight[i];
	d->sig[j].weight[i+2] += ifData->sig[j].weight[i+2];
      }
  }

  for(i=0; i< DATA::LastContentsType; i++)
    d->contents[i] +=
      (ifData->contents[i] ? ifData->contents[i] : elseData->contents[i]);
}


// -----------------------
//  GiveInfo signal stop.
// -----------------------
void SensorIfElse :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  DATA ifData;
  DATA elseData;
  int i;

  for(i=0; i<  DATA::LastSigType;  i++)
    ifData.sig[i].Clear();
  for(i=0; i< DATA::LastContentsType; i++) ifData.contents[i] = 0.0;

  sensorIf->GiveInfo(X,pos,&ifData);

  for(i=0; i<  DATA::LastSigType;  i++)
    elseData.sig[i].Clear();
  for(i=0; i< DATA::LastContentsType; i++) elseData.contents[i] = 0.0;


  sensorElse->GiveInfo(X,pos,&elseData);
  MergeInfo(&ifData,&elseData,d);
}


// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorIfElse :: Plot(DNASeq *X)
{
  sensorIf->Plot(X);
  sensorElse->Plot(X);
}


// ------------------
//  Post analyse
// ------------------
void SensorIfElse :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
  sensorIf->PostAnalyse(pred, MINFO);
  sensorElse->PostAnalyse(pred, MINFO);
}
