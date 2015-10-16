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
// $Id: Sensor.EuStop.cc,v 1.15 2005-08-30 12:38:04 bardou Exp $
// ------------------------------------------------------------------
// File:     Sensor.EuStop.cc
// Contents: Sensor EuStop
// ------------------------------------------------------------------

#include "Sensor.EuStop.h"

extern Parameters PAR;

/*************************************************************
 **                      SensorEuStop                       **
 *************************************************************/
// ----------------------
//  Default constructor.
// ----------------------
SensorEuStop :: SensorEuStop (int n, DNASeq *X) : Sensor(n)
{
  type = Type_Stop;
}

// ----------------------
//  Default destructor.
// ----------------------
SensorEuStop :: ~SensorEuStop ()
{
}

// ----------------------
//  Init stop.
// ----------------------
void SensorEuStop :: Init (DNASeq *X)
{
  stopP = PAR.getD("EuStop.stopP*");
  
  if (PAR.getI("Output.graph")) Plot(X);
}

// -----------------------
//  GiveInfo signal stop.
// -----------------------
void SensorEuStop :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  double theStop;
  
  if ((theStop = X->IsStop(pos-3,1)) != 0.0) {
    d->sig[DATA::Stop].weight[Signal::Forward] += -stopP+log(theStop);
    d->sig[DATA::Stop].weight[Signal::ForwardNo] += log(1.0-theStop); 
  }
      
  if ((theStop = X->IsStop(pos+2,-1)) != 0.0) {
    d->sig[DATA::Stop].weight[Signal::Reverse] += -stopP+log(theStop);
    d->sig[DATA::Stop].weight[Signal::ReverseNo] += log(1.0-theStop); 
  }
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorEuStop :: Plot(DNASeq *X)
{
  double Strength;

  for(int i = 0;  i <= X->SeqLen;  i++) {
    if ((Strength= X->IsStop(i-3,1)) != 0.0) 
      PlotStop(i,(i%3)+1,(Strength == 1.0));

    if((Strength = X->IsStop(i+2,-1)) != 0.0) 
      PlotStop(i,-((X->SeqLen-i)%3)-1,(Strength == 1.0));
  }
}

// ------------------
//  Post analyse
// ------------------
void SensorEuStop :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}
