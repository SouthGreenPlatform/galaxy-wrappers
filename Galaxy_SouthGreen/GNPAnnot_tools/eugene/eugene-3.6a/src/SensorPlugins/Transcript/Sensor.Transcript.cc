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
// $Id: Sensor.Transcript.cc,v 1.6 2005-08-30 12:40:07 bardou Exp $
// ------------------------------------------------------------------
// File:     Sensor.Transcript.cc
// Contents: Sensor Transcript
// ------------------------------------------------------------------

#include "Sensor.Transcript.h"

extern Parameters PAR;

#define NORM(x,n) (((n)+(Max(-(n),x)))/(n))

/*************************************************************
 **                      SensorTranscript                   **
 *************************************************************/

// ----------------------
//  Default constructor.
// ----------------------
SensorTranscript :: SensorTranscript (int n, DNASeq *X) : Sensor(n)
{
  type = Type_Start;
}

// ----------------------
//  Default destructor.
// ----------------------
SensorTranscript :: ~SensorTranscript ()
{
}

// ----------------------
//  Init start.
// ----------------------
void SensorTranscript :: Init (DNASeq *X)
{
  transStart = PAR.getD("Transcript.Start*");
  transStop = PAR.getD("Transcript.Stop*");

  if (PAR.getI("Output.graph")) Plot(X);
}

// -----------------------
//  GiveInfo signal start.
// -----------------------
void SensorTranscript :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  for (int i = Signal::Forward; i <= Signal::Reverse; i++) {
    d->sig[DATA::tStart].weight[i] -= transStart;
    d->sig[DATA::tStop].weight[i] -= transStop;
  }
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorTranscript :: Plot(DNASeq *X)
{
}

// ------------------
//  Post analyse
// ------------------
void SensorTranscript :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}
