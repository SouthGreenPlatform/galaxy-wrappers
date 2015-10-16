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
// $Id: Sensor.Plotter.cc,v 1.5 2005-08-30 12:40:02 bardou Exp $
// ------------------------------------------------------------------
// File:     Sensor.Plotter.cc
// Contents: Sensor Plotter
// ------------------------------------------------------------------

#include "Sensor.Plotter.h"

extern Parameters PAR;

/*************************************************************
 **                       SensorPlotter
 *************************************************************/

// ----------------------
// Default constructor.
// ----------------------
SensorPlotter :: SensorPlotter (int n, DNASeq *X) : Sensor(n)
{
  type = Type_None;

  // Save parameters to limit the map access number
  window  = PAR.getI("Output.window");
  plotGC  = PAR.getI("Plotter.GC");
  plotGC3 = PAR.getI("Plotter.GC3");
  plotAT  = PAR.getI("Plotter.A|T/A+T");
}

// ----------------------
//  Default destructor.
// ----------------------
SensorPlotter :: ~SensorPlotter ()
{
}

// ----------------------
//  Init.
// ----------------------
void SensorPlotter :: Init (DNASeq *X)
{
  if (PAR.getI("Output.graph")) Plot(X);
}

// ------------------------
//  GiveInfo.
// ------------------------
void SensorPlotter :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
}

// ----------------------------
//  Plot Sensor information.
// ----------------------------
void SensorPlotter :: Plot(DNASeq *X)
{
  double GC[3], pGC[3], A_AT, pA_AT, T_AT, pT_AT;
  unsigned int Width, pWidth;
  unsigned int GC1, AT1, A, T;
  int i;
  
  GC[0] = GC[1] = GC[2] = 0;
  A_AT = 0;
  T_AT = 0;
  
  for (i = 0; i < window/2; i++) {
    GC1 = (((*X)(i,0) & CodeG) != 0) + (((*X)(i,0) & CodeC) != 0);
    AT1 = (((*X)(i,0) & CodeA) != 0) + (((*X)(i,0) & CodeT) != 0);
    GC[i%3] += ((double)GC1/(GC1+AT1));
    
    A = (((*X)(i,0) & CodeA) != 0);
    T = (((*X)(i,0) & CodeT) != 0);
    if (AT1 != 0) {
      A_AT += ((double)A/AT1);
      T_AT += ((double)T/AT1);
    }
  }
  
  Width = window/2;
  
  for (i = 0; i < X->SeqLen; i++) {
    pGC[0] = GC[0];
    pGC[1] = GC[1];
    pGC[2] = GC[2];
    pA_AT  = A_AT;
    pT_AT  = T_AT;

    pWidth = Width;
    
    if (i-window/2 >= 0) {
      GC1 = (((*X)(i-window/2,0) & CodeG) != 0) +
	(((*X)(i-window/2,0) & CodeC) != 0);
      AT1 = (((*X)(i-window/2,0) & CodeA) != 0) +
	(((*X)(i-window/2,0) & CodeT) != 0);
      GC[(i-window/2)%3] -= ((double)GC1/(GC1+AT1));
      Width--;

      A = (((*X)(i-window/2,0) & CodeA) != 0);
      T = (((*X)(i-window/2,0) & CodeT) != 0);
      if (AT1 != 0) {
	A_AT -= ((double)A/AT1);
	T_AT -= ((double)T/AT1);
      }
    }

    if (i+window/2 < X->SeqLen) {
      GC1 = (((*X)(i+window/2,0) & CodeG) != 0) +
	(((*X)(i+window/2,0) & CodeC) != 0);
      AT1 = (((*X)(i+window/2,0) & CodeA) != 0) +
	(((*X)(i+window/2,0) & CodeT) != 0);
      GC[(i+window/2)%3] += ((double)GC1/(GC1+AT1));
      Width++;
      
      A = (((*X)(i+window/2,0) & CodeA) != 0);      
      T = (((*X)(i+window/2,0) & CodeT) != 0);
      if (AT1 != 0) {
	A_AT += ((double)A/AT1);
	T_AT += ((double)T/AT1);
      }
    }

    if (plotGC)
      PlotLine(((i == 0) ? 0 : i-1),i, 0, 0,
               ((double)(pGC[0]+pGC[1]+pGC[2])/ pWidth)*1.5 - 0.25,
               ((double)(GC[0] +GC[1] +GC[2]) / Width) *1.5 - 0.25, 5);
    
    if (plotGC3) {
      PlotLine(((i == 0) ? 0 : i-1),i, 1, 1,
               (((double)pGC[2]*3) / pWidth),
               (((double)GC[2] *3) / Width), 5);
      
      PlotLine(((i == 0) ? 0 : i-1),i, 2, 2,
               (((double)pGC[0]*3) / pWidth),
               (((double)GC[0] *3) / Width), 5);
      
      PlotLine(((i == 0) ? 0 : i-1),i, 3, 3,
               (((double)pGC[1]*3) / pWidth),
               (((double)GC[1] *3) / Width), 5);
      
      PlotLine(((i == 0) ? 0 : i-1),i, -1, -1,
               (((double)pGC[X->SeqLen%3]*3) / pWidth),
               (((double)GC[X->SeqLen%3] *3) / Width), 5);
      
      PlotLine(((i == 0) ? 0 : i-1),i, -2, -2,
               (((double)pGC[(X->SeqLen-1)%3]*3) / pWidth),
               (((double)GC[(X->SeqLen-1)%3] *3) / Width), 5);
      
      PlotLine(((i == 0) ? 0 : i-1),i, -3, -3,
               (((double)pGC[(X->SeqLen-2)%3]*3) / pWidth),
               (((double)GC[(X->SeqLen-2)%3] *3) / Width), 5);
    }

    if (plotAT) {
        PlotLine(((i == 0) ? 0 : i-1),i, +4, +4,
	       ((double)(pT_AT) / (pA_AT+pT_AT)),
	       ((double)(T_AT)  / (A_AT+T_AT)), 9);
	PlotLine(((i == 0) ? 0 : i-1),i, -4, -4,
               ((double)(pA_AT) / (pA_AT+pT_AT)),
	       ((double)(A_AT)  / (A_AT+T_AT)), 9);
    }
  }
}

// ------------------
//  Post analyse.
// ------------------
void SensorPlotter :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}
