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
// $Id: Sensor.GCPlot.cc,v 1.4 2005-08-30 12:39:02 bardou Exp $
// ------------------------------------------------------------------
// File:     Sensor.GCPlot.cc
// Contents: Sensor GCPlot  
// ------------------------------------------------------------------

#include "Sensor.GCPlot.h"

extern Parameters PAR;

/*************************************************************
 **                       SensorGCPlot                         **
 *************************************************************/

// ----------------------
// Default constructor.
// ----------------------
SensorGCPlot :: SensorGCPlot (int n, DNASeq *X) : Sensor(n)
{
  Color = PAR.getI("GCPlot.Color");
  Window = PAR.getI("GCPlot.Smooth");
  Zoom = PAR.getD("GCPlot.Zoom");
  Zoom3 = PAR.getD("GCPlot.Zoom3");
}

// ----------------------
//  Default destructor.
// ----------------------
SensorGCPlot :: ~SensorGCPlot ()
{

}

// ----------------------
//  Init.
// ----------------------
void SensorGCPlot :: Init (DNASeq *X)
{
  char *tmp = PAR.getC("GCPlot.Up");
  int len = strlen(tmp);

  Up = 0;
  for (int i =0; i<len;i++)
    Up |= X->Nuc2Code(tmp[i]);
  
  tmp = PAR.getC("GCPlot.Over");
  len = strlen(tmp);
  Over = 0;
  for (int i =0; i<len;i++)
    Over |= X->Nuc2Code(tmp[i]);

  if(PAR.getI("Output.graph"))  Plot(X);
}

// ------------------------
//  GiveInfo.
// ------------------------
void SensorGCPlot :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
}

// ----------------------------
//  Plot Sensor information.
// ----------------------------
void SensorGCPlot :: Plot(DNASeq *X)
{
  int up[3],pup[3];
  int over[3],pover[3];
  int i;
  
  up[0] = up[1] = up[2] = 0;
  over[0] = over[1] = over[2] = 0;

  for (i = 0; i < Window/2; i++) {
    up[i%3] += (((*X)(i,0) & Up) != 0);
    over[i%3] += (((*X)(i,0) & Over) != 0);
  }

  for (i = 0; i < X->SeqLen; i++) {
    pup[0] = up[0];    pup[1] = up[1];    pup[2] = up[2];
    pover[0] = over[0];    pover[1] = over[1];    pover[2] = over[2];
  
    if (i-Window/2 >= 0) {
      up[(i-Window/2)%3] -= (((*X)(i-Window/2,0) & Up) != 0);
      over[(i-Window/2)%3] -= (((*X)(i-Window/2,0) & Over) != 0);
    }
    
    if (i+Window/2 < X->SeqLen) {
      up[(i+Window/2)%3] += (((*X)(i+Window/2,0) & Up) != 0);
      over[(i+Window/2)%3] += (((*X)(i+Window/2,0) & Over) != 0);
    }
    
    PlotLine(((i == 0) ? 0 : i-1),i, 0, 0,
	     0.5+Zoom*(((double)pup[0]+pup[1]+pup[2])/(pover[0]+pover[1]+pover[2])-0.5),
	     0.5+Zoom*(((double)up[0]+up[1]+up[2])/(over[0]+over[1]+over[2])-0.5),Color);
    
    PlotLine(((i == 0) ? 0 : i-1),i, 1, 1,
	     0.5+Zoom3*(((double)pup[2]/pover[2])-0.5),
	     0.5+Zoom3*(((double)up[2]/over[2])-0.5),Color);
    
    PlotLine(((i == 0) ? 0 : i-1),i, 2, 2,
	     0.5+Zoom3*(((double)pup[0]/pover[0])-0.5),
	     0.5+Zoom3*(((double)up[0]/over[0])-0.5),Color);
	     
    PlotLine(((i == 0) ? 0 : i-1),i, 3, 3,
	     0.5+Zoom3*(((double)pup[1]/pover[1])-0.5),
	     0.5+Zoom3*(((double)up[1]/over[1])-0.5),Color);
    
    PlotLine(((i == 0) ? 0 : i-1),i, -1, -1,
	     0.5+Zoom3*(((double)pup[X->SeqLen%3]/pover[X->SeqLen%3]) - 0.5),
	     0.5+Zoom3*(((double)up[X->SeqLen%3]/over[X->SeqLen%3]) - 0.5),Color);
      
    PlotLine(((i == 0) ? 0 : i-1),i, -2, -2,
	     0.5+Zoom3*(((double)pup[(X->SeqLen-1)%3]/pover[(X->SeqLen-1)%3])-0.5),
	     0.5+Zoom3*(((double)up[(X->SeqLen-1)%3]/over[(X->SeqLen-1)%3]) - 0.5),Color);
    
    PlotLine(((i == 0) ? 0 : i-1),i, -3, -3,
	     0.5+Zoom3*(((double)pup[(X->SeqLen-2)%3]/pover[(X->SeqLen-2)%3]) - 0.5),
	     0.5+Zoom3*(((double)up[(X->SeqLen-2)%3]/over[(X->SeqLen-2)%3]) -0.5),Color);
     }
}

// ------------------
//  Post analyse.
// ------------------
void SensorGCPlot :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}
