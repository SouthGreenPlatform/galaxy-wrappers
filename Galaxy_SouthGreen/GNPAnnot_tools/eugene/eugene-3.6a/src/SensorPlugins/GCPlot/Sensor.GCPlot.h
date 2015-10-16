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
// $Id: Sensor.GCPlot.h,v 1.4 2005-08-30 12:39:02 bardou Exp $
// ------------------------------------------------------------------
// File:     Sensor.GCPlot.h
// Contents: Sensor GCPlot  
// ------------------------------------------------------------------

#ifndef  SENSOR_GCPLOT_H_INCLUDED
#define  SENSOR_GCPLOT_H_INCLUDED

#include "../../Sensor.h"

/*************************************************************
 **                      SensorGCPlot
 *************************************************************/
class SensorGCPlot : public Sensor
{
 private:
  int Color;
  int Window;
  int Up;
  int Over;
  double Zoom;
  double Zoom3;

 public:
  SensorGCPlot            (int n, DNASeq *X);
  virtual ~SensorGCPlot   ();
  virtual void Init       (DNASeq *X);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorGCPlot* builder0( int n, DNASeq *X) { return new SensorGCPlot(n, X);}

#endif
