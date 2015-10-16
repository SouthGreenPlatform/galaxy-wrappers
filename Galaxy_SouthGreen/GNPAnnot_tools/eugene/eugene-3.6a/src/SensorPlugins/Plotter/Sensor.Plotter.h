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
// $Id: Sensor.Plotter.h,v 1.4 2005-08-30 12:40:02 bardou Exp $
// ------------------------------------------------------------------
// File:     Sensor.Plotter.h
// Contents: Sensor Plotter
// ------------------------------------------------------------------

#ifndef  SENSOR_PLOTTER_H_INCLUDED
#define  SENSOR_PLOTTER_H_INCLUDED

#include "../../Sensor.h"

/*************************************************************
 **                      SensorPlotter
 *************************************************************/
class SensorPlotter : public Sensor
{
 private:
  int window;
  int plotGC, plotGC3, plotAT;
  
 public:
  SensorPlotter          (int n, DNASeq *X);
  virtual ~SensorPlotter ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorPlotter* builder0( int n, DNASeq *X) { return new SensorPlotter(n, X);}

#endif
