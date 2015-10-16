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
// $Id: Sensor.FrameShift.h,v 1.7 2005-08-30 12:38:39 bardou Exp $
// ------------------------------------------------------------------
// File:     Sensor.FrameShift.h
// Contents: Sensor FrameShift 
// ------------------------------------------------------------------

#ifndef  SENSOR_FRAMESHIFT_H_INCLUDED
#define  SENSOR_FRAMESHIFT_H_INCLUDED

#include "../../Sensor.h"

/*************************************************************
 **                     SensorFrameShift                    **
 *************************************************************/
class SensorFrameShift : public Sensor
{
 private:
  // Probability of deletion/insertion
  double insProb;
  double delProb;

 public:
  SensorFrameShift          (int n, DNASeq *X);
  virtual ~SensorFrameShift ();
  virtual void Init         (DNASeq *X);
  virtual void GiveInfo     (DNASeq *X, int, DATA *);
  virtual void Plot         (DNASeq *X);
  virtual void PostAnalyse  (Prediction *, FILE *);
};

extern "C" SensorFrameShift * builder0(int n, DNASeq *X) {  return new SensorFrameShift(n, X); }

#endif
