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
// $Id: Sensor.EuStop.h,v 1.13 2005-08-30 12:38:04 bardou Exp $
// ------------------------------------------------------------------
// File:     Sensor.EuStop.h
// Contents: Sensor EuStop
// ------------------------------------------------------------------

#ifndef  SENSOR_EUSTOP_H_INCLUDED
#define  SENSOR_EUSTOP_H_INCLUDED

#include "../../Sensor.h"

/*************************************************************
 **                     SensorEuStop                        **
 *************************************************************/
class SensorEuStop : public Sensor
{
 private:
  double   stopP;

 public:
  SensorEuStop            (int n, DNASeq *X);
  virtual ~SensorEuStop   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorEuStop * builder0(int n, DNASeq *X) {  return new SensorEuStop(n, X); }

#endif
