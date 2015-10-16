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
// $Id: Sensor.IfElse.h,v 1.8 2005-08-30 12:39:57 bardou Exp $
// ------------------------------------------------------------------
// File:     Sensor.IfElse.h
// Contents: Sensor IfElse
// ------------------------------------------------------------------

#ifndef  SENSOR_IFELSE_H_INCLUDED
#define  SENSOR_IFELSE_H_INCLUDED

#include "../../Sensor.h"

/*************************************************************
 **                     SensorIfElse                        **
 *************************************************************/
class SensorIfElse : public Sensor
{
 private:
  Sensor *sensorIf;
  Sensor *sensorElse;

 public:
  SensorIfElse  (int n, DNASeq *X);
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot(DNASeq *);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorIfElse * builder0(int n, DNASeq *X) {  return new SensorIfElse(n, X); }

#endif
