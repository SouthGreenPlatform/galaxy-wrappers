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
// $Id: Sensor.PatConst.h,v 1.3 2005-08-30 12:40:01 bardou Exp $
// ------------------------------------------------------------------
// File:     Sensor.PatConst.h
// Contents: Sensor PatConst (Constant Pattern)
// ------------------------------------------------------------------

#ifndef  SENSOR_PATCONST_H_INCLUDED
#define  SENSOR_PATCONST_H_INCLUDED

#include "../../Sensor.h"

/*************************************************************
 **                      SensorPatConst
 *************************************************************/
class SensorPatConst : public Sensor
{
 private:
  double patP;
  double patPNo;
  char*  pattern;
  char*  patType;
  int    newStatePos;
  int    sigTypeIndex;
  int    patLen;
  
 public:
  SensorPatConst          (int);
  virtual ~SensorPatConst ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorPatConst* builder0( int n ) { return new SensorPatConst(n);}

#endif
