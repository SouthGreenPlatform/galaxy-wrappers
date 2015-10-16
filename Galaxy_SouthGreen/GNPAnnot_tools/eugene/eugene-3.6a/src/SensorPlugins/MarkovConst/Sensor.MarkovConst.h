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
// $Id: Sensor.MarkovConst.h,v 1.9 2005-08-30 12:39:57 bardou Exp $
// ------------------------------------------------------------------
// File:     Sensor.MarkovConst.h
// Contents: Sensor MarkovConst
// ------------------------------------------------------------------

#ifndef  SENSOR_MARKOVCONST_INCLUDED
#define  SENSOR_MARKOVCONST_INCLUDED

#include "../../Sensor.h"

/*************************************************************
 **                     SensorMarkovConst                    **
 *************************************************************/
class SensorMarkovConst : public Sensor
{
 private:
  double transCodant, transIntron, transIntronU, transInter, transUTR5, transUTR3;
  double minGC,maxGC;

 public:
  SensorMarkovConst  (int n, DNASeq *X);
  virtual ~SensorMarkovConst   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorMarkovConst * builder0(int n, DNASeq *X) {  return new SensorMarkovConst(n, X); }

#endif
