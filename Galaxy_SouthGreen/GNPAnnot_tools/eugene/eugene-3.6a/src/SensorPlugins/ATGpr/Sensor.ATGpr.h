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
// $Id: Sensor.ATGpr.h,v 1.5 2005-08-30 12:34:19 bardou Exp $
// ------------------------------------------------------------------
// File:     Sensor.ATGpr.h
// Contents: Sensor ATGpr 
// ------------------------------------------------------------------

#ifndef  SENSOR_ATGPR_H_INCLUDED
#define  SENSOR_ATGPR_H_INCLUDED

#include "../../Sensor.h"

/*************************************************************
 **                     SensorATGpr                         **
 *************************************************************/
class SensorATGpr : public Sensor
{
 private:
  int PositionGiveInfo;
  std::vector<int>    vPosF, vPosR;
  std::vector<double> vValF, vValR;
  int indexF, indexR;
  double startP, startB;
  
  void ReadATGprF (char[FILENAME_MAX+1], int);
  void ReadATGprR (char[FILENAME_MAX+1], int);

 public:
  SensorATGpr             (int n, DNASeq *X);
  virtual ~SensorATGpr    (void);
  virtual void Init       (DNASeq *X);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorATGpr* builder0( int n, DNASeq *X) { return new SensorATGpr(n, X);}

#endif
