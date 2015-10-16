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
// $Id: Sensor.NStart.h,v 1.19 2007-05-31 16:42:26 cnoirot Exp $
// ------------------------------------------------------------------
// File:     Sensor.NStart.h
// Contents: Sensor NStart
// ------------------------------------------------------------------

#ifndef  SENSOR_NSTART_H_INCLUDED
#define  SENSOR_NSTART_H_INCLUDED

#include "../../Sensor.h"

/*************************************************************
 **                    SensorNStart                         **
 *************************************************************/
class SensorNStart : public Sensor
{
 private:
  int PositionGiveInfo;

  std::vector<int>    vPosF, vPosR;
  std::vector<double> vValF, vValR;

  int indexF, indexR;
  double startP, startB;
  
  void ReadNStartF (char *, int);
  void ReadNStartR (char *, int);
  void ReadNStartGff3 (char *, int);

 public:
  SensorNStart   (int n, DNASeq *X);
  virtual ~SensorNStart   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorNStart* builder0( int n, DNASeq *X) { return new SensorNStart(n, X);}

#endif
