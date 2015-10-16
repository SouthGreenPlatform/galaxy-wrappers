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
// $Id: Sensor.PepSignal.h,v 1.4 2007-05-31 16:43:15 cnoirot Exp $
// ------------------------------------------------------------------
// File:     Sensor.PepSignal.h
// Contents: Sensor Peptide Signal
// ------------------------------------------------------------------

#ifndef  SENSOR_PEPSIGNAL_H_INCLUDED
#define  SENSOR_PEPSIGNAL_H_INCLUDED

#include "../../Sensor.h"


/*************************************************************
 **                    SensorPepSignal                         **
 *************************************************************/
class SensorPepSignal : public Sensor
{
 private:
  int PositionGiveInfo;

  std::vector<int>    vPosF, vPosR;
  std::vector<double> vValF, vValR;

  int indexF, indexR;
  double startP, startB;
  
  void ReadPepSignalStarts (char *, int);
  void ReadPepSignalGff3 (char *, int);

 public:
  SensorPepSignal   (int n, DNASeq *X);
  virtual ~SensorPepSignal   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorPepSignal* builder0( int n, DNASeq *X) { return new SensorPepSignal(n, X);}

#endif
