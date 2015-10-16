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
// $Id: Sensor.SPred.h,v 1.18 2007-05-31 16:47:06 cnoirot Exp $
// ------------------------------------------------------------------
// File:     Sensor.SPred.h
// Contents: Sensor SPred
// ------------------------------------------------------------------

#ifndef  SENSOR_SPRED_H_INCLUDED
#define  SENSOR_SPRED_H_INCLUDED

#include "../../Sensor.h"

/*************************************************************
 **                     SensorSPred                         **
 *************************************************************/
class SensorSPred : public Sensor
{
 private:
  int PositionGiveInfo;

  std::vector<int>    vPosAccF, vPosAccR, vPosDonF, vPosDonR;
  std::vector<double> vValAccF, vValAccR, vValDonF, vValDonR;

  int iAccF, iAccR, iDonF, iDonR;
  double accP, accB, donP, donB;
  
  void ReadSPredF(char[FILENAME_MAX+1], int);
  void ReadSPredR(char[FILENAME_MAX+1], int);
  void ReadSPredGff3(char[FILENAME_MAX+1], int);
  
 public:
  SensorSPred   (int n, DNASeq *X);
  virtual ~SensorSPred    ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorSPred * builder0( int n, DNASeq *X) {  return new SensorSPred(n, X);}

#endif
