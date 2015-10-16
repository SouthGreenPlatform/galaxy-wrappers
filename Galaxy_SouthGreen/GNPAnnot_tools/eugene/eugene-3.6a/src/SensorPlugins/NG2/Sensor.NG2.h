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
// $Id: Sensor.NG2.h,v 1.18 2007-05-31 16:39:35 cnoirot Exp $
// ------------------------------------------------------------------
// File:     Sensor.NG2.h
// Contents: Sensor NG2
// ------------------------------------------------------------------

#ifndef  SENSOR_NG2_H_INCLUDED
#define  SENSOR_NG2_H_INCLUDED

#include "../../Sensor.h"

/*************************************************************
 **                       SensorNG2                         **
 *************************************************************/
class SensorNG2 : public Sensor
{
 private:
  int PositionGiveInfo;

  std::vector<int>    vPosAccF, vPosAccR, vPosDonF, vPosDonR;
  std::vector<double> vValAccF, vValAccR, vValDonF, vValDonR;

  int iAccF, iAccR, iDonF, iDonR;
  double accB, accP, donB, donP;
  
  void ReadNG2F(char[FILENAME_MAX+1], int);
  void ReadNG2R(char[FILENAME_MAX+1], int);
  void ReadNG2Gff3(char[FILENAME_MAX+1], int);
 public:
  SensorNG2  (int n, DNASeq *X);
  virtual ~SensorNG2      ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorNG2 * builder0( int n, DNASeq *X) {  return new SensorNG2(n, X);}


#endif
