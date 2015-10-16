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
// $Id: Sensor.GSplicer.h,v 1.10 2007-05-31 16:41:02 cnoirot Exp $
// ------------------------------------------------------------------
// File:     Sensor.GSplicer.h
// Contents: Sensor GSplicer  
// ------------------------------------------------------------------


#ifndef  SENSOR_GSPLICER_H_INCLUDED
#define  SENSOR_GSPLICER_H_INCLUDED

#include "../../Sensor.h"

/*************************************************************
 **                      SensorGSplicer
 *************************************************************/
class SensorGSplicer : public Sensor
{
 private:
  int PositionGiveInfo;
  
  std::vector<int>    vPosAccF, vPosAccR, vPosDonF, vPosDonR;
  std::vector<double> vValAccF, vValAccR, vValDonF, vValDonR;

  int iAccF, iAccR, iDonF, iDonR;
  double coefAcc, penAcc, coefDon, penDon;
  
  void   ReadGSplicer(char[FILENAME_MAX+1], int);
  void   ReadGSplicerGff3(char [FILENAME_MAX+1], int );
  double Norm(double, double);

 
 public:
  SensorGSplicer          (int n, DNASeq *X);
  virtual ~SensorGSplicer ();
  virtual void Init       (DNASeq *X);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorGSplicer* builder0( int n, DNASeq *X) { return new SensorGSplicer(n, X);}

#endif
