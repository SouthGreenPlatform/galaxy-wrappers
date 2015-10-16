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
// $Id: Sensor.SMachine.h,v 1.9 2007-05-04 13:28:04 cnoirot Exp $
// ------------------------------------------------------------------
// File:     Sensor.SMachine.h
// Contents: Sensor SMachine
// ------------------------------------------------------------------

#ifndef  SENSOR_SMACHINE_H_INCLUDED
#define  SENSOR_SMACHINE_H_INCLUDED

#include "../../Sensor.h"


/*************************************************************
 **                       SensorSMachine                         **
 *************************************************************/
class SensorSMachine : public Sensor
{
 private:
  int PositionGiveInfo;

  std::vector<int>    vPosAccF, vPosAccR, vPosDonF, vPosDonR;
  std::vector<double> vValAccF, vValAccR, vValDonF, vValDonR;

  std::vector<int>    vPosF, vPosR;
  std::vector<double> vValF, vValR;

  int iAccF, iAccR, iDonF, iDonR;
  double accB, accP, donB, donP;
  double transSpliceB;
  int isScaled;

  int indexF, indexR;
  double startP, startB;

  void ReadSMachineSplices(char *, int);
  void ReadSMachineStarts(char *, int);

  void SpliceMachine();
  void Print (char name[FILENAME_MAX+1]);
  void ReadMachineGff3(char *, int);
 public:
  SensorSMachine  (int n, DNASeq *X);
  virtual ~SensorSMachine ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorSMachine * builder0( int n, DNASeq *X) {  return new SensorSMachine(n, X);}


#endif
