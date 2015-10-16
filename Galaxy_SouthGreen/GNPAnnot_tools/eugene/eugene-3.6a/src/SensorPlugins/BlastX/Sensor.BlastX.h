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
// $Id: Sensor.BlastX.h,v 1.29 2009-09-08 08:30:49 sallet Exp $
// ------------------------------------------------------------------
// File:     Sensor.BlastX.h
// Contents: Sensor BlastX
// ------------------------------------------------------------------


#ifndef  SENSOR_BLASTX_H_INCLUDED
#define  SENSOR_BLASTX_H_INCLUDED

#include "../../Sensor.h"
#include "../../Hits.h"
#include "../../State.h"

/*************************************************************
 **                     SensorBlastX                        **
 *************************************************************/
class SensorBlastX : public Sensor
{
 private:
  float *ProtMatch, *ProtMatchLevel;
  int    *ProtMatchPhase;
  std::vector<int>    vPos,     vPMPhase;
  std::vector<float> vPMatch;
  std::vector<int>::iterator iter;
  char  *levels;
  char  *intronlevels;
  int    index;
  float keyBXLevel[10];
  int    minIn;
  int    blastxM;
  int    ppNumber;
  int    stepid;
  int    N;
  int    sloppy;

  void LoadContentScore (DNASeq *);
  char ph06             (char);

  // For postprocess 2
  Hits **HitTable;
  int  NumProt;
  void ProtSupport (Prediction *, FILE *, int debut, int fin,
		    Hits **HitTable,  int Size,  int NumG);

  int  LenSup      (Hits **HitTable, State* state, unsigned char* Sup,
                    std::vector<int> vSupProtI, int& additionalsup,
                    int index, int beg, int end);

  void Print (char name[FILENAME_MAX+1]);

 public:
  SensorBlastX            (int n, DNASeq *X);
  virtual ~SensorBlastX   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorBlastX * builder0( int n, DNASeq *X) { return new SensorBlastX(n, X); }
#endif
