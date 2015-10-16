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
// $Id: Sensor.NcRNA.h,v 1.3 2010-02-01 16:18:41 sallet Exp $
// ------------------------------------------------------------------
// File:     Sensor.NcRNA.h
// Contents: Sensor NcRNA
// ------------------------------------------------------------------

#ifndef  SENSOR_NCRNA_H_INCLUDED
#define  SENSOR_NCRNA_H_INCLUDED

#include "../../Sensor.h"


/*************************************************************
 **                       SensorNcRNA                         **
 *************************************************************/
class SensorNcRNA : public Sensor
{
  private:
  char  *fileExt; // File name extension
  std::vector <Contents*> vCon;
  std::vector <Signals*>  vSig;
  char tStartNpc[20], tStopNpc[20], npcRna[20];
  int   iSig, iCon;

  void ReadNcRNAGff3(GeneFeatureSet & geneFeatureSet, int len);

 public:
  SensorNcRNA            (int n, DNASeq *X);
  virtual ~SensorNcRNA      ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorNcRNA * builder0( int n, DNASeq *X) { return new SensorNcRNA(n, X);}

#endif
