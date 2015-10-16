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
// $Id: Sensor.AnnotaStruct.h,v 1.10 2010-02-01 16:18:41 sallet Exp $
// ------------------------------------------------------------------

#ifndef  SENSOR_ANNOTASTRUCT_H_INCLUDED
#define  SENSOR_ANNOTASTRUCT_H_INCLUDED

#include <algorithm>
#include "../../Sensor.h"


/*************************************************************
 **                      SensorAnnotaStruct
 *************************************************************/
class SensorAnnotaStruct : public Sensor
{
 private:
  std::vector <Signals*>  vSig;
  std::vector <Contents*> vCon;
  int   PosSigGiveInfo, PosConGiveInfo;
  int   iSig, iCon;
  char  *fileExt;
  string  transFeatName;
  char  startPAR[20], stopPAR[20],   accPAR[20];
  char  donPAR[20],   tStartPAR[20], tStopPAR[20];
  char tStartNpcPAR[20], tStopNpcPAR[20];
  float exonPAR,      intronPAR,     cdsPAR, npcRnaPAR;
  int  exonInline,      intronInline,     cdsInline, npcRnaInline;
  int  startInline,      stopInline,     accInline;
  int  donInline,      tStartInline,     tStopInline;
  int  tStartNpcInline, tStopNpcInline;
  void ReadAnnotaStruct(char[FILENAME_MAX+1], int seqlen);
  void FillOntologyTerm(GeneFeatureSet & geneFeatureSet);
  void ReadAnnotaStructGff3(GeneFeatureSet & geneFeatureSet, int len);
  void PushInCon(int d,   int e,     float *sc,
		 char st, char p[2], int f);
  char * GetScoreC (int type , float scF, bool inlineScore);
 public:
  SensorAnnotaStruct          (int n, DNASeq *X);
  virtual ~SensorAnnotaStruct ();
  virtual void Init           (DNASeq *);
  virtual void GiveInfo       (DNASeq *X, int, DATA *);
  virtual void Plot           (DNASeq *X);
  virtual void PostAnalyse    (Prediction *, FILE *);
};

extern "C" SensorAnnotaStruct* builder0( int n, DNASeq *X) { return new SensorAnnotaStruct(n, X);}

#endif
