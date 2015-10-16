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
// $Id: Sensor.GFF.h,v 1.7 2007-08-17 13:20:04 cnoirot Exp $
// ------------------------------------------------------------------
// File:     Sensor.GFF.h
// Contents: Sensor GFF  
// ------------------------------------------------------------------

#ifndef  SENSOR_GFF_H_INCLUDED
#define  SENSOR_GFF_H_INCLUDED

#include "../../Sensor.h"

/*************************************************************
 **                        GFFObject                        **
 *************************************************************/
class GFFObject
{
 private:
  
 public:
  char  name[30];
  char  feature[30];
  int   start, end;
  char  strand, frame;
  int   A, T, C, G, N;
  float GC;
  
  GFFObject  ();
  GFFObject  (char*, char*,int,int,char,char,int,int,int,int,int,float);
  void Print       ();
  void PrintHeader ();
  ~GFFObject ();
};

/*************************************************************
 **                      SensorGFF
 *************************************************************/
class SensorGFF : public Sensor
{
 private:
  std::vector <GFFObject*> gffList;
  static FILE *ppfile;
  static bool IsInitialized;

  void ReadGFF(char[FILENAME_MAX+1], int seqlen);
  void ReadGFF3(GeneFeatureSet & geneFeatureSet , DNASeq *X );

 public:
  SensorGFF               (int n, DNASeq *X);
  virtual ~SensorGFF      ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorGFF* builder0( int n, DNASeq *X) { return new SensorGFF(n, X);}



FILE * SensorGFF::ppfile;
bool SensorGFF::IsInitialized = false;


#endif
