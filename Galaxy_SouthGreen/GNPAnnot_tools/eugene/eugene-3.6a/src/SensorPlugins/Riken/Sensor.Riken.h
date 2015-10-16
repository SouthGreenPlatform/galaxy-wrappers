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
// $Id: Sensor.Riken.h,v 1.18 2007-07-13 13:37:16 cnoirot Exp $
// ------------------------------------------------------------------
// File:     Sensor.Riken.h
// Contents: Sensor Riken
// RAFL: Riken Arabidopsis Full Length cDNA
// ------------------------------------------------------------------

#ifndef  SENSOR_RIKEN_H_INCLUDED
#define  SENSOR_RIKEN_H_INCLUDED

#include <vector>
#include <algorithm>

#include "../../Sensor.h"

// ************
// * RAFLgene *     // RAFL: Riken Arabidopsis Full Length cDNA 
// ************
class RAFLgene
{
 private:

 public:
  int deb;
  int fin;
  signed char sens;
  char ID[FILENAME_MAX+1];

  RAFLgene  ();
  ~RAFLgene ();
};

/*************************************************************
 **                     SensorRiken                         **
 *************************************************************/
class SensorRiken : public Sensor
{
 private:
  std::vector <RAFLgene> RAFL;
  int RAFLpos;                  // Position par rapport a un gene RAFL
  int RAFLindex;                // Index du RIKEN en cours
  int StrandRespect;
  int MIN_EST_DIFF;         // default = 100;
  int MAX_OVERLAP;          // default = 60;
  int MAX_RIKEN_LENGTH;     // default = 60000;
  int MAX_RIKEN_EST_LENGTH; // default  = 3000;
  int MIN_RIKEN_LENGTH;     // default  = 120; // 2* riken overlap (60)
  int MIN_RIKEN_EST_LENGTH; // default  = 10;
  double RAFLPenalty;       // default  = -120;
  int checkRAFL ( RAFLgene & tmp, int beg5, int end5, int beg3, int end3, int len);
 public:
  SensorRiken   (int n, DNASeq *X);
  void readRiken (char name[FILENAME_MAX+1],DNASeq *X, std::vector <RAFLgene> & RAFLtmp);
  void readRikenGff3 (GeneFeatureSet & geneFeatureSet ,DNASeq *X, std::vector <RAFLgene> & RAFLtmp);
  virtual ~SensorRiken    ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorRiken * builder0( int n, DNASeq *X) {  return new SensorRiken(n, X);}

#endif
