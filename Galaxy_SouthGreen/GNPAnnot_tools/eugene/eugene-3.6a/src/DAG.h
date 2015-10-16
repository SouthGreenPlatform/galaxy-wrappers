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
// $Id: DAG.h,v 1.15 2010-01-25 17:03:05 sallet Exp $
// ------------------------------------------------------------------
// File:     DAG.h
// Contents: Modeling genes as a DAG and exploring it - Header
// ------------------------------------------------------------------

#ifndef DAG_H_INCLUDED
#define DAG_H_INCLUDED
#endif

#include "Const.h"
#include "Prediction.h"
#include "BackP.h"
#include "SensorIF.h"
#include "Param.h"
#include "DNASeq.h"
 
class DAG
{

   static DNASeq *TheSeq;
   static double ExPrior, InPrior, IGPrior, FivePrior, ThreePrior, IntronFivePrior, RnaPrior;
   static double SplicedStopPen;
   static int estuse;
   static double NormalizingPath;
   static double PBest[NbTracks];
   static BackPoint *PrevBP[NbTracks];

private: 

  int StartPosition;
  int EndPosition;
  Track  LBP[NbTracks];
  	   
 public:
  char EvidenceName[FILENAME_MAX+1];//   **
  Prediction *pred; //   **

  DAG();
  DAG (int start, int end, Parameters &PAR, DNASeq *Seq);
  DAG (int start, int end, DAG *RefDag,char* name);
  ~DAG ();
  void WeightThePrior();
  void LoadDistLength();
  void ShortestPathAlgoForward (int position, DATA Data);
  void ShortestPathAlgoBackward (int position, DATA Data, int NoContentsUpdate=0);
  void MarkAndSweep(int,int,int);
  //  int ChooseTheBestTrackIndex (int position, double *bestscore, int sense);
  //  int ChooseTheBestTrackIndex (int position, double *bestscore, DAG* DagR);
  //  void InsertLastBestUsableBP (int position, int sense=1);
  //  void SubstractContents (DATA *data);
  //  double ConnectDag (DAG *RightDag, int position, DATA * data, int datalength);
  double BuildPrediction (int From, int To, int Forward = 1);
  //  void CleanPrediction (DAG* dag, DAG* dagrev);

 private:
  void Normalize();
  void ComputeSigShifts(enum Signal::Edge Strand, DATA Data, int position);
  void ComputeRequired(Signal::Edge, DATA, int);
  void ApplyScore(int position, DATA Data, int NoContentsUpdate);  
  void ApplyLengthPenalty(int position, DATA Data, int NoContentsUpdate);  

  inline int GetStart() { return StartPosition; };
  void Print();
  void StatActive();
  void StatGC();
};

