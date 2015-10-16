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
// $Id: BackP.h,v 1.23 2007-10-02 07:55:03 tschiex Exp $
// ------------------------------------------------------------------
// File:     BackP.h
// Contents: Definitions for a linear time shortest path with constraint alg.
// Comments: Une liste doublement chainee et circulaire de point de retour
// arriere pour le plus court chemin avec contrainte de duree minimum.
// Le premier element Path (dans Track) ne doit pas etre
// efface. Path.Next est le dernier insere, Path.Prev est le plus vieux
// ??? A t'on vraiment besoin de Prev ?
// ------------------------------------------------------------------

#ifndef  BACKP_H_INCLUDED
#define  BACKP_H_INCLUDED


#include "Const.h"
#include "Prediction.h"
#include "PenaltyDist.h"

class BackPoint
{
  friend class Track;
 private:
  char Status;
 public:
  char IsOptimal();
  void SetOptimal();
  void ClearOptimal();
  char IsMarked();
  void SetMark();
  void ClearMark();
  
  signed char State;
  int StartPos;  
  double Cost;
  double Additional;
  BackPoint *Origin;
  
  BackPoint *Next;
  BackPoint *Prev;    
  
  BackPoint ();
  BackPoint  (char state, int pos, double cost);
  ~BackPoint();

  inline void InitState(int state, int start) { State = state; StartPos = start;};
  void Print();
};


class Track 
{
 public: 

  static unsigned int NumBPAlloc; 
  static unsigned int NumBPCollect;

  BackPoint Path;
  PenaltyDist *PenD;
  double Optimal;
  int OptPos;

  inline void LoadPenalty(char* name) { PenD = new PenaltyDist; PenD->LoadPenaltyDist(name); };
  inline void Update(double cost) {  Path.Next->Additional += cost; Optimal += cost;};
  inline void PayTheSlope() {Path.Next->Additional -= PenD->FinalSlope;
                             Optimal  -= PenD->FinalSlope;  };
  void InsertNew(char state,  int pos, double cost, BackPoint *Or);
  void ForceNew(char state, int pos, double cost, BackPoint *Or);
  BackPoint *BestUsable(int pos, double *cost, int pen = 1);
  Prediction* BackTrace(int From, int To, int Forward = 1);
  void Dump();
  void Zap();
  void ClearMark(int);
  void Mark(int);
  void Sweep(int);

  Track ();
  Track(Track *other);
  ~Track();
  
};

#endif
