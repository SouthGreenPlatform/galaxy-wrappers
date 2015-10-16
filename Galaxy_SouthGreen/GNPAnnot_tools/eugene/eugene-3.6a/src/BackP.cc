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
// $Id: BackP.cc,v 1.29 2009-08-10 15:06:20 sallet Exp $
// ------------------------------------------------------------------
// File:     BackP.cc
// Contents: Definitions for a linear time shortest path with constraint alg.
// ------------------------------------------------------------------

#include <stdio.h>
#include <math.h>
#include "assert.h"

#include "BackP.h"
#include "System.h"

// Bit masks in the Status field
const char OptimalMask = 0x01;
const char MarkedMask  = 0x02;
// ----------------------------------------------------------------
// Status handling
// ----------------------------------------------------------------
inline char BackPoint :: IsOptimal() { return (Status & OptimalMask); }
inline char BackPoint :: IsMarked() { return (Status & MarkedMask); }
inline void BackPoint :: SetOptimal() { Status |= OptimalMask; }
inline void BackPoint :: SetMark() { Status |= MarkedMask; }
inline void BackPoint :: ClearOptimal(){ Status &= (0xff^OptimalMask); }
inline void BackPoint :: ClearMark(){ Status &= (0xff^MarkedMask); }
// ----------------------------------------------------------------
// Default constructor.
// ----------------------------------------------------------------
BackPoint :: BackPoint  ()
{
  State = -1;
  StartPos = 0;
  Cost = 0.0;
  Additional = 0.0;
  Next = Prev = this;
  Origin = NULL;
  Status = 0;
}
// ----------------------------------------------------------------
//  Default constructor.
// ----------------------------------------------------------------
BackPoint :: BackPoint  (char state, int pos, double cost)
{
  State = state;
  StartPos = pos;
  Cost = cost;
  Additional = 0.0;
  Next = Prev = Origin = NULL;
  Status = 0;
}

// ----------------------------------------------------------------
//  Default destructor.
// ----------------------------------------------------------------
BackPoint :: ~BackPoint  ()
{
  Prev->Next = Next;
  Next->Prev = Prev;
}

// ----------------------------------------------------------------
// Prints the BackPoint contents
// ----------------------------------------------------------------
void BackPoint :: Print()
{
  printf("pos = %d, state = %d\n", StartPos,State);
}


// ----------------------------------------------------------------
// Class Track
// STATIC INIT
// ----------------------------------------------------------------
unsigned int Track::NumBPAlloc = 0; 
unsigned int Track::NumBPCollect = 0;
// ----------------------------------------------------------------
// Track creator
// ----------------------------------------------------------------
Track :: Track()
{
  Optimal = 0.0;
  OptPos =0;
  PenD = NULL;
}
// ----------------------------------------------------------------
// Track creator: from another Track
// ----------------------------------------------------------------
Track :: Track(Track *other)
{
  Optimal 		= other->Optimal;
  OptPos 		= other->OptPos;
  PenD 		= other->PenD;
}

// ----------------------------------------------------------------
// Track destructor
// ----------------------------------------------------------------
Track :: ~Track()
{
  Zap();
  delete PenD;
}
// ----------------------------------------------------------------
// Zap  the path data of a  whole track
// ----------------------------------------------------------------
void Track :: Zap()
{
  BackPoint *Dead = Path.Next,*Skip;

  while (Dead != &Path) {
    Skip = Dead->Next;
    delete Dead;
    Dead = Skip;
  }
}
// ----------------------------------------------------------------
// Clear Mark of all BP in a whole track
// ----------------------------------------------------------------
void Track :: ClearMark(int pos)
{
  BackPoint *It = Path.Next;

  while ((It != &Path) && (It->StartPos >= pos)) {
    It->ClearMark();
    It = It->Next;
  }
}
// ----------------------------------------------------------------
// Mark all backpoint from a useful one following Origin
// ----------------------------------------------------------------
void MarkTrace(int pos, BackPoint *Source) {

  while ((Source) &&  (Source->StartPos >= pos) && 
	 (!Source->IsMarked())) {
    Source->SetMark();
    Source = Source->Origin;
  }
}
// ----------------------------------------------------------------
// Mark all useful BP in a whole track
// ----------------------------------------------------------------
void Track :: Mark(int pos)
{
  BackPoint *It = Path.Next;

  do {
    MarkTrace(pos,It);
    if (isinf(It->Additional)) break;

    if (It->IsOptimal() &&  
	(abs(pos-It->StartPos) > PenD->MaxLen)) break;

    It = It -> Next;
  } while (It != Path.Next);
}
// ----------------------------------------------------------------
// Sweep for Mark and Sweep
// ---------------------------------------------------------------- We
// have marked any BP accessible from Origin of "young" (not in the
// linear tail). So remaining BP are old BP unaccessible by any Origin
// backtracing from young BP. They can be deleted once the information
// inside is taken into account in existing BP.
// Among Optimal/Mark/State/StartPos/Cost/additional/Origin/Prev/Next
// only Additional will be "saved" by including it in the Additional
// field of the next BackP structure

void Track :: Sweep(int pos) {

  BackPoint *It = Path.Next;

  while ( (It != &Path)){// && (It->StartPos >= pos)) {
    if (!It->IsMarked()) {
      It->Next->Additional += It->Additional;
      It->Next->Prev = It->Prev;
      It->Prev->Next = It->Next;
      delete It;
      NumBPCollect++;
    }

    It = It->Next;
  }
}
// ----------------------------------------------------------------
// Insert  a new backpoint
// ----------------------------------------------------------------
void Track :: InsertNew(char state, int pos, double cost, BackPoint *Or)
{
  BackPoint *It = Path.Next;
  
  //  if (cost > NINFINITY) {
  //  if (cost > Optimal-PenD->MaxDelta) {
  //  if (cost > Optimal-PenD->GetDelta(pos-OptPos)) {
  //  if (cost > Path.Next->Cost+Path.Next->Additional-PenD->GetDelta(pos-Path.Next->StartPos)) {
  if (cost > Optimal-PenD->GetDelta(abs(pos-OptPos)) && 
      cost > Path.Next->Cost+Path.Next->Additional-PenD->GetDelta(abs(pos-Path.Next->StartPos))) {
    NumBPAlloc++;
    It =  new BackPoint(state,pos,cost);
    if (cost > Optimal) { 
      Optimal = cost; 
      It->SetOptimal();
      OptPos = pos;
    }
    It->Next = Path.Next;
    It->Prev = Path.Next->Prev;
    It->Origin = Or;
    Path.Next->Prev = It;
    Path.Next = It;
  }
}
// ----------------------------------------------------------------
// Insert  a new backpoint
// ----------------------------------------------------------------
void Track :: ForceNew(char state, int pos, double cost, BackPoint *Or)
{
  BackPoint *It = Path.Next;
  
  NumBPAlloc++;
  It =  new BackPoint(state,pos,cost);
  if (cost > Optimal) { Optimal =cost; It->SetOptimal();}
  It->Next = Path.Next;
  It->Prev = Path.Next->Prev;
  It->Origin = Or;
  Path.Next->Prev = It;
  Path.Next = It;
}
// ----------------------------------------------------------------
// Returns the best BackPoint taking into account length penalties
// Update value of *cost
// ----------------------------------------------------------------
BackPoint *Track :: BestUsable(int pos, double *cost, int pen)
{
  BackPoint *BestBP = NULL;
  double BestCost = NINFINITY;
  int Len;
  BackPoint *It = Path.Next;
  double LenPen,Add = 0.0;

  do {
    if (isinf(It->Additional)) break;
    Add += It->Additional;
    Len = abs(pos-It->StartPos);

    // when pen == 0 or Origin is NULL, this means we reach the
    // extremities of the sequence. Therefore we account for an
    // optimistic penality (MinPen) given that the length can only be
    // longer than the actual length. In all cases, we discount the
    // Slope penalty on the length. We further discount one on extremities 
    // (because -1 and Data_Len+1 are used).

    if (pen && It->Origin) 
      LenPen = (*PenD)[Len] - (PenD->FinalSlope)*Len;
    else 
      LenPen = PenD->MinPen(Len) - PenD->FinalSlope*(Len-1);

    if ((Add + It->Cost - LenPen) > BestCost) {
      BestCost = Add+It->Cost - LenPen;
      BestBP = It;
    }
    if (It->IsOptimal() && Len >= PenD->MaxLen) break;

    It = It -> Next;
  } while (It != Path.Next);
  

  *cost = BestCost;
  return BestBP;
}
// ----------------------------------------------------------------
// BackTrace and build a prediction object
// ----------------------------------------------------------------
Prediction* Track :: BackTrace (int From, int To, int Forward)
{
  std::vector <int>         vPos;
  std::vector <signed char> vState;
  Prediction *pred;
  BackPoint  *It;
  int  pos;
  char etat;
  int  prevpos, prevstate;

  // put state back on correct transitions for backward predictions
  if (!Forward) {
    It = Path.Next;
    etat = It->State;
    It = It->Origin;
    
    while (It != NULL) {
      prevstate = etat;
      etat = It->State;
      It->State = prevstate;
      It = It->Origin;
    }
  }

  It = Path.Next;

  // initialisation by the terminal state
  pos  = It->StartPos;
  etat = It->State;
  It   = It->Origin;
  if (pos >= From) {
    vPos.push_back  ( pos  );
    vState.push_back( etat );
  }

  prevpos = pos;
  prevstate = etat;

  //  printf("pos %d etat %d CDSlen %d prevpos %d\n", pos,etat,CDSlen,prevpos);

  while (It != NULL) {

    pos  = It->StartPos;
    etat = It->State;
    It   = It->Origin;

    //    printf("pos %d etat %d CDSlen %d prevpos %d\n", pos,etat,CDSlen,prevpos);

    if (pos >= From) {
      vPos.push_back  ( pos  );
      vState.push_back( etat );
    }
    
    prevpos   = pos;
    prevstate = etat;
  }
  
  if (Forward) {
    vPos[0] -=  1;
    reverse(vState.begin(), vState.end());
    reverse(vPos.begin(),   vPos.end());
  }
  
  pred = new Prediction(From, To, vPos, vState);
  vPos.clear();
  vState.clear();

  return pred;
}
// ----------------------------------------------------------------
// Dumps the contents of all the Backpoints in the path
// ----------------------------------------------------------------
void Track :: Dump ()
{
  BackPoint* It = Path.Next;

  printf("Number of BP allocated: %d\n",NumBPAlloc);
  printf("Number of BP garbage collected: %d\n",NumBPCollect);
  printf("Number of active BP: %d\n\n",NumBPAlloc-NumBPCollect);

  do {
    printf("pos = %d, state = %d, cost = %f%s, additional = %f",
	   It->StartPos,It->State,It->Cost,
	   (It->IsOptimal() ? "*" : ""), It->Additional);
    if (It->Origin)
      printf(" Or.state = %d Or.pos = %d",It->Origin->State,It->Origin->StartPos);
    printf("\n");
    It = It->Next;
  }  while (It != Path.Next);
}
