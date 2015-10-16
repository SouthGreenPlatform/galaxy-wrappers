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
// $Id: Sensor.cc,v 1.13 2010-02-01 16:18:42 sallet Exp $
// ------------------------------------------------------------------
// File:     Sensor.cc
// Contents: class Sensor
// ------------------------------------------------------------------

#include "Sensor.h"

/*************************************************************
 **                        Sensor                           **
 *************************************************************/
//length of horizontal section of start/den region delimiters
const int HLen = 30;
//colors for phased similarity levels
const int LevelColor[3] = {6,7,8}; 
// ----------------------
//  Default constructor.
// ----------------------
Sensor :: Sensor (int n)
{
  instanceNumber = n;
  type = Type_None;
}

// ----------------------
//  Default destructor.
// ----------------------
Sensor :: ~Sensor ()
{
}

// ---------------------------------------
//  Sanity check : check that ATG occurs.
// ---------------------------------------
void Sensor :: CheckStart(DNASeq *X, std::vector<int> vPosF, std::vector<int> vPosR)
{
  for (int i = 0; i<(int)vPosF.size(); i++)
    if (!X->IsEStart(vPosF[i],1))
      fprintf(stderr,"WARNING: Non ATG start predicted at %d on + strand!\n", vPosF[i]+1);
  
  for (int i = 0; i<(int)vPosR.size(); i++)
    if (!X->IsEStart(vPosR[i]-1,-1))
      fprintf(stderr,"WARNING: Non ATG start predicted at %d on - strand!\n", vPosR[i]+1);
}

// ---------------------------------------
//  Sanity check : stop site check.
// ---------------------------------------
void Sensor :: CheckStop(DNASeq *X, std::vector<int> vPosF, std::vector<int> vPosR)
{
  for (int i = 0; i<(int)vPosF.size(); i++)
    if (!X->IsStop(vPosF[i]-3,1))
      fprintf(stderr,"WARNING: Non T{AA,AG,GA} stop predicted at %d on + strand!\n", vPosF[i]+1);
  
  for (int i = 0; i<(int)vPosR.size(); i++)
    if (!X->IsStop(vPosR[i]+2,-1))
      fprintf(stderr,"WARNING: Non T{AA,AG,GA} stop predicted at %d on - strand!\n", vPosR[i]+1);
}

// -----------------------------------------------
//  Sanity check; AG/GT or GC  splice site check.
// -----------------------------------------------
void Sensor :: CheckSplices (DNASeq *X,
			     std::vector<int> vPosAccF, std::vector<int> vPosDonF,
			     std::vector<int> vPosAccR, std::vector<int> vPosDonR)
{
  for (int i = 0; i<(int)vPosAccF.size(); i++)
    if (X->IsAcc(vPosAccF[i]-2,1)==0.0)
      fprintf(stderr,"WARNING: Non AG (%c%c) acceptor at %d (+ strand) !\n", 
	      (*X)[vPosAccF[i]-2], (*X)[vPosAccF[i]-1], vPosAccF[i]);
  
  for (int i = 0; i<(int)vPosAccR.size(); i++)
    if(X->IsAcc(vPosAccR[i]+1,-1)==0.0)
      fprintf(stderr,"WARNING: Non AG (%c%c) acceptor at %d (- strand) !\n",
	      (*X)(vPosAccR[i]), (*X)(vPosAccR[i]+1),vPosAccR[i]);
  
  for (int i = 0; i<(int)vPosDonF.size(); i++)
    if(X->IsDon(vPosDonF[i],1)==0.0)
      fprintf(stderr,"WARNING: Non GT/GC (%c%c) donor at %d (+ strand) !\n", 
	      (*X)[vPosDonF[i]],(*X)[vPosDonF[i]+1], vPosDonF[i]);
  
  for (int i = 0; i<(int)vPosDonR.size(); i++)
    if(X->IsDon(vPosDonR[i]-1,-1)==0.0)
      fprintf(stderr,"WARNING: Non GT/GC (%c%c) donor at %d (- strand) !\n", 
	      (*X)(vPosDonR[i]-1),(*X)(vPosDonR[i]-2),vPosDonR[i]);
}
// -----------------------------------------------
//  Plot: STOPS
// -----------------------------------------------
void Sensor :: PlotStop(int pos,int phase, int sure)
{
  PlotBarF(pos,phase,0.1,0.2,(sure ? 1 : 9));
}
// -----------------------------------------------
//  Plot: STARTS
// -----------------------------------------------
void Sensor :: PlotStart(int pos,int phase, double strength)
{
    PlotBarF(pos, phase, 0.5, strength, 2);
}
// -----------------------------------------------
//  Plot: Acceptors
// -----------------------------------------------
void Sensor :: PlotAcc(int pos, int strand, double strength)
{
    PlotBarF(pos, (strand >= 0? 4 :-4), 0.5, strength, 4);
}
// -----------------------------------------------
//  Plot: Donor
// -----------------------------------------------
void Sensor :: PlotDon(int pos, int strand, double strength)
{
  PlotBarF(pos, (strand >= 0? 4 :-4), 0.5, strength, 11);
}
// -----------------------------------------------
//  Plot: transcript similarity
// -----------------------------------------------
void Sensor :: PlotESTHit(int start, int end, int strand, int filtered)
{
  if (filtered)
    for (int i = start; i <= end; i++)
      PlotBarI(i, (strand > 0 ? 4 : -4),0.9,1,7);
  else
    for (int i = start; i <= end; i++)
      PlotBarI(i, (strand > 0 ? 4 : -4),0.6,1,2);
}
// -----------------------------------------------
//  Plot: transcript similarity
// -----------------------------------------------
void Sensor :: PlotESTGap(int start, int end, int strand, int filtered)
{
  if (filtered)
    PlotLine(start,end,(strand > 0 ? 4 : -4),(strand > 0 ? 4 : -4),0.9,0.9,7);
  else
    PlotLine(start,end,(strand > 0 ? 4 : -4),(strand > 0 ? 4 : -4),0.6,0.6,2);
}
// -----------------------------------------------
//  Plot: Phased similarity hit
// -----------------------------------------------
void Sensor :: PlotBlastHit(int start, int end, int phase, int level)
{  
    for (int i = start; i <= end; i++)
      PlotBarI(i, phase,0.6+(level/8.0),1,LevelColor[level]);
}
// -----------------------------------------------
//  Plot: Phased similarity gap
// -----------------------------------------------
void Sensor :: PlotBlastGap(int start, int phase1, int end, int phase2, int level)
{
  PlotLine(start,end,phase1,phase2,0.6+(level/8.0),0.6+(level/8.0),LevelColor[level]);
}
// -----------------------------------------------
//  Plot: Start/End region corners
// -----------------------------------------------
void Sensor :: PlotStartReg(int pos, int phase, int color)
{
    PlotBarF(pos, phase, 0.9, 0.2, color);
    PlotLine(pos, pos+HLen, phase, phase, 1.0, 1.0, color);
}
// -----------------------------------------------
//  Plot: Start/End region corners
// -----------------------------------------------
void Sensor :: PlotEndReg(int pos, int phase, int color)
{
    PlotBarF(pos, phase, 0.9, 0.2, color);
    PlotLine(pos, pos-HLen, phase, phase, 1.0, 1.0, color);
}
// -----------------------------------------------
//  Plot: Repeats
// -----------------------------------------------
void Sensor :: PlotRepeat(int start, int end)
{
  for (int i = start; i<= end; i++)
    PlotBarI(i, 0, 0.25, 2, 6);
}


/*************************************************************
 **                     Signals object                      **
 *************************************************************/
// -------------------------
//  Default constructor.
// -------------------------
Signals :: Signals ()
{
  pos   = -1;
  type  = -1;
  edge  = -1;
  score = NULL; 
}

// -------------------------
//  Constructor.
// -------------------------
Signals :: Signals (int p, int t, int e, char *s)
{
  pos   = p;
  type  = t;
  edge  = e;
  score = s;
}


// -------------------------
//  Default destructor.
// -------------------------
Signals :: ~Signals () {}

// -------------------------
//  Print
//--------------------------
void Signals :: PrintS ()
{
  char t[7];
  char s = '+';

  switch (type) {
  case DATA::tStart    : strcpy(t, "tStart   "); break;
  case DATA::tStop     : strcpy(t, "tStop    "); break;
  case DATA::Start     : strcpy(t, "Start    "); break;
  case DATA::Stop      : strcpy(t, "Stop     "); break;
  case DATA::Don       : strcpy(t, "Acc      "); break;
  case DATA::Acc       : strcpy(t, "Don      "); break;
  case DATA::Ins       : strcpy(t, "Ins      "); break;
  case DATA::Del       : strcpy(t, "Del      "); break;
  case DATA::tStartNpc : strcpy(t, "tStartNpc"); break;
  case DATA::tStopNpc  : strcpy(t, "tStopNpc "); break;
  }
  if (edge) s = '-';
  fprintf(stdout, "%d\t%s %c %s\n", pos, t, s, score);
}


/*************************************************************
 **                    Contents object                      **
 *************************************************************/
// -------------------------
//  Default constructor.
// -------------------------
Contents :: Contents ()
{
  start = -1;
  end   = -1;
  type  = -1;
  score = NULL ;
}

// -------------------------
//  Default constructor.
// -------------------------
Contents :: Contents (int sta, int e, int t, float *s)
{
  start = sta;
  end   = e;
  type  = t;
  score = s;
}

// -------------------------
//  Default destructor.
// -------------------------
Contents :: ~Contents () {}

// -------------------------
//  Print
//--------------------------
void Contents :: PrintC ()
{
  char t[11];

  switch (type) {
  case DATA::ExonF1     : strcpy(t, "ExonF1    "); break;
  case DATA::ExonF2     : strcpy(t, "ExonF2    "); break;
  case DATA::ExonF3     : strcpy(t, "ExonF3    "); break;
  case DATA::ExonR1     : strcpy(t, "ExonR1    "); break;
  case DATA::ExonR2     : strcpy(t, "ExonR2    "); break;
  case DATA::ExonR3     : strcpy(t, "ExonR3    "); break;
  case DATA::IntronF    : strcpy(t, "IntronF   "); break;
  case DATA::IntronR    : strcpy(t, "IntronR   "); break;
  case DATA::InterG     : strcpy(t, "InterG    "); break;
  case DATA::UTR5F      : strcpy(t, "UTR5F     "); break;
  case DATA::UTR5R      : strcpy(t, "UTR5R     "); break;
  case DATA::UTR3F      : strcpy(t, "UTR3F     "); break;
  case DATA::UTR3R      : strcpy(t, "UTR3R     "); break;
  case DATA::IntronUTRF : strcpy(t, "IntronUTRF"); break;
  case DATA::IntronUTRR : strcpy(t, "IntronUTRR"); break;
  case DATA::RNAF       : strcpy(t, "RnaF      "); break;
  case DATA::RNAR       : strcpy(t, "RnaR      "); break;
  }
  fprintf(stdout, "%d\t%d\t%s %f\n", start, end, t, *score);
}


