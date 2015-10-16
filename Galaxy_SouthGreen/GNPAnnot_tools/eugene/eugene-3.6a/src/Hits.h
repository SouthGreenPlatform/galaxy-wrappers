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
// $Id: Hits.h,v 1.13 2009-04-29 15:06:18 sallet Exp $
// ------------------------------------------------------------------
// File:     Hits.h
// Contents: Definitions for a class representing alignements with genomic seq.
// ------------------------------------------------------------------

#ifndef  HITS_H_INCLUDED
#define  HITS_H_INCLUDED

#include <stdio.h>
#include "GeneFeatureSet.h"
#include "DNASeq.h"
class Block
{
 private:
  
 public:
  Block ();
  Block(int Start, int End, int LStart, int LEnd, int Ph, int Scr);
  
  ~Block ();
  
  void AddBlockAfter(int Start,int End,int LStart,int LEnd,int Ph,int Scr,char *HSP = NULL);
  void AddBlockAfter(Block* newBlock, char *HSP = NULL);
  int Start;
  int End;
  int Phase;
  int Score;
  int LStart;
  int LEnd;
  char *HitSeq;
  Block *Prev,*Next;
};

class  Hits
{
 private:
  
 public:
  Hits ();
  Hits (char* name, int length, char strand, int deb, int fin, int ldeb,
	int lfin, int Ph, int Scr, double prob, int level, int sup);
  Hits* ReadFromFile(FILE* HitFile, int *NumHits, int level, int margin, int maxPos);
  Hits* ReadFromGeneFeatureSet(GeneFeatureSet& HitSet , int *NumHits, int level, int margin, DNASeq *X);
  ~Hits ();


  char   *Name;
  char   Strand;
  int    Length;
  int    NGaps;
  int    Start;
  int    End;
  double Evalue;
  char   Rejected;
  int    Level;         // For blastx level 0-9 
  int    Support;

  Block *Match;
  Hits  *Next;
  
  char * getName ();
  
};
#endif
