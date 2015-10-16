// ------------------------------------------------------------------
// Copyright (C) 2005 INRA <eugene@ossau.toulouse.inra.fr>
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
// $Id: AltEst.h,v 1.5 2009-08-28 12:48:45 sallet Exp $
// ------------------------------------------------------------------
// File:     AltEst.h
// Contents: class Alternative Est
// ------------------------------------------------------------------

#ifndef ALTEST_H_INCLUDED
#define ALTEST_H_INCLUDED

#include <vector>
#include <algorithm>
#include <string.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif

#include "Const.h"
#include "System.h"
#include "Param.h"
#include "Prediction.h"
#include "SensorIF.h"
#include "State.h"
/*************************************************************
 **                      OneAltEst
 *************************************************************/
class OneAltEst
{
  friend class AltEst;

 private:
  char id[FILENAME_MAX+1];
  int  start, end, index;
  int  exonsNumber, totalLength;
  bool altSplicingEvidence;
  std::vector <int> vi_ExonStart;
  std::vector <int> vi_ExonEnd;
  void Penalize(int pos, DATA *Data, double altPenalty);
  
 public:
  OneAltEst  ();
  OneAltEst  (char* id, int i, int j);
  ~OneAltEst ();
  void Reset             ();
  void AddExon           (int, int);
  void RemoveExon        (int  i);
  void ExtremitiesTrim   (int  exonucleasicLength);
  int  IsFiltered        (bool unspliced, bool extremelen, bool verbose,
			  int  minIn,     int  maxIn,      int  maxEx,   int minEx,
			  int  minEstLen, int  maxEstLen);
  bool IsInconsistentWith(OneAltEst*);
  bool CompatibleWith(Prediction *pred);
  void Print           ();
  inline void  UpdateBoundaries() { start = vi_ExonStart[0]; end = vi_ExonEnd[vi_ExonEnd.size()-1]; };
  inline char* GetId()            { return id;  };
  inline int   GetEnd()           { return end; };
  inline int GetStart()           {return start; };
  inline int   GetAltSplE()       { return altSplicingEvidence; };
  inline void  PutAltSplE(bool b) { altSplicingEvidence = b;    };
  inline void PutIndex(int i) { index = i;};
};

/*************************************************************
 **                      AltEst
 *************************************************************/
class AltEst
{
 private:
  bool altEstDisplay,  verbose;
  int  minIn, maxIn, maxEx, minEx;
  int  minEstLength, maxEstLength, exonucleasicLength;

  double altPenalty;
  int includedEstFilter;
  int compatibleEstFilter;
  int unsplicedEstFilter;
  int extremeLengthFilter;
  int nextAdd, nextRemove;

  int  ReadAltFile (char[FILENAME_MAX+1], int &nbUnspliced, int &nbExtremLen);
  void Compare     (int &nbIncomp, int &nbNoevidence, int &nbIncluded);

 public:
  int  totalAltEstNumber;
  std::vector<OneAltEst> voae_AltEst;

  AltEst  (DNASeq *X);
  ~AltEst ();
  void Penalize(int i, int pos, DATA *Data);
  int convertHitsToAltEst (Hits * AllEST, int &nbUnspliced, int &nbExtremLen, int & NumEST);
};

#endif
