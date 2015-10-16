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
// $Id: WAM.h,v 1.4 2008-07-02 11:19:44 tschiex Exp $
// ------------------------------------------------------------------
// File:     WAM.h
// Contents: class WAM
// ------------------------------------------------------------------

#ifndef WAM_INCLUDED
#define WAM_INCLUDED

#include <vector>

#include "markov.h"

#ifndef FILENAME_MAX
#define FILENAME_MAX        1024
#endif

#define TPFILESUFFIX ".TP."
#define FPFILESUFFIX ".FP."
#define SUFFIXLENGTH 4

/*************************************************************
 **                 WAM                    **
 *************************************************************/

class WAM
{
 private:
  int MarkovianOrder;
  int MotifLength;
  Chaine* Alphabet;

  // two vectors of markovian models  
  std::vector<TabChaine<Chaine, unsigned short int>*> TPMOD; 
  std::vector<TabChaine<Chaine, unsigned short int>*> FPMOD;

 public:
  WAM();
  WAM(int order, int length, const char* alphabet, char* prefixfilename);
  ~WAM();
  double ScoreTheMotif(char* oligont);  // sum of likelihood ratio for each position of the entire motif
};

#endif
