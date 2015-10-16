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
// $Id: DNASeq.h,v 1.15 2009-07-07 15:51:38 tschiex Exp $
// ------------------------------------------------------------------
// File:     DNASeq.cc
// Contents: Definitions for a class representing DNA seq with ambiguous data
//  and 6 possible unambiguous versions (6 phases possible completions)
// ------------------------------------------------------------------


#ifndef  DNASEQ_H_INCLUDED
#define  DNASEQ_H_INCLUDED

#include <stdio.h>

const unsigned short MASKSEQ = 0x000F;

const unsigned short BitT = 3;
const unsigned short BitG = 2;
const unsigned short BitC = 1;
const unsigned short BitA = 0;


const unsigned short CodeT = (1 << BitT);
const unsigned short CodeG = (1 << BitG);
const unsigned short CodeC = (1 << BitC);
const unsigned short CodeA = (1 << BitA);

class  DNASeq
{
 private:
  int  Size;
  unsigned short *Sequence;

  void UpdateMarkov();

 public:
  DNASeq ();
  DNASeq (int size);
  DNASeq (char* filename);

  ~DNASeq ();

  int  SeqLen;
  char *Name;

  double Markov0[4];
  
  void  Print(FILE *);
  void  PrintTr(FILE*,int, int,signed char);

  // Degeneracy returns the number of possible completely known sequence
  // represented by a degenerated subsequence
  unsigned char Degeneracy(int i, int sens, int len);

  // bit vector values for spliced stop detection
  static const int isTf  = 1;
  static const int isTGf = 2;
  static const int isTAf = 4;
  static const int isTr  = 8;
  static const int isTGr = 16;
  static const int isTAr = 32;

  static const int isGf  = 1;
  static const int isGAf = 2;
  static const int isAf  = 4;
  static const int isARf = 8;
  static const int isGr  = 16;
  static const int isGAr = 32;
  static const int isAr  = 64;
  static const int isARr = 128;

  // see comments in .cc, used to detect possible spliced stops
  // returns a bit vector inside an int.
  int    IsStartStop(int i);
  int    IsStopStop(int i);

  // returns the fraction of completely known sequences that represent the
  // corresponding signal. Example: IsStop on TGN returns 0.25 (only TGA is a Stop)
  double IsAcc(int pos,int strand);
  double IsDon(int pos,int strand);
  double IsStop(int pos,int strand);
  double IsStart(int pos,int strand); // prokaryotic case
  double IsEStart(int pos,int strand);
  
  // Computes the markov probability of emission of the nuc. at position pos
  double Markov(int pos);
  // same for the reverse strand.
  double MarkovR(int pos);
  // returns the GC or AT% depending on the nuc. at position pos (GC or AT).
  double GC_AT(int pos);

  unsigned short Nuc2Code(char Ch);

  // copy from Pos to Len to To. mode is ambiguous or not. Len < 0 means reverse complement.
  void Transfer(int Pos, int Len, char *To, int mode);
  // compute the Frame (as in blast). Strand is '+' or anything else.
  int Pos2Frame(int pos, char strand);
  // remove ambiguity in the nuc. at pos
  unsigned short Unambit (int pos);

  char operator [] (int i); // reads DNA char at position i (forward)
  char operator () (int i); // reads DNA char at position i (complement)
  unsigned short operator () (int i, int mode); // reads code at position i (mode: 0/forward, 1/complement, 2/revcomp)
  char AA(int i, int mode); // translate the codon at position i, using mode as above.
  char nt(int i, int mode); // dead code ?
};
#endif
