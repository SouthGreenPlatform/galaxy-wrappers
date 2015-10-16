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
// $Id: DNASeq.cc,v 1.19 2009-07-07 15:51:38 tschiex Exp $
// ------------------------------------------------------------------
// File:     DNASeq.cc
// Contents: Definitions for a class representing DNA seq with ambiguous data
//  and 6 possible unambiguous versions (6 phases possible completions)
// ------------------------------------------------------------------



#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#include <stdlib.h>
#include <assert.h>
#ifdef STDC_HEADERS
#include <string.h>
#else
#include <strings.h>
#endif

#include "DNASeq.h"
#include "Const.h"
#include "System.h"

// ---------------------------------------------------------------------
// Warning. These tables depend on the numerical representation of 
// Code: numeric code of a nucleotide (bit vector)
// Nuc:  corresponding char in the IUPAC representation, lower case
// CNuc: complement IUPAC
// UNuc: unambiguous IUPAC (arbitrary choice, preferring C to avoid stop codons)
// Bit: bit position of the unambiguous code (T is 8 (code) and 3 (bit))
// CCode: complement code
// ---------------------------------------------------------------------

//                            0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
const char  Code2Nuc[16]  = {'-','a','c','m','g','r','s','v','t','w','y','h','k','d','b','n'}; // numeric code to DNA char
const char  Code2CNuc[16] = {'-','t','g','k','c','y','s','b','a','w','r','d','m','h','v','n'}; // numeric code to complement DNA char
const char  Code2UNuc[16] = {'-','a','c','c','g','g','c','c','t','t','c','c','t','g','c','c'}; // numeric code to unambiguous DNA char


const unsigned short Code2UCode[16] = {0    ,CodeA,CodeC,CodeC,
				       CodeG,CodeG,CodeC,CodeC,
				       CodeT,CodeT,CodeC,CodeC,
				       CodeT,CodeG,CodeC,CodeC};

const unsigned short Code2Bit[16]  = {0    ,BitA,BitC,BitC,
				       BitG,BitG,BitC,BitC,
				       BitT,BitT,BitC,BitC,
				       BitT,BitG,BitC,BitC};

const unsigned short Code2CCode[16] = {0    ,CodeT,CodeG,CodeT|CodeG,
				       CodeC,CodeT|CodeC,CodeG|CodeC,CodeT|CodeG|CodeC,
				       CodeA,CodeT|CodeA,CodeG|CodeA,CodeT|CodeG|CodeA,
				       CodeC|CodeA,CodeT|CodeC|CodeA,CodeG|CodeC|CodeA,CodeT|CodeC|CodeA|CodeG};

const unsigned short Code2Code[16] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};

//                            0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15
const int Code2NumOne[16]  = {0,  1,  1,  2,  1,  2,  2,  3,  1,  2,  2,  3,  2,  3,  3,  4}; // number of possible different nucleotides

// Translation table for codons
const unsigned char Trad[64] = {'K','N','K','N','T','T','T','T',
				'R','S','R','S','I','I','M','I',
				'Q','H','Q','H','P','P','P','P',
				'R','R','R','R','L','L','L','L',
				'E','D','E','D','A','A','A','A',
				'G','G','G','G','V','V','V','V',
				'*','Y','*','Y','S','S','S','S',
				'*','C','W','C','L','F','L','F'};

// ---------------------------------------------------------------------
//  Default constructor.
// ---------------------------------------------------------------------
DNASeq :: DNASeq  ()
{
  SeqLen = 0;
  Size = 0;
  Sequence = NULL;
}

// ---------------------------------------------------------------------     
//  Construct a  DNASeq  of size L
// ---------------------------------------------------------------------
DNASeq :: DNASeq(int L)
{
  int  i;

  Sequence = (unsigned short int *) Safe_malloc (L * sizeof(unsigned short int));
  Size = L;
  SeqLen = 0;

  Markov0[BitA] = Markov0[BitT] = Markov0[BitG] = Markov0[BitC] = 0.0;

  for  (i = 0;  i < L;  i ++)  Sequence [i] = 0;
}
// ---------------------------------------------------------------------
// Construct a DNASeq from a filename and store the contents of the
// file into Sequence, Allocate memory as needed. Assumes FASTA
// first. If no > is found assumes raw DNA every unknown DNA char is
// replaced by N
// ---------------------------------------------------------------------
DNASeq :: DNASeq (char *filename)
{
  const unsigned int  INCR_SIZE = 10000;
  const unsigned int  INIT_SIZE = 10000;

  char  *P, Line [MAX_LINE];
  int  Len;
  int  Ch;
  FILE   *fp;  

  fp = (*filename ? FileOpen (NULL, filename, "r") : stdin);
  
  if (fp == NULL) {
    fprintf(stderr, "Cannot open fasta file %s\n", filename);
    exit(3);
  }

  // reach the ">"
  while  ((Ch = fgetc (fp)) != EOF && Ch != '>')
    ;
  
  if  (Ch != EOF) {
    fgets (Line, MAX_LINE, fp);
    Len = strlen (Line);
    assert (Line [Len - 1] == '\n');
    P = strtok (Line, " \t\n");
    if (!P) P = BaseName(filename);
    Len = strlen (P);
    Name = (char *)Safe_malloc(sizeof(char)*(Len+1));
    strcpy (Name, P);
  }
  else {
    // Assume the file is not FASTA but RAW DNA - will not work on stdin
    rewind(fp);
  }
  
  Sequence = (unsigned short int *) Safe_malloc (INIT_SIZE * sizeof(unsigned short int));
  Size = INIT_SIZE;
  
  Len = 0;
  while  ((Ch = fgetc (fp)) != EOF && Ch != '>') {
    if (isspace (Ch))
      continue;
    
    if  (Len >= Size) {
      Size += INCR_SIZE;
      Sequence = (unsigned short int *) Safe_realloc (Sequence, sizeof(unsigned short int)*Size);
    }
    Sequence [Len++] = Nuc2Code (Ch);
  }
  SeqLen = Len;
  
  UpdateMarkov();

  if (fp != stdin) fclose(fp);
  
  return;
}

// ---------------------------------------------------------------------
//  Destroy this DNASeq by freeing its memory.
// ---------------------------------------------------------------------
DNASeq :: ~ DNASeq  ()
{
  free(Sequence);
  free(Name);
}

// ---------------------------------------------------------------------
//  Return the char at position i
//  On the sequence ATG, (*this)[2] returns 'g' 
// ---------------------------------------------------------------------

char DNASeq :: operator [] (int i)
{
  if ((i >= SeqLen) || (i < 0))
    return Code2Nuc[CodeT|CodeC|CodeA|CodeG];
  else
    return  Code2Nuc[Sequence[i] & MASKSEQ];
}
// ---------------------------------------------------------------------
//  Return the complement char at position i
//  On the sequence ATG, (*this)(2) returns 'c' 
// ---------------------------------------------------------------------

char DNASeq :: operator () (int i)
{
  if ((i >= SeqLen) || (i < 0))
    return Code2CNuc[CodeT|CodeC|CodeA|CodeG];
  else
    return  Code2CNuc[Sequence[i] & MASKSEQ];
}
// ---------------------------------------------------------------------
// Reads the numerical ccode at position i.
// 	mode 0 is for forward
// 	mode 1 is for complement
// 	mode 2 is reverse complement (starts from the end)
// ---------------------------------------------------------------------
unsigned short DNASeq :: operator () (int i, int mode)
{
  if ((i >= SeqLen) || (i < 0))
    return CodeT|CodeC|CodeA|CodeG;
  
  switch (mode) {
  case 0:
    return Code2Code[Sequence[i] & MASKSEQ];
  case 1:
    return  Code2CCode[Sequence[i] & MASKSEQ];
  case 2:
    return Code2CCode[Sequence[SeqLen-i-1] & MASKSEQ];
  default:
    printf("Error, call to DNASeq :: operator () with mode >2\n");
    exit(1);
    return 255;
  }
}
// ---------------------------------------------------------------------
//  Print this DNASeq 
// ---------------------------------------------------------------------
void  DNASeq :: Print (FILE * fp)
{
  int i;

  for (i = 0; i < SeqLen; i++) {
    fprintf(fp,"%c", (*this)[i]);
  }
}
// ---------------------------------------------------------------------
// Translate the codon at position i (from the 5' end) in the given sense and 
// return the corresponding amino acid
// ---------------------------------------------------------------------
char DNASeq :: AA(int i, int mode)
{
  int codon,c1,c2,c3;

  if (mode == 1) {
    c3 = (*this)(i-1,1);
    c2 = (*this)(i-2,1);
    c1 = (*this)(i-3,1);
  }
  else{
    c1 = (*this)(i+2,mode);
    c2 = (*this)(i+1,mode);
    c3 = (*this)(i,mode);
  }

  if (Code2NumOne[c1]*Code2NumOne[c2]*Code2NumOne[c3] != 1) 
    return('X');
  else {
    codon = Code2Bit[c1]+Code2Bit[c2]*4+Code2Bit[c3]*16; 
    return(Trad[codon]);
  }
}
// ---------------------------------------------------------------------
// Print a subseq of the DNASeq translated.  Ambiguous codons are
// translated to X. The from-to is not necessarily the start/end of a
// codon (partial gene) and the frame of the gene is given.
// ---------------------------------------------------------------------
void  DNASeq :: PrintTr (FILE* fp, int from, int to, signed char phase)
{
  int loffset,roffset,i,codon,frame,mode = 0;
  int c1,c2,c3;
  int len = -1;
  fprintf(fp,"> %s-%d-%d %d\n",Name,from,to,phase);

  frame = abs(phase) -1;
  from--;
  to--;

  if (from < 0) from = 0;
  if (to >= SeqLen) to = SeqLen-1;

  if (phase < 0) {
    i= SeqLen - 1 - to;
    to = SeqLen - 1 - from;
    from = i;
    mode = 2;
  }
  
  loffset = (3+ frame - from) % 3;
  if (loffset) { //codon partiel a gauche
    from += loffset;
  }
  roffset = (to + 4 - frame) %3;
  if (roffset) { //codon partiel a droite
    to -= roffset;
  }
  
  for (i = from; i < to; i += 3)
    {
      c1 = (*this)(i+2,mode);
      c2 = (*this)(i+1,mode);
      c3 = (*this)(i,mode);
      
      if (Code2NumOne[c1]*Code2NumOne[c2]*Code2NumOne[c3] != 1)
	fprintf(fp,"X");
      else {
	codon = Code2Bit[c1]+Code2Bit[c2]*4+Code2Bit[c3]*16; 
	fprintf(fp,"%c",Trad[codon]);
      }
      len++;
      if ((len%FASTA_Len) == FASTA_Len-1) fprintf(fp,"\n");
    }
  fprintf(fp,"\n");
  return;
}
// ---------------------------------------------------------------------
// Update 0th order Markov model probabilities
// ---------------------------------------------------------------------
void DNASeq :: UpdateMarkov()
{
  int i;
  double Sum;

  Markov0[BitA] = Markov0[BitT] = Markov0[BitG] = Markov0[BitC] = 0.0;

  for (i = 0; i < SeqLen; i++) 
    {
      if (Sequence[i] & CodeA) Markov0[BitA] += 1.0;
      if (Sequence[i] & CodeT) Markov0[BitT] += 1.0;
      if (Sequence[i] & CodeG) Markov0[BitG] += 1.0;
      if (Sequence[i] & CodeC) Markov0[BitC] += 1.0;
    }
  Sum = Markov0[BitA] + Markov0[BitT] + Markov0[BitG] + Markov0[BitC];
  
  Markov0[BitA] /= Sum;
  Markov0[BitT] /= Sum;
  Markov0[BitG] /= Sum;
  Markov0[BitC] /= Sum;

  return;
}
// ---------------------------------------------------------------------
// Compute 0th order Markov model probabilities on nuc i
// ---------------------------------------------------------------------
double DNASeq :: Markov(int i)
{

  if (i >= SeqLen) return 0.0;
  
  return Markov0[Code2Bit[Sequence[i] & MASKSEQ]];
}
// ---------------------------------------------------------------------
// Compute 0th order Markov model probabilities on nuc i
// ---------------------------------------------------------------------
double DNASeq :: MarkovR(int i)
{

  if (i >= SeqLen) return 0.0;
  
  return Markov0[Code2Bit[Code2CCode[Sequence[i] & MASKSEQ]]];
}
// ---------------------------------------------------------------------
// Compute GC/AT probabilities on nuc i
// ---------------------------------------------------------------------
double DNASeq :: GC_AT(int i)
{

  if (i >= SeqLen) return 0.0;
  
  if (Code2UCode[Sequence[i] & MASKSEQ] & (CodeA|CodeT))
    return (Markov0[BitA]+Markov0[BitT])/2.0;
  else
    return (Markov0[BitG]+Markov0[BitC])/2.0;
}
// ---------------------------------------------------------------------
// transforms the char for a nucleotide to its corresponding code
// ---------------------------------------------------------------------

unsigned short DNASeq :: Nuc2Code(char Ch)
{
  Ch = tolower (Ch);
  switch  (Ch)
    {
    case  'a' :
      return CodeA;
      
    case  'c' :
      return CodeC;
      
    case  'g' :
      return CodeG;
      
    case  't' :
      return CodeT;
      
    case  's' :   // c ou g
      return CodeC | CodeG;
      
    case  'w' :   // a ou t
      return CodeA | CodeT;
      
    case  'r' :   // a ou g
      return CodeA | CodeG;
      
    case  'y' :   // c ou t
      return CodeC | CodeT;
      
    case  'm' :   // a ou c
      return CodeA | CodeC;
      
    case  'k' :   // g ou t
      return CodeG | CodeT;
      
    case  'b' :   // c, g ou t
      return CodeC | CodeG | CodeT;
      
    case  'd' :   // a, g ou t
      return CodeA | CodeG | CodeT;
      
    case  'h' :   // a, c ou t
      return CodeA | CodeC | CodeT;
      
    case  'v' :   // a, c ou g
      return CodeA | CodeC | CodeT;
      
    case  'n' :
      return CodeA | CodeC | CodeG | CodeT;
      
    default :
      fprintf (stderr, "Not a DNA character `%c\' in string %s\n", Ch, Name);
      fprintf (stderr, "Replaced by `N'.\n");
      return CodeA | CodeC | CodeG | CodeT;
    }
}


// ---------------------------------------------------------------------
// Compute the code for the unambiguous version at position pos
// ---------------------------------------------------------------------
unsigned short  DNASeq :: Unambit (int pos)
{
  return Code2Bit[Sequence[pos]&MASKSEQ];
}

// ---------------------------------------------------------------------
// Compute the code for the complement 
// ---------------------------------------------------------------------
inline unsigned short Complement (unsigned short code)
{
  return Code2CCode[code & MASKSEQ];
}
  /* version calculatoire
{
  unsigned short res = (code & 1);
  int i;
  
  for (i = 0 ;i < 3 ; i++) {
    res <<= 1;
    code >>= 1;
    res |= (code & 1);
  }

  return res;
}
*/
// ---------------------------------------------------------------------
// Transfer |Len| characters from Pos...  to *To and add null
// terminator.  If Len > 0 go in forward direction; otherwise, go in
// reverse direction and use complements.  The mode indicates the 
// type of information used: ambiguous or unambiguous
// ---------------------------------------------------------------------
void DNASeq :: Transfer(int Pos, int Len, char *To, int mode)
{
  long int  i;
  const char * decode;

  decode = (mode ? Code2UNuc : Code2Nuc);
  
  if  (Len > 0) {
    // We go forward
    for  (i = 0;  i < Len;  i++) 
      To[i] = decode[Sequence[Pos + i] & MASKSEQ];
   
    To [i] = '\0';
  }
  else  {
    for  (i = 0;  i < -Len;  i++) 
      To[i] = decode[Complement(Sequence[Pos - i] & MASKSEQ)];
    
    To [i] = '\0';
  }
  
  return;
}


/* (bl0b) define IsCodeN predicate */
#define CodeN (CodeA|CodeC|CodeG|CodeT)
#define _IsN(_x) (((_x) & CodeN) == CodeN)

// ---------------------------------------------------------------------
// test if an acceptor consensus site starts at position i. Returns true if
// an AG occurs precisely
// ---------------------------------------------------------------------
double DNASeq :: IsAcc(int i,int sens)
{
  int mode = 0;

  if (sens < 0) {
    i = SeqLen -i -1;
    mode = 2;
  }
  
  /* (bl0b) Recognize degenerated sites */
  /*if ((*this)(i,mode) != CodeA) return false;*/
  /*return ((*this)(i+1,mode) == CodeG);*/
  /* bool return version */
  //if ( !( (*this)(i,mode) & CodeA) ) return false;
  //return !!( (*this)(i+1,mode) & CodeG );


	int count = Code2NumOne[(*this)(i,mode) & CodeA] * Code2NumOne[(*this)(i+1,mode) & (CodeG)];
	double degen = (double) Degeneracy(i,mode,2);
	/*fprintf(stderr, "[DEBUG] IsAcc : count=%i degen=%lf\n", count, degen);*/
	return ((double)count) / degen;

#if 0  
  if( ((*this)(i,mode) & CodeA) && ((*this)(i+1,mode)&CodeG) ) {
  	/*int num_N = _IsN((*this)(i,mode)) + _IsN((*this)(i+1,mode));*/
	/*if(num_N) {*/
		/*int N_factor = 1 << (num_N<<1);*/
		/* N_factor is 4 or 16 */
		/*return 1.0 / N_factor;*/
	/*}*/
	/*return 1.0;*/
	double degen = (double) Degeneracy(i,mode,2);
	fprintf(stderr, "[DEBUG] IsAcc : degen=%lf\n", degen);
	return 1.0 / degen;
  }
  return 0.0;
#endif
}
// ---------------------------------------------------------------------
// test if a donor consensus site starts at position i. Returns true if
// a GT/GC may occur
// ---------------------------------------------------------------------
double DNASeq :: IsDon(int i,int sens)
{
  int mode = 0;

  if (sens < 0) {
    i = SeqLen -i -1;
    mode = 2;
  }
  
  /* (bl0b) Recognize degenerated sites */
  
  /*if ((*this)(i,mode) != CodeG) return false;*/
  /*return (((*this)(i+1,mode) == CodeT) || ((*this)(i+1,mode) == CodeC));*/
  /* bool return version */
  //return !!( (*this)(i,mode) & (CodeC|CodeT) );
  /*if( ((*this)(i,mode) & CodeG) && ((*this)(i+1,mode) & (CodeC|CodeT)) ) {*/
  	/*int num_N = _IsN((*this)(i,mode)) + _IsN((*this)(i+1,mode));*/
	/*if(num_N) {*/
		/*int N_factor = 1 << (num_N<<1);*/
		/* N_factor is 4 or 16 */
		/*return 1.0 / N_factor;*/
	/*}*/
	/*return 1.0;*/
  /*}*/
  /*return 0.0;*/
  
  /*if ( !((*this)(i,mode) & CodeG) ) return 0.0;*/

  int count = Code2NumOne[(*this)(i,mode) & CodeG] * Code2NumOne[(*this)(i+1,mode) & (CodeT|CodeC)];
  double degen = (double) Degeneracy(i,mode,2);
  /*fprintf(stderr, "[DEBUG] IsDon : count=%i degen=%lf\n", count, degen);*/
  return ((double)count)/degen;
}

// ---------------------------------------------------------------------
// test if a Stop Codon starts at position i. Returns the number of
// possible matches with a STOP codon i.e. T{AA,AG,GA}
// ---------------------------------------------------------------------
double DNASeq :: IsStop(int i,int sens)
{
  int count = 0;
  int mode = 0;

  if (sens < 0) {
    i = SeqLen -i -1;
    mode = 2;
  }
  
  if (((*this)(i,mode) & CodeT) == 0) return 0.0;

  // le cas des TA{A,G}
  if (((*this)(i+1,mode)) & CodeA)
    count += Code2NumOne[(*this)(i+2,mode) & (CodeA|CodeG)]; 
  
  // le cas du TGA
  if (((*this)(i+1,mode) & CodeG) && ((*this)(i+2,mode) & CodeA))
    count ++;
  
  return (double)count/ Degeneracy(i,mode,3);
}
// ---------------------------------------------------------------------
// test if a spliced Stop Codon may start just before position i. 
// Sets specific bits depending on the occurrence:
//   - bit values 1,2,4 are used to report if the current pos is just 
//     after (FORWARD strand)
//     - 1: a T  
//     - 2: a TG 
//     - 4: a TA 
//   - bit values 8,16,32 are used to report if the seq at hand can generate 
//     a STOP if it appears (REVERSE)
//     - 8 : after G (AT)
//     - 16: after A (GT|AT)
//     - 32: after AG|AA|GA (T) 
// ---------------------------------------------------------------------
int DNASeq :: IsStartStop(int i)
{
  int stop = 0;

  // le T
  if (((*this)(i-1,0)) & CodeT) stop |= isTf; //StopAfterT;

  // le cas du TG/TA
  if (((*this)(i-2,0)) & CodeT) {
    if (((*this)(i-1,0)) & CodeG) stop |= isTGf; //StopAfterTG
    if (((*this)(i-1,0)) & CodeA) stop |= isTAf; //StopAfterTA
  }

  i = SeqLen -i -1;

  // le T
  if (((*this)(i,2)) & CodeT) stop |= isTr; //StopAfterAG

  // le cas du GT/AT
  if (((*this)(i-1,2)) & CodeT) {
    if (((*this)(i,2)) & CodeG) stop |= isTGr; //StopAfterA
    if (((*this)(i,2)) & CodeA) stop |= isTAr; // StopAfterG+StopAfterA
  }

  return stop;
}
// ---------------------------------------------------------------------
// test if a spliced Stop Codon may stop just after position i. 
// Sets specific bits depending on the occurrence.
//   - bit values 1/2/4  are used to report if the seq at hand can generate 
//     a STOP if it appears (FORWARD):
//     - 1:  after T (GA|AG|AA)
//     - 2:  after TG (A)
//     - 4:  after TA (G|A)
//   - bit values 8/16/32 
// REVERSE (RC) are used to report if the current pos is just (C BACKWARD):
//     - 8:  after G (AT)
//     - 16: after A (GT|AT)
//     - 32: after AG|AA|GA (T) 
// ---------------------------------------------------------------------
int DNASeq :: IsStopStop(int i)
{
  int stop = 0;

  if (((*this)(i,0)) & CodeG) {
    stop |= isGf; //StopAfterTA
    if (((*this)(i+1,0)) & CodeA)
      stop |= isGAf; //StopAfterT
  }

  if (((*this)(i,0)) & CodeA) {
    stop |= isAf; //StopAfterTA+StopAfterTG
    if (((*this)(i+1,0)) & (CodeA|CodeG))
      stop |= isARf; //StopAfterT
  }

  i = SeqLen -i -1;

  if (((*this)(i+1,2)) & CodeG) {
    stop |= isGr; //StopAfterG
    if (((*this)(i+2,2)) & CodeA)
      stop |= isGAr;// StopAfterAG
  }
  
  if (((*this)(i+1,2)) & CodeA) {
    stop |= isAr; //StopAfterA
    if (((*this)(i+2,2)) & (CodeA|CodeG))
      stop |= isARr; //StopAfterAG
  }

  return stop;
}
// ---------------------------------------------------------------------
// Returns simply A,T,C,G or N depending on the nucleotide at position i
// ---------------------------------------------------------------------
char DNASeq :: nt(int i, int mode)
{
  if ((*this)(i,mode) & CodeT) return 'T';
  if ((*this)(i,mode) & CodeA) return 'A';
  if ((*this)(i,mode) & CodeG) return 'G';
  if ((*this)(i,mode) & CodeC) return 'C';
  return 'N';
}
// ---------------------------------------------------------------------
// returns a penalty for infrequent start depending on the code of the
// first nuc. (summed up in case of degeneracy ).
// ---------------------------------------------------------------------
double StartTypePenalty(unsigned short code)
{
  double pen = 0.0;

  if (code & CodeA) pen += 0.8;
  if (code & CodeG) pen += 0.1;
  if (code & CodeT) pen += 0.1;
  if (code & CodeC) pen += 0.00001;

  return pen;
}
// ---------------------------------------------------------------------
// test is a Start Codon appears at position i ({ATG}TG) returns the
// probability that the seq is a start codon considering a iid
// model with uniform distribution (1/4) and taking into account 
// start frequencies
// ---------------------------------------------------------------------
double DNASeq :: IsStart(int i,int sens)
{
  unsigned short First;
  int mode = 0;

  if (sens < 0) 
    {
      i = SeqLen -i -1;
      mode = 2;
    }

  if (((*this)(i+2,mode) & CodeG) == 0) return 0.0;
  if (((*this)(i+1,mode) & CodeT) == 0) return 0.0;

  First = (*this)(i,mode);

  if ((First & (CodeA | CodeT | CodeG)) == 0) return 0.0;

  return StartTypePenalty(First) / Degeneracy(i,mode,3);
}

// ---------------------------------------------------------------------
// Version EUCA
// ---------------------------------------------------------------------

double DNASeq :: IsEStart(int i,int sens)
{
  int mode = 0;
  
  if (sens < 0)  {
    i = SeqLen -i -1;
    mode = 2;
  }
  
  if (((*this)(i+2,mode) & CodeG) == 0) return 0.0;
  if (((*this)(i+1,mode) & CodeT) == 0) return 0.0;
  if (((*this)(i,mode) & CodeA) == 0) return 0.0;

  return 1.0 / Degeneracy(i,mode,3);
}
// ---------------------------------------------------------------------
// Degeneracy : returns the number of possible completely known sequence
// represented by a degenerated subsequence
// ---------------------------------------------------------------------
unsigned char DNASeq :: Degeneracy(int i, int mode, int len)
{
  unsigned char ret=1;
  for(int k=0;k<len;k++) {
  	ret *= Code2NumOne[(*this)(i+k, mode)];
  }
  return ret;
}

// ---------------------------------------------------------------------
// Returns Blast Frame (1,2,3 -1, -2, -3)
// ---------------------------------------------------------------------
int DNASeq :: Pos2Frame(int pos, char strand)
{
  int frame=0;
  if (strand == '+') 
  { 
      frame = (pos )% 3; 
      if ( frame == 0 )
	frame = 3 ;
  }
  else
  { 
      frame = (SeqLen - pos -1 ) % 3;
      if ( frame == 0 )
	frame = 3 ;
      frame = frame * ( -1 ) ;
  }
  return frame;
}
