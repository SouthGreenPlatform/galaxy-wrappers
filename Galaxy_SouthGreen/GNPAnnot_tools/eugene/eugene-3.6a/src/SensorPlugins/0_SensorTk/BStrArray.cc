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
// $Id: BStrArray.cc,v 1.3 2004-09-16 13:13:05 cros Exp $
// ------------------------------------------------------------------
// File:     BStrArray.cc
// Contents: Access to binary IMM model
// ------------------------------------------------------------------


#include <math.h>

#include "../../Const.h"
#include "../../System.h"
#include "BStrArray.h"
#include "EndianConv.h"

// ---------------------------------------------------------------------
//  Default constructor.
// ---------------------------------------------------------------------
BString_Array :: BString_Array  ()
{
  Val = NULL;
  Num_Entries = Max_Str_Len = Alphabet_Size = 0;
}

// ---------------------------------------------------------------------     
//  Construct a  BString_Array  that can be indexed by strings of length
//  up to  L  over an alphabet of size  S .
// ---------------------------------------------------------------------
BString_Array :: BString_Array(int L, int S)
{
  int  i;
  
  myassert (L > 0 && S > 0);
  
  Max_Str_Len = L;
  Alphabet_Size = S;
  
  Num_Entries = int (pow (S, L + 1) - 1) / (S - 1);
  Val = (unsigned short *) Safe_malloc (Num_Entries *
					sizeof (unsigned short));
  
  Offset = new int[L+1];
  for (i = 0; i < 10; i++) {
    Offset[i] = int (pow (Alphabet_Size, i) - 1) / (Alphabet_Size - 1);
  }
  
  Val [0] = 65535U;
  for  (i = 1;  i < Num_Entries;  i ++)
    Val [i] = 0;
}

// ---------------------------------------------------------------------
//  Destroy this  BString_Array  by freeing its memory.
// ---------------------------------------------------------------------
BString_Array :: ~ BString_Array  ()
{
  delete [] Offset;
  free (Val);
}


// ---------------------------------------------------------------------
// ead a BString_Array from file. Returns non zero if a problem occurs
// ---------------------------------------------------------------------
int  BString_Array :: Read (FILE * fp)
{
  int  i, A, N;
  unsigned int M;
  char endian;
  int Ok;
  
  endian = 0;
  M = A = N = 0;

  Ok = fread (&M, sizeof(int), 1, fp);
  if (!Ok) return 1;

  if (M == LEndianReverse(Max_Str_Len))
    {
      endian = 1;
      M = Max_Str_Len;
    }
  
  Ok = fread (&A, sizeof(int), 1, fp);
  if (!Ok) return 1;

  if (endian) A = LEndianReverse(A);
  
  Ok = fread (&N, sizeof(int), 1, fp);
  if (!Ok) return 1;

  if (endian) N = LEndianReverse(N);
  
  myassert (M == Max_Str_Len && A == Alphabet_Size && N == Num_Entries);
  
  Ok = fread (Val, sizeof(unsigned short), Num_Entries, fp);
  if (Ok != Num_Entries) return 1;
  
  if (endian)
    for (i = 0;  i < Num_Entries;  i ++)
      Val[i] = SEndianReverse(Val[i]);

  return 0;     
}

// ---------------------------------------------------------------------
//  Convert string  S [0 .. L-1]  to a subscript in the  Val  array.
// ---------------------------------------------------------------------

int  BString_Array :: String_To_Sub  (DNASeq *S, unsigned int P, unsigned int L)
{
  unsigned int  i, Sub;
  
  myassert (L <= Max_Str_Len);
  
  if  (L == 0) return  0;
  
  Sub = 0;
  for  (i = 0;  i < L;  i ++)
    Sub = Sub * Alphabet_Size + S->Unambit(i+P);
  
  Sub += Offset[L];
  
  return  Sub;
}

// ---------------------------------------------------------------------
//  Convert the antiparallel of string S [P .. P+L-1] to a subscript in
//  the Val array.
// ---------------------------------------------------------------------
int  BString_Array :: AntiString_To_Sub  (DNASeq *S, unsigned int P, unsigned int L)
{
  int  i, Sub;
  
  myassert (L <= Max_Str_Len);
  
  if  (L == 0) return  0;
  
  Sub = 0;
  for  (i = L-1;  i >= 0;  i --) 
    Sub = Sub * Alphabet_Size + (3-S->Unambit(i+P));
  
  Sub += Offset[L];
  
  return  Sub;
}

// ---------------------------------------------------------------------
//  Return a reference to the value associated with string  S [0 .. L-1] .
// ---------------------------------------------------------------------
unsigned short &  BString_Array :: operator ()  (DNASeq *S, unsigned int P, unsigned int L)
{
  int  i;
  
  i = String_To_Sub(S, P,L);
  return  Val [i];
}


