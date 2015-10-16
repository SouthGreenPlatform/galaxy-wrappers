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
// $Id: BStrArray.h,v 1.3 2004-09-16 13:13:05 cros Exp $
// ------------------------------------------------------------------
// File:     BStrArray.h
// Contents: Access to binary IMM model
// ------------------------------------------------------------------

#ifndef  BSTRARRAY_H_INCLUDED
#define  BSTRARRAY_H_INCLUDED

#include "../../DNASeq.h"

class  BString_Array
{
 private:
  unsigned int  Max_Str_Len;
  int  Alphabet_Size;
  int  Num_Entries;
  unsigned short  * Val;
  
 public:
  int *Offset;
  BString_Array  ();
  BString_Array  (int, int);
  ~ BString_Array  ();
  int  Read  (FILE *);    
  int  String_To_Sub  (DNASeq  *, unsigned int, unsigned int);
  int  AntiString_To_Sub  (DNASeq  *, unsigned int, unsigned int);

  unsigned short operator [] (int i) {
    return  Val[i];
  }
  unsigned short &  operator ()  (DNASeq *, unsigned int, unsigned int);
};

#endif
