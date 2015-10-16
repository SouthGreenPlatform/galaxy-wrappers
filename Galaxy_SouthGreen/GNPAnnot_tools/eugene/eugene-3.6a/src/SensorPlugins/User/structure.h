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
// $Id: structure.h,v 1.3 2004-09-16 13:36:14 cros Exp $
// ------------------------------------------------------------------
// File:     structure.h
// Contents: definition of structures for the User sensor
// ------------------------------------------------------------------


#ifndef STRUCTURE_H_INCLUDED
#define STRUCTURE_H_INCLUDED

typedef struct UTIL *ptUTIL ;

typedef struct UTIL
{
  int tab;
  int   n1;
  int   n2;
  double  delta;
  char* rais;
  bool check;
  ptUTIL suiv;
} UTIL;

extern int Utilisateur(char *nom_fich,ptUTIL *a, ptUTIL *b);
void Util(int i, ptUTIL ut, DATA *d);
void WriteUtils(ptUTIL,FILE *);

#endif 
