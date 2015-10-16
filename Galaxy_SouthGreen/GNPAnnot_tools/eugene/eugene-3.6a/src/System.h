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
// $Id: System.h,v 1.16 2006-11-29 13:34:14 tschiex Exp $
// ------------------------------------------------------------------
// File:     System.h
// Contents: utilitary functions
// ------------------------------------------------------------------

#ifndef SYSTEM_H
#define SYSTEM_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sstream>

#define myassert(x) assert(x)
//#define myassert(x)

#if defined (__sun)
int isinf(double x);
#endif

char * BaseName(char *path);
int ProbeFile(const char *defdir,const char *path);


template<class T>
inline T Min(T x, T y)
{
  if (x < y) return x;
  else
    return y;
}

template<class T>
inline T Max(T x, T y)
{
  if (x < y) return y;
  else
    return x;
}

FILE *FileOpen    (const char *defdir, const char *filename,
		   const char *mode, int sloppy=0);
void *Safe_malloc (size_t Len);
void *Safe_realloc(void * Q, size_t Len);
double cpuTime();
void GetStrDate (char* d);

/**
 *	\fn	to_string
 *	\brief	Convertit n'importe quoi en string
 *	\date	07 novembre 2006
 *	\author	Aurelien.Regat-Barrel ?
 *	Cette fonction a ete copier depuis le site
 *	http://c.developpez.com/faq/cpp/?page=strings#STRINGS_convertform
 *	"L'utilisation de ostringstream permet de convertir en une string n'importe
 *	quel objet pour lequel l'opérateur ostream::operator <<() a été défini."
**/
template<typename T>
std::string to_string( const T & Value )
{
    // utiliser un flux de sortie pour créer la chaîne
    std::ostringstream oss;
    // écrire la valeur dans le flux
    oss << Value;
    // renvoyer une string
    return oss.str();
}

#endif

