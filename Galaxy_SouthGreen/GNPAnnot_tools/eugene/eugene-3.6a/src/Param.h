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
// $Id: Param.h,v 1.30 2008-07-02 10:05:44 tschiex Exp $
// ------------------------------------------------------------------
// File:     Param.h
// Contents: class Parameters
// ------------------------------------------------------------------

#ifndef  PARAM_H_INCLUDED
#define  PARAM_H_INCLUDED

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <map>
#include <vector>
#include <string>
#include <string.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif
// MacOS-X has getopt() defined is stdlib and the library in the libSystem
#ifndef __APPLE__
#ifdef HAVE_GETOPT_H
#include <getopt.h>
#else
#ifndef HAVE_GETOPT
#include "getopt.h"
#endif
#endif
#endif

#include "Const.h"
#include "System.h"

extern char   *optarg;   
extern int     optind;


class ltstr
{
  public :
    bool operator() (const char* s1, const char* s2) const
    { return (strcmp(s1, s2) <0); };
};

class Parameters
{
 private: 
  std::map <const char*, const char*, ltstr> m;
  std::map <const char*, const char*, ltstr>::iterator iter;

  void UpdateParametersFileName (int argc, char* argv[]);
  void ReadArg   (int, char *[]);
  void ReadPar   (std::string para_file_name);
  void ShowUsage (void);
  
 public:
  FILE *fp;
    
  ~Parameters ();
  void   initParam (int, char *[]);
  int    count(char *key);
  bool   probeKey(char *key, int index = 0);
  char*  getC    (const char *key, int index = 0, int sloppy = 0);
  double getD    (const char *key, int index = 0, int sloppy = 0);
  int    getI    (const char *key, int index = 0, int sloppy = 0);
  int    getUseSensor (char **, int*);
  void   set  (const char *key, const char *value);
  void   setD (const char *key, double n);
  void   setI (const char *key, int n);
  std::string WriteParam (std::vector<std::string> para_name, 
			  std::vector<double> para_val);

  void ResetIter(void);
};

#endif
