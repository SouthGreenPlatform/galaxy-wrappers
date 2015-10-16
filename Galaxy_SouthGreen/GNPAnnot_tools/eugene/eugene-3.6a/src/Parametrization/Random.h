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
// $Id: Random.h,v 1.4 2008-10-27 09:27:00 tschiex Exp $
// ------------------------------------------------------------------
// File:     Random.h
// Contents: Class for generators of numbers
// ------------------------------------------------------------------

#ifndef RANDOM_H_INCLUDED
#define RANDOM_H_INCLUDED


#include <string>
#include <climits>
#include <cstdlib>
#include <cstring>

class Random {

 public:
  std::string Seed;

  Random(std::string seed);
  Random(void);
  int RandBoundedInt (int bound);
  double RandUniform(void);
  double RandGauss(void);
  int Flip(double prob);

 private:
  bool IsInitialized;
  unsigned short int xsubi[3];
  int iset;
  double gset;

};

#endif
