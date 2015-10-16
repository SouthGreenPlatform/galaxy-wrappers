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
// $Id: Random.cc,v 1.3 2004-09-16 13:06:35 cros Exp $
// ------------------------------------------------------------------
// File:     Random.cc
// Contents: Class for generators of numbers
// ------------------------------------------------------------------

#include <iostream>
#include <unistd.h>
#include <math.h>

#include "Random.h"


//-------------------------------------------------------
// constructor: Initialise the seed of the generator
// if the argument seed is "ALEA" then use time and pid to compute the seed,
// else use the integer defined in seed
//-------------------------------------------------------
Random::Random(std::string seed)
{
  Seed = seed;
  int s = atoi(seed.c_str());
  
  if (seed == "ALEA")
    srandom((unsigned) ((int) time(NULL)) * getpid());
  else
    srandom( atoi(seed.c_str()) );
  
  xsubi[0]=s&0xffff;
  xsubi[1]=(s&0xffff0000)>>4;
  xsubi[2]=0;

  iset=0;
  gset=0;

  IsInitialized = true;  
}


//-------------------------------------------------------
// return an integer between 0 and 1
//-------------------------------------------------------
int Random::RandBoundedInt (int bound)
{
  int n;
  if (IsInitialized)
    n = random()% bound;
  else {
    std::cerr <<"ERROR: Call to Random::RandBoundedInt method without having initialized the generator."
	      <<std::endl; 
    exit(100);
  }

    return n;
}
   

//-------------------------------------------------------
//-------------------------------------------------------
double Random::RandUniform(void)
{
  return erand48(xsubi);
}


//-------------------------------------------------------
//-------------------------------------------------------
double Random::RandGauss(void)
{
  double fac,r,v1,v2;

  if (iset == 0) {
    do {
      v1=2.0*RandUniform()-1.0;
      v2=2.0*RandUniform()-1.0;
      r=v1*v1 + v2*v2;
    } while (r >= 1.0);
    fac=sqrt(-2.0*log(r)/r);
    gset = v1*fac;
    iset = 1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}


//-------------------------------------------------------
// Flip a biased coin - true if heads
//-------------------------------------------------------
int Random::Flip(double prob)
{
  if((RandUniform()) <= prob)
    return(1);
  else
    return(0);
}

