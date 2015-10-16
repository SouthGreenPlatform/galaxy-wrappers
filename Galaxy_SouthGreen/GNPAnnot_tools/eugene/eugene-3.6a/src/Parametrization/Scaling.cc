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
// $Id: Scaling.cc,v 1.5 2004-09-16 13:06:35 cros Exp $
// ------------------------------------------------------------------
// File:     Scaling.cc
// Contents: Class realizing the scaling for the genetic algorithm
// ------------------------------------------------------------------

#include <iostream>
#include <math.h>

#include "Scaling.h"

#include "ParaOptimization.h"
#include "Genetic.h"

extern ParaOptimization OPTIM;

//-------------------------------------------------------
// Constructor
//-------------------------------------------------------
Scaling::Scaling(int type)
{
  Type = type;
}


//-------------------------------------------------------
//-------------------------------------------------------
void Scaling::Scalestat (int gen)
{
  Genetic* algo =  (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];

  int i;
  double moy=0.0, moy2=0.0, mini=DBL_MAX, maxi= -DBL_MAX, sigma, f;
  double s=0.1, s0=0.1, p1=0.05, p2=0.1, alpha=0.1;
  double t1, t2, t3, t4, k;
  double C = 2.0;
  
  algo->BestChromosome = NULL;

  /* First, we compute the maximum fitness (maxi), the average fitness moy
     and the quadratic average fitness moy2, and then the standard deviation
     sigma. */
  for (i=0; i<algo->NbElement; i++) {
    if (algo->Population[i]->RawFitness > maxi) {
      maxi = algo->Population[i]->RawFitness;
      algo->BestChromosome = algo->Population[i];
    }
    if (algo->Population[i]->RawFitness < mini) 
      mini = algo->Population[i]->RawFitness;
    moy += algo->Population[i]->RawFitness;
    moy2 += algo->Population[i]->RawFitness*algo->Population[i]->RawFitness;
  }  
  moy /= algo->NbElement;
  algo->AvgFitness = moy;
  moy2 /= algo->NbElement;
  sigma = moy2 - moy*moy;
  if (sigma<0.0) sigma=0.0;
  else sigma=sqrt(sigma);
  algo->SigmaFitness = sigma;

  /*******************************************************/
  /*                         SCALING                     */
  /*******************************************************/
  switch (Type) {
  case 0 : /* no scaling, scaled fitness is equal to raw fitness */
    for (i=0; i< algo->NbElement; i++)
	  algo->Population[i]->ScaledFitness = algo->Population[i]->RawFitness;
    break;

  case 1 : /* Sigma Truncation scaling */
    for (i=0; i<algo->NbElement; i++) {
      f=algo->Population[i]->RawFitness;
      f=f-(moy-C*sigma);
      if (f<0.0000001) f=0.0000001;
      algo->Population[i]->ScaledFitness=f;
    }
    break;

  case 2 : /* Power Law scaling */
    t1=p2*pow((s0/s),alpha);
    t2=((double)gen)/(double)(algo->NbGeneration)*(M_PI/2.0);
    t3=pow((tan(t2)),t1);
    t4=pow((s/s0),p1);
    k=t4*t3;
    for (i=0;i<algo->NbElement;i++) {
      f=algo->Population[i]->RawFitness;
      if (f<0.0)
	{std::cerr <<"ERROR: Negative fitness in Scaling::Scalestat."<<std::endl; exit(100);}
      f=pow(f,k);
      algo->Population[i]->ScaledFitness=f;
    } 
    break;
    
  default :
    {std::cerr <<"ERROR: Invalid type of scaling."<<std::endl; exit(100);}
  } 
}
