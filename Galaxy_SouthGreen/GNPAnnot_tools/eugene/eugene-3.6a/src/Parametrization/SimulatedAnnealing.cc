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
// $Id: SimulatedAnnealing.cc,v 1.3 2004-09-16 13:06:35 cros Exp $
// ------------------------------------------------------------------
// File:     SimulatedAnnealing.cc
// Contents: Class realizing the simulated annealing for the genetic algorithm
// ------------------------------------------------------------------

#include <iostream>
#include <math.h>

#include "SimulatedAnnealing.h"
#include "ParaOptimization.h"
#include "Genetic.h"

extern ParaOptimization OPTIM;

//-------------------------------------------------------
// Constructor
//-------------------------------------------------------
SimulatedAnnealing::SimulatedAnnealing(void)
{
  CP=1;
  CC1=0.8;
  CC2=0.95;
}

//-------------------------------------------------------
//-------------------------------------------------------
void SimulatedAnnealing::SAInit(void)
{
  Genetic* algo = (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];

  double Ts,Tf,DE_MAX,DE_MIN,k,z1,z2;
  double diff;
  int i,j;
  int z;

  /* DE_MAX : expected maximum fitness deviation between two generations */
  /*          at the bigining of the annealing process                   */
  /* DE_MIN : expected minimum fitness deviation between two generations */
  /*          at the end of the annealing process                        */ 
  
  DE_MIN=DBL_MAX;
  DE_MAX=0.0;
  for (i=0;i<algo->NbElement;i++)
    for (j=0;j<algo->NbElement;j++)
      {
	diff=fabs((algo->Population[i])->RawFitness - (algo->Population[j])->RawFitness);
	if ((diff!=0.0)&&(diff<DE_MIN))
	  DE_MIN=diff;
        if (diff>DE_MAX)
          DE_MAX=diff;
      }

  if (DE_MAX==0.0) {
    DE_MAX=0.1;
    DE_MIN=DE_MAX/2.0;
  }

  if (algo->IsTracing) std::cout << "DE_MAX=" << DE_MAX << "\tDE_MIN=" << DE_MIN <<std::endl;
  k=0.75;
  Ts=-DE_MAX/log(1/k-1); /* Starting Temperature */
  k=0.99;
  Tx=-DE_MAX/log(1/k-1); /* Transition Temperature */
  k=0.99;
  Tf=-DE_MIN/log(1/k-1); /* Final Temperature */
      

  z1=(log(Tx/Ts))/log(CC1);/* number of temperature levels during the first */
                           /* re-annealing scheme                           */

  z2=(log(Tf/Tx))/log(CC2);/* number of temperature levels during the second */
                           /* re-annealing scheme                            */
      
  z=(int)(z1+z2);
  CC=CC1;
  /* CP : number of generations in each temperature level */
  /*      (as large as possible)                          */
  
  CP=(int)((double)algo->NbGeneration/z+0.5);
  if (CP==0) CP=1;
  algo->NbGeneration=CP*z;
  T=Ts; 
  
  std::cout << "Ts=" <<Ts<< "\tTx=" <<Tx<< "\tTf=" <<Tf <<std::endl;
  std::cout << "z=" <<z<< "\tNbGeneration=" << algo->NbGeneration << "\tCP=" <<CP <<std::endl;
}


//-------------------------------------------------------
/* re-annealing scheme */
//-------------------------------------------------------
void SimulatedAnnealing::SAUpdate(int gen)
{	  
  if (gen%CP==0)
    {
      T=CC*T;
      if (T<Tx) CC=CC2;
    }
}


//-------------------------------------------------------
/* tournament selection between children and parents*/
/* Child is selected if
   Its fitness is better than the fitness of its parent
   OR its fitness is lower but the differences between parent and 
   child fitness are "not too high" according to simulated annealing */
//-------------------------------------------------------
void SimulatedAnnealing::SATournament4(Chromosome *C1,Chromosome *C2, 
				       Chromosome *P1, Chromosome *P2)
{
  Genetic* algo = (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];

  double r,Ei,Ef;
  double admission_proba,delta;

  Ei=P1->RawFitness;
  Ef=C1->RawFitness;
  delta=(Ei-Ef)/T;
  if (delta<-100.0)
    admission_proba=1.0;
  else
    if (delta>100.0)
      admission_proba=0.0;
    else
      admission_proba=1/(1+exp(delta));
  r=algo->Rand->RandUniform();
  if (r<admission_proba)
    {
      P1->RawFitness=C1->RawFitness;
      P1->IsNew=true;
      P1->data = C1->data;       
    }

  Ei=P2->RawFitness;
  Ef=C2->RawFitness;
  delta=(Ei-Ef)/T;
  if (delta<-100.0)
    admission_proba=1.0;
  else
    if (delta>100.0)
      admission_proba=0.0;
    else
      admission_proba=1/(1+exp(delta));
  r=algo->Rand->RandUniform();
  if (r<admission_proba)
    {
      P2->RawFitness=C2->RawFitness;
      P2->IsNew=true;
      P2->data->P = C2->data->P;       
    }
}


//-------------------------------------------------------
//-------------------------------------------------------
void SimulatedAnnealing::SATournament2(Chromosome* c1, Chromosome* c2)
{
  Genetic* algo = (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];

  double r, Ei, Ef;
  double admission_proba, delta;
  
  /* tournament selection between children and parents*/
  Ei=c2->RawFitness;
  Ef=c1->RawFitness;
  delta=(Ei-Ef)/T;
  if (delta<-100.0)
    admission_proba=1.0;
  else
    if (delta>100.0)
      admission_proba=0.0;
    else
      admission_proba=1/(1+exp(delta));
  r=algo->Rand->RandUniform();
  if (r<admission_proba){
    c2->RawFitness = c1->RawFitness;
    c2->IsNew = true;
    c2->data->P = c1->data->P;
  }
}  

