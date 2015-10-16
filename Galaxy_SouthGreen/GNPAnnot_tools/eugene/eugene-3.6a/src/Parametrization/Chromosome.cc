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
// $Id: Chromosome.cc,v 1.8 2004-09-16 13:06:35 cros Exp $
// ------------------------------------------------------------------
// File:     Chromosome.cc
// Contents: Class Data and Chromosome of the Genetic algorithm
// ------------------------------------------------------------------



#include <stdio.h>
#include <iostream>
#include <math.h>

#include "Chromosome.h"

#include "ParaOptimization.h"
#include "Genetic.h"

extern ParaOptimization OPTIM;


//-------------------------------------------------------
//-------------------------------------------------------
Data::Data(int n)
{
  Genetic* algo = (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];

  for (int i=0; i<n; i++)
    P.push_back( algo->ParMin[i] + 
		 algo->Rand->RandUniform()*(algo->ParMax[i] - algo->ParMin[i]) );
}




//-------------------------------------------------------
// Data distance for clustering
//-------------------------------------------------------
double Data::GetDistance(Data* d)
{
  Genetic* algo = (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];

  double dist=0.0;

  for (unsigned int i=0; i< algo->ParName.size(); i++)
    dist += ( (algo->ParMax[i]==algo->ParMin[i])
	      ? 0 
	      : pow( (d->P[i] - P[i])/(algo->ParMax[i] - algo->ParMin[i]), 2)  );
 
  return dist;
}



//-------------------------------------------------------
//-------------------------------------------------------
void Data::Recenter(void)
{
  Genetic* algo = (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];

  for (unsigned int i=0; i<P.size(); i++) {
    if (P[i] < algo->ParMin[i]) P[i] = algo->ParMin[i];
    if (P[i] > algo->ParMax[i]) P[i] = algo->ParMax[i];
  }
}


//-------------------------------------------------------
//-------------------------------------------------------
void Data::PutBarycenter(int nb1, Data *data2, int nb2) 
{
  if((nb1+nb2)==0){
    nb1=1;nb2=1;
  }

  for (unsigned int i=0; i<P.size(); i++) 
    P[i] = (P[i]*nb1 + data2->P[i]*nb2)/(nb1+nb2);

  Recenter();
}


//-------------------------------------------------------
// Constructor
//-------------------------------------------------------
Chromosome::Chromosome(int n)
{
  RawFitness = 0;
  ScaledFitness = 0;
  IsEvaluated = false;
  IsNew = true;
  IsBestInCluster = false;
  NoCluster = 0;
  data = new Data(n);
}

//-------------------------------------------------------
// Destructor
//-------------------------------------------------------
Chromosome::~Chromosome(void)
{
  delete data;
}


//-------------------------------------------------------
// Evaluate chromosome
//-------------------------------------------------------
void Chromosome::Evaluate(void)
{
  Genetic* algo = (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];
  
  algo->Par = data->P;
  for (unsigned int k=0; k<algo->Para.size(); k++) algo->Para[k] = algo->Par[algo->ParaPar[k]];
  RawFitness = OPTIM.ParaEvaluate(false);

  IsEvaluated = true;
}

//-------------------------------------------------------
// Print values of parameters and RawFitness
//-------------------------------------------------------
void Chromosome::Print(std::string msg)
{
  for (unsigned int i=0; i<data->P.size(); i++) 
    std::cout << data->P[i] << "\t"; 
  std::cout << "Fitness="<<RawFitness <<"\t"<<msg;
  fflush(stdout); fflush(stderr);
}

//-------------------------------------------------------
// Constructor
//-------------------------------------------------------
Cluster::Cluster(int n)
{
  CentralData = new Data(n);
}


//-------------------------------------------------------
// Destructor
//-------------------------------------------------------
Cluster::~Cluster(void)
{
  delete CentralData;
}
