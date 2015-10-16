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
// $Id: Mutating.cc,v 1.6 2004-09-16 13:06:35 cros Exp $
// ------------------------------------------------------------------
// File:     Mutating.cc
// Contents: Class realizing mutation for the genetic algorithm
// ------------------------------------------------------------------

#include "Mutating.h"

#include "ParaOptimization.h"
#include "Genetic.h"

extern ParaOptimization OPTIM;


//-------------------------------------------------------
// Constructor
//-------------------------------------------------------
Mutating::Mutating(double proba)
{
  Chromosome* c;
  Genetic* algo =  (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];

  Proba = proba;
  if (algo->IsSAMutating) 
      for (int i=0; i<2; i++) {
	c = new Chromosome(algo->ParName.size());
	ChildrenPopulation.push_back(c);
      }
}


//-------------------------------------------------------
// Destructor
//-------------------------------------------------------
Mutating::~Mutating(void)
{
  for (unsigned int i=0; i<ChildrenPopulation.size(); i++) 
    delete ChildrenPopulation[i];
}





//-------------------------------------------------------
//-------------------------------------------------------
void Mutating::Muteval (void)
{
  Genetic* algo =  (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];

  int i,j,n,k;

  n = (int)(algo->NbElement*Proba + 0.5);  

  if (!algo->IsSAMutating)
    for (i=0; i<n; i++) {
      ChooseChromosome(&j);
      Mutation(algo->Population[j]) ;
      algo->Population[j]->IsNew = true;
      algo->Population[j]->IsEvaluated = false;
    }
  else {
    for (i=0; i<n; i++){
      ChooseChromosome(&j) ;
      ChildrenPopulation[0]->data->P = algo->Population[j]->data->P;
      Mutation(ChildrenPopulation[0]) ;
      algo->Par = ChildrenPopulation[0]->data->P;
      for (k=0; k<(int)algo->Para.size(); k++) algo->Para[k] = algo->Par[algo->ParaPar[k]];
      ChildrenPopulation[0]->RawFitness = OPTIM.ParaEvaluate(false);
      /* Here everything as been calculated, we just have to crossover */
      algo->SA->SATournament2(ChildrenPopulation[0], algo->Population[j]);
      algo->Population[j]->IsEvaluated = true;
    }
  }
} 



//-------------------------------------------------------
//-------------------------------------------------------
void Mutating::ChooseChromosome(int* n)
{
  Genetic* algo =  (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];

  int k;

  /* If we use an elitist strategy, we can not mutate the overall */
  /* best chromosom (we keep at least one best overall) */
  /* Anyway, we do not mutate an element which has just been added */
  /* in the population */
  if (algo->IsElitist) {
    k=0;
    do *n = (int)(algo->Rand->RandUniform()*(double)algo->NbElement);
    while ( ((algo->Population[*n] == algo->BestChromosome)
	     ||(algo->Population[*n]->IsNew)
	     ||(algo->Population[*n]->IsBestInCluster) )
	    && (k++<1000) );
    /* If we went over limit then we must choose one but we */
    /* refuse best_chrom, we accept best_of_cluster to be   */
    /* choosen and then deleted */
    if (k>1000) { 
      do *n = (int)(algo->Rand->RandUniform()*(double)algo->NbElement);
      while ( (algo->Population[*n] == algo->BestChromosome)
	      ||(algo->Population[*n]->IsNew) ) ;
    }
  } else
    do *n = (int)(algo->Rand->RandUniform()*algo->NbElement);
    while (algo->Population[*n]->IsNew);
}


//-------------------------------------------------------
//-------------------------------------------------------
void Mutating::Mutation(Chromosome* c)
{
  Genetic* algo =  (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];

  double alpha=0.2;
  // 68% des valeurs mutees tombent dans l'intervalle 
  // [ Valeur initiale-alpha*(Max-Min) -> Valeur initiale+alpha*(Max-Min) ]
  // et 95% entre Val.init-2alpha*(Max-Min) et Val.init+2alpha*(Max-Min)
  // amelioration: diminuer alpha au cours du temps?
  
  /* Evalue la donnée bruitée avec un bruit gaussien*/
  for (unsigned int i=0; i<c->data->P.size(); i++)
    c->data->P[i] = c->data->P[i] + 
      algo->Rand->RandGauss()*alpha*(algo->ParMax[i] - algo->ParMin[i]);

  c->data->Recenter();
}
