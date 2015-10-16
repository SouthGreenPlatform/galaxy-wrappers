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
// $Id: CrossingOver.cc,v 1.6 2004-09-16 13:06:35 cros Exp $
// ------------------------------------------------------------------
// File:     CrossingOver.cc
// Contents: Class for cross over of the Genetic algorithm
// ------------------------------------------------------------------


#include "CrossingOver.h"

#include "ParaOptimization.h"
#include "Genetic.h"

extern ParaOptimization OPTIM;


//-------------------------------------------------------
// Constructor
//-------------------------------------------------------
CrossingOver::CrossingOver(double proba)
{
  Genetic* algo =  (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];
  Chromosome* c;

  Proba = proba;

  if (algo->IsSACrossingOver) 
      for (int i=0; i<2; i++) {
	c = new Chromosome(algo->ParName.size());
	ChildrenPopulation.push_back(c);
      }
}


//-------------------------------------------------------
// Destructor
//-------------------------------------------------------
CrossingOver::~CrossingOver(void)
{
  for (unsigned int i=0; i<ChildrenPopulation.size(); i++) 
      delete ChildrenPopulation[i];
}


//-------------------------------------------------------
//-------------------------------------------------------
void CrossingOver::Crosseval (void)
{
  Genetic* algo = (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];

  /* Number of elements to crossover */
  int max_elems = (int)((double) algo->NbElement * Proba) ;
  int i, n1, n2, k;

  if (!algo->IsSACrossingOver)
    for (i=0; i<max_elems; i+=2) {
      /* step 2, because we change 2 elements each time */
      ChooseCouple(&n1, &n2);
      /* Children will override parents so no need to copy them */
      Crossover(algo->Population[n1], algo->Population[n2]) ;
      /* No need to evaluate parents yet, all evaluations will */
      /* be done later */
      algo->Population[n1]->IsNew = true; algo->Population[n1]->IsEvaluated = false;
      algo->Population[n2]->IsNew = true; algo->Population[n2]->IsEvaluated = false;
    } 
  else
    for (i=0; i<max_elems; i+=2){ 
      ChooseCouple(&n1,&n2);
      ChildrenPopulation[0]->data->P = algo->Population[n1]->data->P;
      ChildrenPopulation[1]->data->P = algo->Population[n2]->data->P;
      Crossover(ChildrenPopulation[0], ChildrenPopulation[1]);
      algo->Par = ChildrenPopulation[0]->data->P;
      for (k=0; k<(int)algo->Para.size(); k++) algo->Para[k] = algo->Par[algo->ParaPar[k]];
      ChildrenPopulation[0]->RawFitness = OPTIM.ParaEvaluate(false);
      algo->Par = ChildrenPopulation[1]->data->P;
      for (k=0; k<(int)algo->Para.size(); k++) algo->Para[k] = algo->Par[algo->ParaPar[k]];
      ChildrenPopulation[1]->RawFitness = OPTIM.ParaEvaluate(false);
      /* Here everything as been calculated, we just have to crossover */
      algo->SA->SATournament4(ChildrenPopulation[0], ChildrenPopulation[1],
			      algo->Population[n1], algo->Population[n2]);
      algo->Population[n1]->IsEvaluated = true;
      algo->Population[n2]->IsEvaluated = true;
    }
}



//-------------------------------------------------------
//-------------------------------------------------------
void CrossingOver::Crossover(Chromosome* c1, Chromosome* c2)
{
  Genetic* algo = (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];

  double alpha;

  Chromosome tmpc1(algo->ParName.size());
  Chromosome tmpc2(algo->ParName.size());

  for (unsigned int i=0; i<c1->data->P.size(); i++) {
    alpha = algo->Rand->RandUniform() * 2.0 - 0.5;
    tmpc1.data->P[i] = (1.0-alpha)*c1->data->P[i] + alpha*c2->data->P[i] ;
    tmpc2.data->P[i] = alpha*c1->data->P[i] + (1.0-alpha)*c2->data->P[i] ;
  }

  tmpc1.data->Recenter();
  tmpc2.data->Recenter();

  c1->data->P = tmpc1.data->P;
  c2->data->P = tmpc2.data->P; 
}

//-------------------------------------------------------
//-------------------------------------------------------
void CrossingOver::ChooseCouple(int* n1, int* n2)
{
  Genetic* algo = (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];

  int k;

  if (algo->IsElitist) { 
    /* If we use an elitist strategy, we can not mutate the overall */
    /* best chromosom (we keep at least one best overall) */
    /* Anyway, we do not mutate an element which has just been added */
    /* in the population */
    k=0;
    do *n1 = (int)( algo->Rand->RandUniform() * algo->NbElement ); 
    while ( ((algo->Population[*n1]->IsNew) ||
	     (algo->Population[*n1]->IsBestInCluster)||
	     (algo->Population[*n1] == algo->BestChromosome)) 
	    && (k++<1000) );
    /* If we went over limit then we must choose one but we */
    /* refuse best_chrom, we accept best_of_cluster to be   */
    /* choosen and then deleted */
    if (k>1000) {
      do *n1 = (int)( algo->Rand->RandUniform() * (double)algo->NbElement ); 
      while ( (algo->Population[*n1] == algo->BestChromosome) || 
	      (algo->Population[*n1]->IsNew) ) ;  // on en sort pas
    } 
    k=0;
    do *n2 = (int)( algo->Rand->RandUniform() * algo->NbElement ); 
    while ( ((algo->Population[*n2]->IsNew) 
	     || (*n2 == *n1)
	     ||(algo->Population[*n2]->IsBestInCluster)
	     ||(algo->Population[*n2] == algo->BestChromosome))
	    && (k++<1000));
    /* If we went over limit then we must choose one but we */
    /* refuse best_chrom, we accept best_of_cluster to be   */
    /* choosen and then deleted */
    if (k>1000) {
      do *n2 = (int)( algo->Rand->RandUniform() * (double)algo->NbElement );
      while ( (algo->Population[*n2] == algo->BestChromosome) || 
	      (algo->Population[*n2]->IsNew) ) ;
    }  
  } else {
    do  *n1 = (int)( algo->Rand->RandUniform() * algo->NbElement );
    while (algo->Population[*n1]->IsNew);
    do  *n2 = (int)( algo->Rand->RandUniform() * algo->NbElement );
    while ( (algo->Population[*n2]->IsNew) || 
	    (*n2 == *n1) );
  } 
}
