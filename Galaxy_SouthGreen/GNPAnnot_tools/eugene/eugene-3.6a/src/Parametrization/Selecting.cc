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
// $Id: Selecting.cc,v 1.4 2004-09-16 13:06:35 cros Exp $
// ------------------------------------------------------------------
// File:     Selecting.cc
// Contents: Class realizing the selection for the genetic algorithm
// ------------------------------------------------------------------

#include <iostream>
#include <vector>

#include "Selecting.h"

#include "ParaOptimization.h"
#include "Genetic.h"
#include "Chromosome.h"

extern ParaOptimization OPTIM;


//-------------------------------------------------------
// Constructor
//-------------------------------------------------------
Selecting::Selecting(int type)
{
  Type = type;

  tab = NULL;
  fraction = NULL;
  choices = NULL;
}


//-------------------------------------------------------
// Destructor
//-------------------------------------------------------
Selecting::~Selecting(void)
{
  if (tab!=NULL)
    free(tab);
  else {
    free(fraction);
    free(choices);
  }
}


//-------------------------------------------------------
//-------------------------------------------------------
void Selecting::Selection(void)
{
  Genetic* algo =  (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];

  int i, j, k, l, p, n, jassign, jpick;
  double a, expected, avg_s_fitness, maxo;
  double total_fitness=0.0;
  std::vector<Chromosome*> buf_pop;
 
  if (Type==0) {
    /*************************************************************/
    /*       Roulette Wheel Selection                            */
    /*************************************************************/
    if (tab == NULL)
      tab = (double*)malloc(algo->NbElement*sizeof(double));
    
    for (i=0; i<algo->NbElement; i++)
      total_fitness += algo->Population[i]->ScaledFitness;
    
    tab[0] = algo->Population[0]->ScaledFitness/total_fitness;
    for (i=1; i<algo->NbElement; i++)
      tab[i]=algo->Population[i]->ScaledFitness/total_fitness+tab[i-1];
    
    for (i=0; i<algo->NbElement; i++) {
      a=algo->Rand->RandUniform();
      n=algo->NbElement/2;
      p=(n+1)/2;
      while(true) {
	if (a<tab[n]) {
	  if ((n==0) || (a> tab[n-1])) break;
	  n=n-p;
	  if (n<0) n=0;
	} else {
	  if ((n==algo->NbElement-1)||(a<tab[n+1])) break;
	  n=n+p;
	  if (n>=algo->NbElement) n=algo->NbElement-1;
	} 
	p=(p+1)/2;
      } /*while*/

      algo->NewPopulation[i]->IsNew = false;
      algo->NewPopulation[i]->data->P = algo->Population[n]->data->P;
      algo->NewPopulation[i]->RawFitness = algo->Population[n]->RawFitness;
      algo->NewPopulation[i]->IsEvaluated = algo->Population[n]->IsEvaluated;
      algo->NewPopulation[i]->IsBestInCluster = false;
    } /*for*/
  } else {
    /*************************************************************/
    /*  Stochastic Remainder without Replacement Selection       */
    /*************************************************************/
    if (choices==NULL) {
      choices = (int*)malloc(algo->NbElement*sizeof(int));
      fraction = (double*)malloc(algo->NbElement*sizeof(double));
    }

    /* calcul de la moyenne de la scaled fitness */
    /*-------------------------------------------*/
    avg_s_fitness=0.0;
    for (i=0; i<algo->NbElement; i++)
      avg_s_fitness = avg_s_fitness + algo->Population[i]->ScaledFitness;
    avg_s_fitness = avg_s_fitness/(double)algo->NbElement;

    if(avg_s_fitness == 0)
      for(j=0; j<algo->NbElement; j++) choices[j] = j;
    else {
      j = 0;
      k = 0;
      
      /* Assign whole numbers */
      do { 
	expected = algo->Population[j]->ScaledFitness /avg_s_fitness;
	jassign = (int) expected; 
	/* note that expected is automatically truncated */
	fraction[j] = expected - jassign;
	while(jassign > 0) {
	  jassign--;
	  if ( (k>=algo->NbElement) || (k<0) )
	    std::cout << k << " bang...." <<std::endl;
	  else
	    choices[k] = j;
	  k++;
	} /*while*/
	j++;
      } while (j < algo->NbElement);
	      
      j = 0;
      /* Assign fractional parts */
      while(k < algo->NbElement) { 
	if(j >= algo->NbElement) j = 0;
	if(fraction[j] > 0.0) {
	  /* A winner if true */
	  if(algo->Rand->Flip(fraction[j])) {
	    choices[k] = j;
	    fraction[j] = fraction[j] - 1.0;
	    k++;
	  } /*if*/
	} /*if*/
	j++;
      } /*while*/
    } /*if*/
    nremain = algo->NbElement - 1;
    
    for (i=0; i<algo->NbElement; i++) {
      jpick =(int)(algo->Rand->RandUniform()*nremain);
      n = choices[jpick];
      choices[jpick] = choices[nremain];
      nremain--;
      
      algo->NewPopulation[i]->IsNew = false;
      algo->NewPopulation[i]->data->P = algo->Population[n]->data->P;
      algo->NewPopulation[i]->RawFitness = algo->Population[n]->RawFitness;
      algo->NewPopulation[i]->IsEvaluated = algo->Population[n]->IsEvaluated ;
      algo->NewPopulation[i]->IsBestInCluster = false;
    } /*for*/
    
  } /*if*/


  /* Petit addendum pour la strategie elitiste */
  /* On elimine au hasard un element de la population et on  
     le remplace par le  meilleur de la generation precedente, pour 
     etre sur de ne pas le perdre. */
  if (algo->IsElitist) {
    i=(int)(algo->Rand->RandUniform()*algo->NbElement);
    algo->NewPopulation[i]->data->P = algo->BestChromosome->data->P;
    algo->NewPopulation[i]->RawFitness = algo->BestChromosome->RawFitness;
    algo->NewPopulation[i]->IsEvaluated = true;
    algo->BestChromosome = algo->NewPopulation[i];
    /* On reutilise le flag IsNew pour ne pas supprimer dans la population les
       meilleurs elements de chaque cluster alors que l'on vient precisement de
       les sauvegarder. C'est pour cette raison qu'il faudra remettre tous
       les elements a FALSE a la fin de la boucle suivant. */
    algo->NewPopulation[i]->IsNew = true;
    /* Lorsque l'on fait du sharing, on etend l'elitisme a chaque cluster.
       L'idee est de conserver le meilleur element de chaque cluster,
       a condition qu'il soit meilleur qu'une certaine valeur de 
       reference, controlable par un parametre (maxo) a partir de la moyenne
       et du sigma de la fitness */
    if (algo->IsSharing) {
      /* On va pour chaque cluster trouver le meilleur element et verifier
	 que sa fitness est superieure a maxo. */
      maxo = algo->BestChromosome->RawFitness * algo->Elitism ;
      for (j=0; j<(int)algo->NbCluster; j++) { 
	/* Get index of best element from cluster j and put it in j */
	/* if this cluster is optimal, else put -1 in it */
	l = algo->GetOptimalInCluster(j, maxo) ;
	if (l!=-1) {
	  k=0;
	  do {
	    i=(int)(algo->Rand->RandUniform() * algo->NbElement);
	  } while ( (algo->NewPopulation[i]->IsNew)&& (k++<1000) ) ;
	  if (k<1000) {
	    algo->NewPopulation[i]->data->P = algo->Population[l]->data->P;
	    algo->NewPopulation[i]->RawFitness = algo->Population[l]->RawFitness;
	    algo->NewPopulation[i]->IsEvaluated = true;
	    /* Ici simple reutilisation du flag new pour ne pas ecraser */
	    /* cet element dans la prochaine iteration de la boucle     */
	    algo->NewPopulation[i]->IsBestInCluster = true;
	    algo->NewPopulation[i]->IsNew = true;
	  }
	}
      }
    }
    
    /* On a utilise les flags new, on remet tout a false */      
    for (i=0; i<algo->NbElement; i++)
      algo->NewPopulation[i]->IsNew = false;
  } /* Fin de l'addendum */

  buf_pop = algo->Population;
  algo->Population = algo->NewPopulation;
  /* Sert a ne pas avoir a re-reserver l'espace memoire la fois suivante */
  /* on sauvegarde ainsi l'adresse d'une zone allouee a la bonne taille  */
  algo->NewPopulation = buf_pop; 
}
