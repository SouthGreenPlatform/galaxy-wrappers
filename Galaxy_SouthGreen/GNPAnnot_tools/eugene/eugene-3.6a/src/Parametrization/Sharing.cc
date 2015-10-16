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
// $Id: Sharing.cc,v 1.5 2004-09-16 13:06:35 cros Exp $
// ------------------------------------------------------------------
// File:     Sharing.cc
// Contents: Class realizing the sharing for the genetic algorithm
// ------------------------------------------------------------------

#include <iostream>
#include <math.h>

#include "Sharing.h"

#include "ParaOptimization.h"
#include "Genetic.h"


extern ParaOptimization OPTIM;


//-------------------------------------------------------
// Constructor
//-------------------------------------------------------
Sharing::Sharing(double value)
{
  Value = value;
}



//-------------------------------------------------------
//-------------------------------------------------------
void Sharing::Share (void)
{
  Genetic* algo =  (Genetic*) OPTIM.Algorithms[OPTIM.AlgoIndex];

  int i,j ;
  double t1,t2,t3 ;
  
  /* Use clustering, the complexity of the algorithm is then nlogn */
  if (algo->IsClustering) {
    for (i=0; i<algo->NbElement; i++) {
	if (algo->Clusters[algo->Population[i]->NoCluster]->NbChrom>1) {
	  t1 = algo->Population[i]->data->GetDistance(
	    algo->Clusters[algo->Population[i]->NoCluster]->CentralData ); 
	  if (algo->ClusterDmax() != 0.0)
	    t2 = t1/2.0/algo->ClusterDmax();
	  else 
	    t2 = t1/0.0000001 ;
	  if (t2>=1.0) 
	    std::cerr <<"Warning: bad Dmax in Sharing.cc."<<std::endl; 
	  else {
	    t3 = ((double) algo->Clusters[algo->Population[i]->NoCluster]->NbChrom)
	      *(1-t2*t2*t2);
	    algo->Population[i]->ScaledFitness = algo->Population[i]->ScaledFitness/t3;
	  }
	} 
    }
  } else {
    /* don't use clustering so sharing complexity is n^2 */
    for (i=0; i<algo->NbElement; i++) {
      t1=0.0;
      for (j=0; j<algo->NbElement; j++)
	t1 += EvalShare( 
			algo->Population[i]->data->GetDistance(
					   algo->Population[j]->data) ) ;
      algo->Population[i]->ScaledFitness /= t1;
    }
  }
}


//-------------------------------------------------------
//-------------------------------------------------------
double Sharing::EvalShare(double distance) 
{
  double sigma_share = 1;
  double alpha_share = 1;

  if (distance >= sigma_share)
    return(0) ;
  else
    return(1 - pow(distance/sigma_share,alpha_share)) ;
}
