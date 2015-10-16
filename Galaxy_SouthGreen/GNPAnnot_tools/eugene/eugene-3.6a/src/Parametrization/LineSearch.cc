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
// $Id: LineSearch.cc,v 1.15 2005-03-17 14:22:26 cros Exp $
// ------------------------------------------------------------------
// File:     LineSearch.cc
// Contents: Class of the Line Search algorithm
// ------------------------------------------------------------------


#include <iostream>
#include <string>

#include "LineSearch.h"

#include "../Param.h"
#include "ParaOptimization.h"
#include "Random.h"

extern Parameters PAR;
extern ParaOptimization OPTIM;

//-------------------------------------------------------
// Constructor
//-------------------------------------------------------
LineSearch::LineSearch (void) : OptiAlgorithm()
{
  double step;
  int i;
  std::vector <bool> FirstEltCluster;

  NbMaxCycle = PAR.getI("LineSearch.NbMaxCycle");
  NbMinCycle = PAR.getI("LineSearch.NbMinCycle");
  NbMaxStab = PAR.getI("LineSearch.NbMaxStab");
  DivInter = PAR.getI("LineSearch.DivInter");
  Alpha = PAR.getD("LineSearch.Alpha");
  EvolutionMini = PAR.getD("LineSearch.EvolutionMini");
  
  for (int i=0; i<PAR.getI("ParaOptimization.NbParameter"); i++) {
    Para.push_back( PAR.getD("LineSearch.Para.Init", i) );
    ParaMinStep.push_back( PAR.getD("LineSearch.Para.Step", i) );
    ParaMinInter.push_back( PAR.getD("LineSearch.Para.MinInit", i) );
    ParaMaxInter.push_back( PAR.getD("LineSearch.Para.MaxInit", i) );
    ParaLInter.push_back( ParaMaxInter[i] - ParaMinInter[i] );
  }
    
  for (unsigned int i=0; i< Para.size(); i++) {
    step = (ParaMaxInter[i] - ParaMinInter[i]) / DivInter;
    ParaStep.push_back( Max(step, ParaMinStep[i]) );
  }
  
  Rand = new Random(PAR.getC("LineSearch.Seed"));

  // Update de names of parameters with just the first elt of clusters IDENTICAL
  for (i=0; i<(int)ParaClusters.size(); i++) FirstEltCluster.push_back( true );
  MsgParaNames="";
  for (i=0; i<(int)ParaName.size(); i++) 
    if (ParaClusterRelations[ParaCluster[i]] == IDENTICAL) {
      if ( FirstEltCluster[ParaCluster[i]] ) {
	MsgParaNames += ReduceName(ParaName[i]) + "\t";
	FirstEltCluster[ParaCluster[i]] = false;
      }
    } else
      MsgParaNames += ReduceName(ParaName[i]) + "\t";
}



/*======================================================================*/
/* Optimize: launches optimization                                      */
/*======================================================================*/
void LineSearch :: Optimize(bool is_chaining) 
{
  int n;
  int STOP;
  double FitOptPrec = 0;
  unsigned int i; int c;
  std::string warning_message;
  const double semi_interval_reduction_coef = 0.10; // +/- 10% of previous min-max interval

  // if it chains another algo then takes the optimum found as initial point
  // and adapt intervals
  if (is_chaining) {
    double para_interval;
    Para = OPTIM.Algorithms[OPTIM.AlgoIndex-1]->Para;
    for (i=0; i<Para.size(); i++) {
      para_interval = (ParaMax[i] - ParaMin[i]) * semi_interval_reduction_coef;
      ParaMin[i] = Max(ParaMin[i], Para[i] - para_interval);
      ParaMax[i] = Min(ParaMax[i], Para[i] + para_interval);
      ParaMinInter[i] = ParaMin[i];
      ParaMaxInter[i] = ParaMax[i];
      ParaLInter[i] = ParaMaxInter[i] - ParaMinInter[i];
      ParaStep[i] = Max((ParaMaxInter[i] - ParaMinInter[i]) / DivInter, ParaMinStep[i]);
    }
  }


  // Display parametrization of the algorithm
  std::cout <<std::endl << "---------------------------------------------------------------"<<std::endl;
  std::cout << "Optimization of eugene parameters with the LineSearch algorithm"<<std::endl<<std::endl;
  std::cout << "---------------------------------------------------------------"<<std::endl;
  std::cout << "Parametrisation of the algorithm:"<<std::endl<<std::endl;
  std::cout << "NbMaxCycle: " << NbMaxCycle << "\tNbMinCycle: " << NbMinCycle 
       << "\tNbMaxStab: " << NbMaxStab <<std::endl;
  std::cout << "DivInter: " << DivInter << "\tAlpha: " << Alpha 
       << "\tEvolutionMini: " << EvolutionMini <<std::endl;
  std::cout << "Seed: " << Rand->Seed <<std::endl;
  std::cout << "Trace: " << ((IsTracing) ? "TRUE" : "FALSE") <<std::endl<<std::endl;

  std::cout << "NbParameter: " << (int) Para.size() << "\tNbCluster: " << NbParaCluster <<std::endl;
  std::cout << "Param  : \t"; for (i=0; i<ParaName.size(); i++) std::cout << ReduceName(ParaName[i]) << "\t"; std::cout <<std::endl;
  std::cout << "Init   : \t"; for (i=0; i<Para.size(); i++) std::cout << Para[i] << "\t"; std::cout <<std::endl;
  std::cout << "Min    : \t"; for (i=0; i<ParaMin.size(); i++) std::cout << ParaMin[i] << "\t"; std::cout <<std::endl;
  std::cout << "Max    : \t"; for (i=0; i<ParaMax.size(); i++) std::cout << ParaMax[i] << "\t"; std::cout <<std::endl;
  std::cout << "MinInit: \t"; for (i=0; i<ParaMinInter.size(); i++) std::cout << ParaMinInter[i] << "\t"; std::cout <<std::endl;
  std::cout << "MaxInit: \t"; for (i=0; i<ParaMaxInter.size(); i++) std::cout << ParaMaxInter[i] << "\t"; std::cout <<std::endl;
  std::cout << "MinStep: \t"; for (i=0; i<ParaMinStep.size(); i++) std::cout << ParaMinStep[i] << "\t"; std::cout <<std::endl;

  std::cout <<std::endl;
  for (c=0; c<NbParaCluster; c++) {
    std::cout << "Cluster[" << c << "]: \t";
    for (unsigned int j=0; j<ParaClusters[c].size(); j++) 
      std::cout << ParaName[((ParaClusters[c])[j])] <<" ";
    if (ParaClusterRelations[c] == IDENTICAL) 
      std::cout << "   IDENTICAL"<<std::endl;
    else
      std::cout << "   LINKED"<<std::endl;
  }

  Fitness = OPTIM.ParaEvaluate(false);
  std::cout <<std::endl<< "Fitness of initial point: " << Fitness <<std::endl<<std::endl;
  FitOpt = Fitness; Optimums.clear(); Optimums.push_back(Para); 


  // general loop
  n=0; STOP = 0;
  while ((n < NbMaxCycle) && ((n < NbMinCycle) || (STOP < NbMaxStab))) {
    std::cout <<std::endl <<"Cycle: " << (n + 1)
	 << " -----------------------------------------------------\n";
    if (IsTracing) std::cout << MsgParaNames << std::endl;
    
    for (int k = 0; k < NbParaCluster; k++) { 
      ScanCluster(k);
      ChooseOptimal();

      if (IsTracing) {std::cout <<"Local optimal point:"<<std::endl; PrintParam();std::cout <<std::endl;}
    }
    
    if ( (n!=0) && ((FitOpt - FitOptPrec) < EvolutionMini) )
      STOP += 1; 
    else { 
      STOP = 0;
      FitOptPrec = FitOpt; 
    }
     
    ReduceSearch();
    n++;

    std::cout <<std::endl << "Optimal point of the cycle:" <<std::endl; 
    for (i=0; i<ParaName.size(); i++) std::cout << ReduceName(ParaName[i]) << "\t"; std::cout <<std::endl;
    for (i=0; i<ParaName.size(); i++) std::cout << Para[i] << "\t"; 
    std::cout <<"Fitness="<<FitOpt<<std::endl;
    if (IsTracing) {OPTIM.ParaEvaluate(true); std::cout <<OPTIM.DetailedEvaluation;}
    fflush(stdout); fflush(stderr);
  }
  
  std::cout <<std::endl<< "---------------------------------------------------------------"<<std::endl;
  std::cout << "LineSearch stops after " << n << " cycles."<<std::endl;
  if (n==NbMaxCycle) std::cout <<"Maximum number of cycles achieved."<<std::endl;
  if (STOP==NbMaxStab) std::cout <<"Fitness is stable since "<<NbMaxStab<<" cycles." <<std::endl;
  std::cout <<std::endl << "Final Optimal Point:" <<std::endl; 
  for (i=0; i<ParaName.size(); i++) std::cout << ReduceName(ParaName[i]) << "\t"; std::cout <<std::endl;
  for (i=0; i<ParaName.size(); i++) std::cout << Para[i] << "\t"; 
  std::cout <<"Fitness="<<FitOpt<<std::endl;
  std::cout << "---------------------------------------------------------------"<<std::endl;
  std::cerr <<std::endl<< warning_message ;
}


/*======================================================================*/
/* UpdateOpt: Update list of points which fitness is equal to current   */
/*            optimal fitness                                           */
/*======================================================================*/ 
void LineSearch :: UpdateOpt(void) 
{  
  if (FitOpt <= Fitness) { 
    if (FitOpt < Fitness) {
      FitOpt = Fitness;
      Optimums.clear(); 
    }
    
    Optimums.push_back(Para); 
  }
}


/*======================================================================*/
/* ScanCluster: scans parameters of current cluster                     */
/*======================================================================*/
void LineSearch :: ScanCluster(int k) 
{  
  int n, m;
  unsigned int i;
  double t;

  if (ParaClusterRelations[k] == IDENTICAL) { 
    // We assume that identical parameters do have identical intervals
    n = (ParaClusters[k])[0];
	  
    for (int q=0; q<=DivInter; q++) {
      t = ParaMinInter[n] + q*ParaStep[n];
      // all the parameters of the cluster are set the same value
      for (i=0; i<ParaClusters[k].size(); i++) {
	m = (ParaClusters[k])[i];
	Para[m] = t;
      }
      
      Fitness = OPTIM.ParaEvaluate(false);
      if (IsTracing) PrintParam();
      UpdateOpt(); 
    }
  } else
    CartesianProduct(0,k);
}


/*======================================================================*/
/* CartesianProduct: used in ScanCluster(k) when cluster[k] is linked   */
/*======================================================================*/
void LineSearch :: CartesianProduct(int j, int k) 
{
  int m = (ParaClusters[k])[j];
  
  if ( j == ((int) ParaClusters[k].size()) - 1 ) 
    for (int q=0; q<=DivInter; q++) {
      Para[m] = ParaMinInter[m] + q*ParaStep[m];
      Fitness = OPTIM.ParaEvaluate(false);
      UpdateOpt();
      if (IsTracing) PrintParam();
    }
  else {
    j++;
    for (int q=0; q<=DivInter; q++) {
      Para[m] = ParaMinInter[m] + q*ParaStep[m];
      CartesianProduct(j,k);
    }
  }
}

	  
/*======================================================================*/
/* ChooseOptimal: chooses randomly one optimal point from Optimums vector    */
/*======================================================================*/
void LineSearch :: ChooseOptimal(void) 
{  
  int n = (int) Rand->RandBoundedInt((int) Optimums.size());
  
  for (unsigned int i=0; i<Para.size(); i++) Para[i] = (Optimums[n])[i];
  
  Optimums.clear();
  Optimums.push_back(Para);
  Fitness = FitOpt;
}


/*======================================================================*/
/* ReduceSearch: Reduces interval search for each parameter             */
/*======================================================================*/
void LineSearch :: ReduceSearch(void) {
  
  double  marge, step; 
  
  for (unsigned int i=0; i<Para.size(); i++) {
    // Reduce interval if possible
    ParaLInter[i] *= Alpha;
    marge = ParaLInter[i]/2.;    
    if (ParaLInter[i] > ParaMinStep[i]) {
      ParaMinInter[i] = Max(ParaMin[i],Para[i] - marge);
      ParaMaxInter[i] = Min(ParaMax[i],Para[i] + marge);
    }
    // Update the step to explore the interval

    step = (ParaMaxInter[i] - ParaMinInter[i]) / DivInter;
    ParaStep[i] = Max(step, ParaMinStep[i]);
  }
}



//-------------------------------------------------------
// Print Para values and associated fitness
// BEWARE: just the first para of an IDENTICAL cluster is printed.
//-------------------------------------------------------
void LineSearch::PrintParam(void)
{
  unsigned int i;
  std::vector <bool> FirstEltCluster;

  // parameters value with just the first elt of clusters IDENTICAL
  for (i=0; i<ParaClusters.size(); i++) FirstEltCluster.push_back( true );
  for (i=0; i<ParaName.size(); i++) 
    if (ParaClusterRelations[ParaCluster[i]] == IDENTICAL) {
      if ( FirstEltCluster[ParaCluster[i]] ) {
	std::cout << Para[i] << "\t";
	FirstEltCluster[ParaCluster[i]] = false;
      }
    } else
      std::cout << Para[i] << "\t";

  std::cout << "Fitness = " << Fitness << "\n";
  fflush(stdout); fflush(stderr);
}

