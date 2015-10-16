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
// $Id: OptiAlgorithm.h,v 1.5 2004-09-16 13:06:35 cros Exp $
// ------------------------------------------------------------------
// File:     OptiAlgorithm.h
// Contents: The mother class of optimization algorithms classes
// ------------------------------------------------------------------

#ifndef OPTI_ALGORITHM_H_INCLUDED
#define OPTI_ALGORITHM_H_INCLUDED


#include <vector>
#include <string>

#include "Random.h"


// Definition of cluster relations
enum RELATION {LINKED, IDENTICAL, NO_RELATION};
typedef enum RELATION CLUSTER_RELATION;


class OptiAlgorithm {

 public:
  std::vector <double>      Para;  
  std::vector <std::string> ParaName;
  std::vector <double>      ParaMax;
  std::vector <double>      ParaMin; 
  std::vector <int>         ParaCluster;
  std::vector < std::vector <int> > ParaClusters;
  std::vector <CLUSTER_RELATION>    ParaClusterRelations;
  int     NbParaCluster;
  bool    IsTracing;
  Random* Rand;

  OptiAlgorithm(void);
  virtual ~OptiAlgorithm(void) {};
  virtual void Optimize(bool is_chaining) = 0;
  std::string ReduceName(std::string s);
};




#endif
