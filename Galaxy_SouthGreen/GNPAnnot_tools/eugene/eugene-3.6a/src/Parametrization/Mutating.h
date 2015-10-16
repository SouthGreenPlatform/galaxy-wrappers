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
// $Id: Mutating.h,v 1.4 2004-12-03 15:41:44 tschiex Exp $
// ------------------------------------------------------------------
// File:     Mutating.h
// Contents: Class realizing mutation for the genetic algorithm
// ------------------------------------------------------------------

#ifndef MUTATING_H_INCLUDED
#define MUTATING_H_INCLUDED


#include <vector>
#include <float.h>

#include "Chromosome.h"


class Mutating {

 public:
  double Proba;

  Mutating(double proba);
  ~Mutating(void);
  void Muteval (void); 

 private:
  std::vector<Chromosome*> ChildrenPopulation;
  void ChooseChromosome(int* n);
  void Mutation(Chromosome* c);
};

#endif
