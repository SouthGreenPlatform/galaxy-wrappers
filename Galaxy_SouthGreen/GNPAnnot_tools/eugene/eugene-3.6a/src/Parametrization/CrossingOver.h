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
// $Id: CrossingOver.h,v 1.3 2004-09-16 13:06:35 cros Exp $
// ------------------------------------------------------------------
// File:     CrossingOver.cc
// Contents: Class for cross over of the Genetic algorithm
// ------------------------------------------------------------------

#ifndef CROSSOVER_H_INCLUDED
#define CROSSOVER_H_INCLUDED

#include <vector>

#include "Chromosome.h"


class CrossingOver {
 public:
  double Proba;

  CrossingOver(double proba);
  ~CrossingOver(void);
  void Crosseval(void); 

 private:
  std::vector<Chromosome*> ChildrenPopulation;

  void Crossover(Chromosome* c1, Chromosome* c2);
  void ChooseCouple(int* n1, int* n2);

};

#endif
