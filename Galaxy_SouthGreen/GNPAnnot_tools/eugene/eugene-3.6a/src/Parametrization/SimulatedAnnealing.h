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
// $Id: SimulatedAnnealing.h,v 1.2 2004-09-16 13:06:35 cros Exp $
// ------------------------------------------------------------------
// File:     SimulatedAnnealing.h
// Contents: Class realizing the simulated annealing for the genetic algorithm
// ------------------------------------------------------------------

#ifndef SIMULTED_ANNEALING_H_INCLUDED
#define SIMULTED_ANNEALING_H_INCLUDED


#include "Random.h"
#include "Chromosome.h"


class SimulatedAnnealing {

 public:
  SimulatedAnnealing(void);
  void SAInit(void);
  void SAUpdate(int gen);
  void SATournament4(Chromosome *C1,Chromosome *C2, 
		     Chromosome *P1, Chromosome *P2);
  void SATournament2(Chromosome* c1, Chromosome* c2);

 private:
  int CP;
  double CC1;
  double CC2;
  double CC;  /* CC : weight value for temperature modification during */
              /*      the annealing scheme : T(k+1)=CC*T(k) where k    */
              /*      is the current generation number. Two annealing  */
              /*      schemes are used during the temperature decrease:*/
              /*      one with CC1 and the other with CC2.             */
  double Tx;  /* we change the annealing scheme when we reach */
              /* the temperature Tx                           */
  double T; /* re-annealing temperature */
};

#endif
