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
// $Id: Chromosome.h,v 1.4 2004-09-16 13:06:35 cros Exp $
// ------------------------------------------------------------------
// File:     Chromosome.cc
// Contents: Class Data and Chromosome of the Genetic algorithm
// ------------------------------------------------------------------

#ifndef CHROMOSOME_H_INCLUDED
#define CHROMOSOME_H_INCLUDED


#include <vector>
#include <string>

class Data {
 public:
  std::vector <double> P;

  Data(int n);
  double GetDistance(Data* d);
  void PutBarycenter(int nb1, Data *data2, int nb2);
  void Recenter(void);
};


class Chromosome {
 public:
  double RawFitness;    // Raw fitness
  double ScaledFitness;    
  bool IsEvaluated;     // true if chromosome already evaluated
  bool IsNew;           // true if new chromosome, false if copy of old one
  int NoCluster;        // Number of the cluster this element belongs to
  bool IsBestInCluster; // true if best element in cluster, false otherwise.  
			// Used only with elitist selection for each cluster
  Data *data;

  Chromosome(int n);
  ~Chromosome(void);
  void Evaluate(void);
  void Print(std::string msg);
};



class Cluster {
 public:
    Data* CentralData;  /* Pointer on the class center of this cluster */
    int NbChrom; /* number of elements in this cluster */

    Cluster(int n);
    ~Cluster(void);
};


#endif
