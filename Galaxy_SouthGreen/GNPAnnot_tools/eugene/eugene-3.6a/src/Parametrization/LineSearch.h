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
// $Id: LineSearch.h,v 1.4 2004-09-16 13:06:35 cros Exp $
// ------------------------------------------------------------------
// File:     LineSearch.h
// Contents: Class of the Line Search algorithm
// ------------------------------------------------------------------

#ifndef LINESEARCH_H_INCLUDED
#define LINESEARCH_H_INCLUDED


#include "OptiAlgorithm.h"


class LineSearch : public OptiAlgorithm {

 private:
  std::vector <double> ParaMinInter;
  std::vector <double> ParaMaxInter;
  std::vector <double> ParaLInter;
  std::vector <double> ParaStep;
  std::vector <double> ParaMinStep;
  int NbMaxCycle;
  int NbMinCycle;
  int NbMaxStab;
  int DivInter;
  double Alpha;
  double EvolutionMini;
  double FitOpt;           // fitness of the current optimal point
  double Fitness;          // fitness of the point described in Para
  std::vector < std::vector <double> > Optimums;  // list of optimal points after
                                                  // the scan of a cluster
  std::string MsgParaNames;

  void ReduceSearch(void);
  void ScanCluster(int k);
  void CartesianProduct(int j, int k);
  void ChooseOptimal(void);
  void UpdateOpt(void);
  void PrintParam(void);
  
 public:

  LineSearch(void);

  void Optimize(bool is_chaining);
};


#endif
