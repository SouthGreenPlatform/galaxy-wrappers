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
// $Id: ParaOptimization.h,v 1.9 2009-04-20 14:11:52 sallet Exp $
// ------------------------------------------------------------------
// File:     ParaOptimization.h
// Contents: The ParaOptimization class optimizes eugene parameters
// ------------------------------------------------------------------

#ifndef PARA_OPTIMIZATION_H_INCLUDED
#define PARA_OPTIMIZATION_H_INCLUDED


#include <vector>
#include <string>
#include "OptiAlgorithm.h"

class MasterSensor;
class DNASeq;
class Prediction;

class ParaOptimization {

 public:
  std::vector <OptiAlgorithm*> Algorithms;
  int                          AlgoIndex;     // index of running algorithm
  std::string  DetailedEvaluation;            // text of detailed evaluation of Para
                                              // updated by ParaEvaluate

  ~ParaOptimization(void);
  void ParaOptimize (int argc, char* argv[]);
  double ParaEvaluate (bool is_detail_required);  

 private:
  std::vector <MasterSensor*>  MSensors;
  std::vector <DNASeq*>        Sequences;
  std::vector <Prediction*>    References;
  std::vector<std::string>     SeqNames;
  std::string                  TrueCoordFile;
  std::string                  ExecutableName; // first part of parameter file
  bool                         IsTest;         // if true then test mode: 
                                      // do not read sequences and related info
                                      // evaluate para with the method NormalLaw
  double	Regularizer;
  void Init(int argc, char* argv[]);
  double NormalLaw(double x);
  void ReadReferences();
};


#endif
