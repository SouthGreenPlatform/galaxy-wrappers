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
// $Id: Sensor.MarkovIMM.h,v 1.11 2006-11-29 14:12:18 tschiex Exp $
// ------------------------------------------------------------------
// File:     Sensor.MarkovIMM.h
// Contents: Sensor MarkovIMM
// ------------------------------------------------------------------

#ifndef  SENSOR_MARKOVIMM_H_INCLUDED
#define  SENSOR_MARKOVIMM_H_INCLUDED

#include "../../Sensor.h"
#include "../0_SensorTk/BStrArray.h"
#include "../0_SensorTk/BStrArray.cc"

/****************************************************************
 **                     SensorMarkovIMM                        **
 ****************************************************************/
class SensorMarkovIMM : public Sensor
{
 private:
  // to avoid to reread the IMM matrices, they are stored in the list IMMatrixList
  // associated to a matrix: 
  //    - the interval of GC% [minGC maxGC] where it is defined
  //    - the request of using model M0 for IG
  // for a given instance, IMMatrix_index allows to know in lists which matrice 
  // and relative information to consider
  static std::vector < std::vector <BString_Array*> > IMMatrixList;  
  static std::vector <std::string> matrixNameList;
  static std::vector <int> refCount;
  const static int  MODEL_LEN = 9;
  const static int  SIMPLE_MODEL_LEN = 6;
  const static int  ALPHABET_SIZE = 4;
  int IMMatrix_index;
  double minGC;
  double maxGC;
  int maxOrder;
  int IntergenicModel;


 public:
  SensorMarkovIMM         (int n, DNASeq *X);
  virtual ~SensorMarkovIMM();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorMarkovIMM * builder0( int n, DNASeq *X) { return new SensorMarkovIMM(n, X);}

// reserve memory for static variables
std::vector < std::vector <BString_Array*> > SensorMarkovIMM::IMMatrixList;
std::vector <std::string>                    SensorMarkovIMM::matrixNameList;
std::vector <int>                            SensorMarkovIMM::refCount;

#endif
