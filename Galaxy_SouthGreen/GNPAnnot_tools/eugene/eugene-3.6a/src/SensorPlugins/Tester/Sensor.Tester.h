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
// $Id: Sensor.Tester.h,v 1.14 2008-07-02 11:20:31 tschiex Exp $
// ------------------------------------------------------------------
// File:     Sensor.Tester.h
// Contents: Sensor Tester
// ------------------------------------------------------------------

#ifndef  SENSOR_TESTER_H_INCLUDED
#define  SENSOR_TESTER_H_INCLUDED


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif

#include "../../Sensor.h"
#include "../../Dll.h"
#include "../../Prediction.h"


/*************************************************************
 **                     SensorTester                        **
 *************************************************************/
enum EVALUATED_TYPE {EVALUATED_START, EVALUATED_STOP, EVALUATED_SPLICE};

class SensorTester : public Sensor
{
 private:

  static EVALUATED_TYPE EvaluatedType;
  static bool IsSPSN;
  static int MinNumbers;
  static std::string SensorName;
  static std::string Todo;
  static int SensorInstance;

  static bool  IsInitialized;

  static int NbToDoSequence;
  static int NbDoneSequence;
  static std::vector<double> Thresholds; // list of found scores 
  static std::vector<std::vector <std::vector <double> > > Scores; // Score[k][i][j] is 
  // the score of sequence k, at position i, for signal-strand j (SIGNAL_STRAND_NO)
  static std::vector<std::vector <std::vector <bool> > > IsAnnotateds;
  // IsAnnotateds[k][i][j] is true if for sequence k, at position i,
  // there is an annotation of type j (SIGNAL_STRAND_NO)
  static std::vector<std::vector <std::vector <bool> > > IsAPositions;
  // IsAPositions[k][i][j] is true if for sequence k, at position i
  // there is a conventional coding of a site j (SIGNAL_STRAND_NO)
  static std::vector <std::vector <int> > TP, FP, TN, FN;  // True/False Positive/Negative
  // TP[n][l] is the number of TP considering Thresholds[n] as threshold for positive
  // all for l=0 or acc only for l=0, don only for l=1
  static std::vector <std::vector <int> > Nb; // number of detections
  // Nb[n][l] is the number of score egal Thresholds[n] 
  // all for l=0 or acc only for l=0, don only for l=1
  Sensor       *sensor;
  Prediction   *gene;
  char         seqName[FILENAME_MAX+1];
  static FILE  *fp;
  enum SIGNAL_STRAND_NO {START_F=0, START_R=1, STOP_F=0, STOP_R=1, ACC_F=0, ACC_R=1, DON_F=2, DON_R=3};
  unsigned char SensorType;

  void ReadGFFAnnotation(char name[FILENAME_MAX+1], DNASeq *x);
  void ReadGFF3Annotation(GeneFeatureSet & geneFeatureSet , DNASeq *X);
  void ReadAnnotation(DNASeq *);
  const char* SigType_TF(int, int, const char **);
  const char* State(int);

  void AnalyzeSPSN(void);
  void UpdateTP_FP_TN_FN(int no_threshold, bool is_don, bool is_annotated, bool is_detected);
  void UpdateThreshold(double t);
  void print(void);
 public:
  SensorTester (int n, DNASeq *X);
  virtual ~SensorTester   (void);
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorTester * builder0(int n, DNASeq *X) { return new SensorTester(n, X); }


// declare static variables
EVALUATED_TYPE SensorTester::EvaluatedType;
bool SensorTester::IsSPSN;
int  SensorTester::MinNumbers;
std::string SensorTester::SensorName;
std::string SensorTester::Todo;
int  SensorTester::SensorInstance;
bool SensorTester::IsInitialized;
int  SensorTester::NbToDoSequence;
int  SensorTester::NbDoneSequence;
std::vector<double> SensorTester::Thresholds;
std::vector<std::vector <std::vector <double> > > SensorTester::Scores; 
std::vector<std::vector <std::vector <bool> > > SensorTester::IsAnnotateds;
std::vector<std::vector <std::vector <bool> > > SensorTester::IsAPositions;
std::vector<std::vector <int> > SensorTester::TP, SensorTester::FP, SensorTester::TN, SensorTester::FN;  
std::vector<std::vector <int> > SensorTester::Nb;
FILE *SensorTester::fp;



#endif
