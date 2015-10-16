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
// $Id: Sensor.StartWAM.h,v 1.6 2005-08-30 12:40:06 bardou Exp $
// ------------------------------------------------------------------
// File:     Sensor.StartWAM.h
// Contents: Sensor StartWAM
// Definition of a start codon detection sensor based on a Weight Array Model
// ------------------------------------------------------------------

#ifndef  SENSOR_STARTWAM_H_INCLUDED
#define  SENSOR_STARTWAM_H_INCLUDED

#include "../../Sensor.h"
#include "../0_SensorTk/WAM.cc"

/*************************************************************
 **                 SensorStartWAM                    **
 *************************************************************/
class SensorStartWAM : public Sensor
{
 private:
  double ScaleCoef; // coefficient for the WAM score scaling
  double ScalePenalty; //  penality for the WAM score scaling
  static int NbNtBeforeATG;
  static int NbNtAfterATG;
  static int MotifLength;
  static int MarkovianOrder; // order of the markov models in the WAM
  static double PlotScoreIncrease;
  static WAM* WAModel;
  static bool IsInitialized;

  double ScaleWAMScore (double WAMScore);  //  (score= coef * WAMscore - pen)
  inline double NormalizePlot (double x, double n);  // normalization for the plot

 public:
  SensorStartWAM  (int n, DNASeq *X);
  virtual ~SensorStartWAM   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorStartWAM * builder0(int n, DNASeq *X) {  return new SensorStartWAM(n, X); }

int SensorStartWAM::NbNtBeforeATG;
int SensorStartWAM::NbNtAfterATG;
int SensorStartWAM::MotifLength;
int SensorStartWAM::MarkovianOrder; 
double SensorStartWAM::PlotScoreIncrease;
WAM* SensorStartWAM::WAModel;
bool SensorStartWAM::IsInitialized = false;


#endif
