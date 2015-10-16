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
// $Id: Sensor.SpliceWAM.h,v 1.5 2005-08-30 12:40:04 bardou Exp $
// ------------------------------------------------------------------
// File:     Sensor.SpliceWAM.h
// Contents: Sensor SpliceWAM
// A simple splice site detection sensor based on a Weight Array Model
// ------------------------------------------------------------------


#ifndef  SENSOR_SPLICEWAM_H_INCLUDED
#define  SENSOR_SPLICEWAM_H_INCLUDED

#include "../../Sensor.h"
#include "../0_SensorTk/WAM.cc"

/*************************************************************
 **                 SensorSpliceWAM                    **
 *************************************************************/
class SensorSpliceWAM : public Sensor
{
 private:
  double AccScaleCoef;
  double AccScalePenalty;
  double DonScaleCoef; // coef for the WAM score
  double DonScalePenalty; //  penalty for the WAM score (score= ScaleCoef * WAMscore - ScalePenalt
  static int NbNtBeforeGT, NbNtAfterGT, DonorSiteLength;
  static int NbNtBeforeAG, NbNtAfterAG, AcceptorSiteLength;
  static int MarkovianOrder;    // MarkovianOrder of the markov models in the Weight Array Model
  static WAM* DonWAModel;
  static WAM* AccWAModel;
  static bool IsInitialized;

 public:
  SensorSpliceWAM  (int n, DNASeq *X);
  virtual ~SensorSpliceWAM   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorSpliceWAM * builder0(int n, DNASeq *X) {  return new SensorSpliceWAM(n, X); }

int SensorSpliceWAM::NbNtBeforeGT, SensorSpliceWAM::NbNtAfterGT, SensorSpliceWAM::DonorSiteLength;
int SensorSpliceWAM::NbNtBeforeAG, SensorSpliceWAM::NbNtAfterAG, SensorSpliceWAM::AcceptorSiteLength;
int SensorSpliceWAM::MarkovianOrder;   
WAM* SensorSpliceWAM::DonWAModel;
WAM* SensorSpliceWAM::AccWAModel;
bool SensorSpliceWAM::IsInitialized = false;


#endif
