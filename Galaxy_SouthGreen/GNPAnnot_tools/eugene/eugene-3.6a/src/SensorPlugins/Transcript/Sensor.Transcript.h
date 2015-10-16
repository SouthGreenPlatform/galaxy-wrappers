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
// $Id: Sensor.Transcript.h,v 1.7 2005-08-30 12:40:07 bardou Exp $
// ------------------------------------------------------------------
// File:     Sensor.Transcript.h
// Contents: Sensor Transcript
// ------------------------------------------------------------------

#ifndef  SENSOR_TRANSCRIPT_H_INCLUDED
#define  SENSOR_TRANSCRIPT_H_INCLUDED

#include "../../Sensor.h"

/*************************************************************
 **                     SensorTranscript                    **
 *************************************************************/
class SensorTranscript : public Sensor
{
 private:
  // proba. of transcription Start/Stop
  double transStart;
  double transStop;

 public:
  SensorTranscript  (int n, DNASeq *X);
  virtual ~SensorTranscript   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *X, int, DATA *);
  virtual void Plot       (DNASeq *X);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorTranscript * builder0(int n, DNASeq *X) {  return new SensorTranscript(n, X); }

#endif
