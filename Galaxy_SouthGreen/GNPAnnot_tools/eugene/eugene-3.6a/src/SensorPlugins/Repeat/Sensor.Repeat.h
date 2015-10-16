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
// $Id: Sensor.Repeat.h,v 1.17 2007-04-19 09:35:32 cnoirot Exp $
// ------------------------------------------------------------------
// File:     Sensor.Repeat.h
// Contents: Sensor Repeat
// ------------------------------------------------------------------

#ifndef  SENSOR_REPEAT_H_INCLUDED
#define  SENSOR_REPEAT_H_INCLUDED

#include "../../Sensor.h"


/*************************************************************
 **                     SensorRepeat                        **
 *************************************************************/
class SensorRepeat : public Sensor
{
 private:
  int PositionGiveInfo;

  std::vector<int>           vDeb;
  std::vector<int>           vFin;
  int index;
  double intronPenalty;
  double UTRPenalty;
  double exonPenalty;
  void ReadRepeat (char *name, int SeqLen);
  void ReadRepeatGff3 (char *name, int SeqLen);
 public:
  SensorRepeat  (int n, DNASeq *X);
  virtual ~SensorRepeat   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorRepeat * builder0( int n, DNASeq *X) { return new SensorRepeat(n, X);}

#endif
