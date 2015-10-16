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
// $Id: Sensor.h,v 1.21 2010-01-25 17:03:07 sallet Exp $
// ------------------------------------------------------------------
// File:     Sensor.h
// Contents: class Sensor
// ------------------------------------------------------------------

#ifndef  SENSOR_H_INCLUDED
#define  SENSOR_H_INCLUDED

#include <vector>
#include <algorithm>
#include <string.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif

#include "Param.h"
#include "DNASeq.h"
#include "Prediction.h"
#include "Const.h"
#include "SensorIF.h"
#include "GeneFeatureSet.h"
    
extern "C" {
#include "GDIF/gdIF.h"
}

/*************************************************************
 **                        Sensor                           **
 *************************************************************/
class Sensor
{
 private:
  int instanceNumber;

 protected:
  void CheckStart   (DNASeq *,std::vector<int> a, std::vector<int> b);
  void CheckStop    (DNASeq *,std::vector<int> a, std::vector<int> b);
  void CheckSplices (DNASeq *,std::vector<int> a, std::vector<int> b, std::vector<int> c, std::vector<int> d);
  void PlotStop     (int pos,   int phase,  int sure);
  void PlotStart    (int pos,   int phase,  double strength);
  void PlotAcc      (int pos,   int strand, double strength);
  void PlotDon      (int pos,   int strand, double strength);
  void PlotESTHit   (int start, int end,    int strand, int filtered);
  void PlotESTGap   (int start, int end,    int strand, int filtered);
  void PlotBlastHit (int start, int end,    int phase,  int level);
  void PlotBlastGap (int start, int phase1, int end,    int phase2, int level);
  void PlotStartReg (int pos,   int phase,  int color);
  void PlotEndReg   (int pos,   int phase,  int color);
  void PlotRepeat   (int start, int end);
  int  GetNumber    (void) { return instanceNumber; }
 //CN
  std::string inputFormat_;
 public:
  unsigned char type;
  
  Sensor  (int n);
  virtual ~Sensor (void);
  virtual void Init       (DNASeq *X) = 0;
  virtual void GiveInfo   (DNASeq *, int, DATA *) = 0;
  virtual void Plot       (DNASeq *) = 0;
  virtual void PostAnalyse(Prediction *, FILE *) = 0;
};

/*************************************************************
 **                       Signals Object                    **
 *************************************************************/
class Signals
{
 private:
  
 public:
  int   pos;
  int   type;
  int   edge;
  char  *score;
    
  Signals  ();
  Signals  (int, int, int, char*);
  ~Signals ();
  void PrintS ();
  bool operator < (int i) { if (pos < i) return true; else return false; }
};

/*************************************************************
 **                      Contents Object                    **
 *************************************************************/
class Contents
{
 private:
  
 public:
  int    start;
  int    end;
  int    type;
  float *score;
    
  Contents  ();
  Contents  (int, int, int, float*);
  ~Contents ();
  void PrintC ();
};


#endif

