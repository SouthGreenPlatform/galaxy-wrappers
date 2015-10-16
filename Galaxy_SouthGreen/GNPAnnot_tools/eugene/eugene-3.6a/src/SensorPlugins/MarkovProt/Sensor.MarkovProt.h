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
// $Id: Sensor.MarkovProt.h,v 1.8 2005-08-30 12:39:58 bardou Exp $
// ------------------------------------------------------------------
// File:     Sensor.MarkovProt.h
// Contents: Sensor MarkovProt
// ------------------------------------------------------------------

#ifndef  SENSOR_MARKOVPROT_H_INCLUDED
#define  SENSOR_MARKOVPROT_H_INCLUDED

#include "../../Sensor.h"
#include "../0_SensorTk/BStrArray.h"
#include "../0_SensorTk/markov.h"
#include "../0_SensorTk/markov.cc"

/*************************************************************
 **                     SensorMarkovProt                        **
 *************************************************************/
class SensorMarkovProt : public Sensor
{
 private:
  static int maxorder;
  static int order;
  //  char[4] type;
  //  TabChaine<ChainePROT21,double>* ModeleProt;
  static   TabChaine<ChainePROT21,unsigned short>* ModeleProt; 
  static   TabChaine<ChaineADN,double>* ProbacodonGeneral;
  TabChaine<ChaineADN,double>* Probacodon;
  double GCrate;
  static double minGC;
  static double maxGC;
  static bool IsInitialized;

 public:
  SensorMarkovProt  (int n, DNASeq *X);
  virtual ~SensorMarkovProt   ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorMarkovProt * builder0( int n, DNASeq *X) { return new SensorMarkovProt(n, X);}

int SensorMarkovProt::maxorder;
int SensorMarkovProt::order;
TabChaine<ChainePROT21,unsigned short>* SensorMarkovProt::ModeleProt; 
TabChaine<ChaineADN,double>* SensorMarkovProt::ProbacodonGeneral;
double SensorMarkovProt::minGC;
double SensorMarkovProt::maxGC;
bool SensorMarkovProt::IsInitialized = false;

#endif
