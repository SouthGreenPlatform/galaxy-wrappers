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
// $Id: MSensor.h,v 1.17 2005-08-30 08:56:13 bardou Exp $
// ------------------------------------------------------------------
// File:     MSensor.h
// Contents: classes UseSensor, MasterSensor
// ------------------------------------------------------------------

#ifndef  MSENSOR_H_INCLUDED
#define  MSENSOR_H_INCLUDED

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstdio>
#include <vector>
#include <algorithm>
#include <string>
#include <string.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif

#include "Sensor.h"
#include "DNASeq.h"
#include "Param.h"
#include "Const.h"
#include "Dll.h"

/*************************************************************
 **                        UseSensor                        **
 *************************************************************/
class UseSensor
{
 private:
  
 public:
  int  Priority;
  char Name[FILENAME_MAX+1];

  UseSensor  ();
  UseSensor  (int, char[FILENAME_MAX+1]);
};

/*************************************************************
 **                      MasterSensor                       **
 *************************************************************/
class MasterSensor
{
 private:
  std::vector <Sensor*>    theSensors;          // list of created instances of sensors
  static bool              IsInitialized;       // when true MSSensorsList is initialized
  static std::string       PluginsDir;          // path of the plugins dir
  static std::vector <std::string>   MSSensorsList;    // list of sensors defined at top level
  static std::vector <std::string>   LoadedSensorsList;// list of all loaded sensors (included 
                                                       // sensors loaded for IfElse, Tester)
  static std::vector <SensorLoader*> dllList;   // list of dll (loaded .so)
  int  LoadSensor   (std::string name);         // load the .so if it is not yet loaded
  
 public:
  ~MasterSensor (void);
  void InitMaster   (DNASeq *X);
  void InitSensors  (DNASeq *X);
  void GetInfoAt    (DNASeq *X, int pos, DATA *d);
  void PrintDataAt  (DNASeq *X, int pos, DATA *d,   FILE *OUT=stdout);
  int  GetInfoSpAt  (unsigned char type, DNASeq *X, int pos, DATA *d);
  void PostAnalyse  (Prediction *pred,   FILE *MISC_INFO);
  Sensor* MakeSensor (std::string name, int n, DNASeq *X);  // create an instance of a sensor and 
                                      // load before the corresponding .so if it is not yet loaded
};



#endif
