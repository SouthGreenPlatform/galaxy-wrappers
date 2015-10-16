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
// $Id: Sensor.User.h,v 1.14 2005-08-30 12:40:07 bardou Exp $
// ------------------------------------------------------------------
// File:     Sensor.User.h
// Contents: Sensor User
// ------------------------------------------------------------------

#ifndef  SENSOR_USER_H_INCLUDED
#define  SENSOR_USER_H_INCLUDED

#include "../../Sensor.h"
#include "structure.h"

/*************************************************************
 **                       SensorUser                        **
 *************************************************************/
class SensorUser : public Sensor
{
 private:
  
  ptUTIL Signals;
  ptUTIL Contents;
  
 public:
  SensorUser  (int n, DNASeq *X);
  virtual ~SensorUser     ();
  virtual void Init       (DNASeq *);
  virtual void GiveInfo   (DNASeq *, int, DATA *);
  virtual void Plot       (DNASeq *);
  virtual void PostAnalyse(Prediction *, FILE *);
};

extern "C" SensorUser * builder0( int n, DNASeq *X) {  return new SensorUser(n, X);}

#endif
