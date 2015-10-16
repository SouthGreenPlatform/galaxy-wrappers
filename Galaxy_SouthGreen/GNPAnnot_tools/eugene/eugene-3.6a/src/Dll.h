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
// $Id: Dll.h,v 1.6 2004-09-16 12:49:31 cros Exp $
// ------------------------------------------------------------------
// File:     Dll.h
// Contents: class SensorLoader
// ------------------------------------------------------------------

#ifndef __DLL_H
#define __DLL_H

#include "Sensor.h"

class SensorLoader
{
 protected:
  void *h;
  const char *err;
  Sensor *(*builder_func)(int n, DNASeq *X);	

 public:
  SensorLoader (const char *fname, const char *func_name=0);
  ~SensorLoader ();
  Sensor *MakeSensor(int n, DNASeq *X);
  const char *LastError () { return err; }
};

#endif
