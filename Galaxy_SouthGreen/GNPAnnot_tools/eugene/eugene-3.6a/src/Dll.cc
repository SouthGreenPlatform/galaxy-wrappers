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
// $Id: Dll.cc,v 1.8 2004-09-16 12:49:31 cros Exp $
// ------------------------------------------------------------------
// File:     Dll.cc
// Contents: class SensorLoader
// ------------------------------------------------------------------

#include <dlfcn.h>

#include "Dll.h"

SensorLoader :: SensorLoader(const char *fname, const char *builder)
{
  // Try to open the library now and get any error message.
  h = dlopen( fname, RTLD_NOW );
  err = dlerror();
  if (err) printf("dlerror: %s\n",err);
  
  // Try get the creation function if there is no error yet
  builder_func=0;

  if( LastError()==0 && h) {
    builder_func = (Sensor *(*)(int, DNASeq *))dlsym( h, builder ? builder : "builder0" );
    err=dlerror();
  }
}

SensorLoader :: ~SensorLoader()
{
  // Close the library if it isn't null
  if( h != 0 )
    dlclose(h);
}

Sensor* SensorLoader :: MakeSensor(int n, DNASeq *X)
{
  return (* builder_func)(n, X);
}
