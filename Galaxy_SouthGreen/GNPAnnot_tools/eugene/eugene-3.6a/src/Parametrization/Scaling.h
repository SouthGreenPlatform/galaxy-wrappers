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
// $Id: Scaling.h,v 1.2 2004-09-16 13:06:35 cros Exp $
// ------------------------------------------------------------------
// File:     Scaling.h
// Contents: Class realizing the scaling for the genetic algorithm
// ------------------------------------------------------------------

#ifndef SCALING_H_INCLUDED
#define SCALING_H_INCLUDED

#include "OptiAlgorithm.h"


class Scaling {

 public:
  int Type;

  Scaling(int type);

  void Scalestat(int gen);
};

#endif
