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
// $Id: Sharing.h,v 1.2 2004-09-16 13:06:35 cros Exp $
// ------------------------------------------------------------------
// File:     Sharing.h
// Contents: Class realizing the sharing for the genetic algorithm
// ------------------------------------------------------------------

#ifndef SHARING_H_INCLUDED
#define SHARING_H_INCLUDED


#include "OptiAlgorithm.h"


class Sharing {

 public:
  double Value;

  Sharing(double value);
  void Share (void);

 private:
  double EvalShare(double distance);

};

#endif
