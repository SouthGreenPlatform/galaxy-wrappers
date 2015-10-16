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
// $Id: Selecting.h,v 1.2 2004-09-16 13:06:35 cros Exp $
// ------------------------------------------------------------------
// File:     Selecting.h
// Contents: Class realizing the selection for the genetic algorithm
// ------------------------------------------------------------------

#ifndef SELECTING_H_INCLUDED
#define SELECTING_H_INCLUDED


class Selecting {

 public:
  int Type;

  Selecting(int type);
  ~Selecting(void);
  void Selection(void);

 private:
  double* tab;
  double* fraction;
  int* choices;
  int nremain;

};

#endif
