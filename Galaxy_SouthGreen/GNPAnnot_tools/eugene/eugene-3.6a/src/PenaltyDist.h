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
// $Id: PenaltyDist.h,v 1.6 2006-10-06 22:09:25 tschiex Exp $
// ------------------------------------------------------------------
// File:     PenaltyDist.h
// Contents: Definitions for length penalties with bounded arbitrary dist. followed by linear
// cost, includes a minimum and maximum len.
// A class to represent length penalties. The penalty can be arbitrary
// up to a given length and then linearly interpolated for further
// length.  Considering prob. distributions this means we have an
// exponential tail. INFINITE costs are allowed.
// ------------------------------------------------------------------

#ifndef  PENALTYD_H_INCLUDED
#define  PENALTYD_H_INCLUDED

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <math.h>
#include <cstdio>

class PenaltyDist
{
  friend class Track;
  friend class DAG;
 private:

  int MinLen;
  int MaxLen;
  std::vector<double> Distribution;
  std::vector<double> MinDist;
  std::vector<double> Delta;
  double FinalSlope;
  double MaxDelta;

 public:
  PenaltyDist();
  ~PenaltyDist();

  void   UpdateMinDist();
  void   UpdateMaxDelta ();

  inline double GetDelta(int len) {
    if (len < (int)Delta.size()) return Delta[len];
    return MaxDelta;
  }
  double MinPen(int len);
  void   LoadPenaltyDist(char * FileName);
  void   PlotExp(FILE* h,int a,int b);
  void   Dump(FILE* h);
  double Integrate();
  void   Normalize();
  inline double operator [] (int len)
    {
      if (len < MinLen) return -log(0.0);
      if (len - MinLen < (int)Distribution.size())
	return Distribution[len-MinLen];
      return Distribution.back()+(len-MinLen-Distribution.size()+1)*FinalSlope;
    };
};


#endif
