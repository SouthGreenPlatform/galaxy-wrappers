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
// $Id: PenaltyDist.cc,v 1.4 2005-03-17 14:13:18 cros Exp $
// ------------------------------------------------------------------
// File:     PenaltyDist.cc
// Contents: Definitions for length penalties with bounded arbitrary dist. followed by linear
// cost, includes a minimum and maximum len.
// ------------------------------------------------------------------


#include <assert.h>
#include <stdio.h>

#include "PenaltyDist.h"

#include "Param.h"
#include "System.h"

extern Parameters   PAR;

// ----------------------------------------------------------------
// Default constructor.
// Allocate a flat zero penalty
// ----------------------------------------------------------------
PenaltyDist :: PenaltyDist()
{
  MinLen = 0;
  MaxLen = 0;
  FinalSlope = 0.0;
  Distribution.clear();
  Distribution.push_back(0.0);
  UpdateMinDist();
  UpdateMaxDelta();
}
// ----------------------------------------------------------------
// Computes the minimum penalty over all sequences LONGER than a given
// length
// ----------------------------------------------------------------
void PenaltyDist :: UpdateMinDist()
{
  double MinVal = -NINFINITY;

  MinDist.resize(Distribution.size());

  for (int i = Distribution.size()-1; i >= 0; i--) {
    MinVal = Min(MinVal,Distribution[i]);
    MinDist[i] = MinVal;
  }
}
// ----------------------------------------------------------------
// returns the maximum penalty decrease one can get between one point
// and a farther point.
// ----------------------------------------------------------------
void PenaltyDist :: UpdateMaxDelta ()
{
  // The computation is carried on the corrected penalty (correcting
  // for FinalSlope). Therefore we only need to look to the explicit
  // Distribution (not the tail neither the forbidden region).

  MaxDelta = 0.0;
  Delta.resize(Distribution.size());

  for (int i = 0; i < (int)Distribution.size()-1; i++)
    for (int j = i+1; j<(int)Distribution.size(); j++) {
      MaxDelta = Max(MaxDelta,
		       (Distribution[i]-i*FinalSlope)-
		       (Distribution[j]-j*FinalSlope));

      Delta[j-i] = Max(Delta[j-i],
		       (Distribution[i]-i*FinalSlope)-
		       (Distribution[j]-j*FinalSlope));
    }
}
// ----------------------------------------------------------------
// Returns the minimum penalty over all sequences LONGER than a given
// length
// ----------------------------------------------------------------
double PenaltyDist :: MinPen(int len)
{
  if (len < MinLen) return MinDist[0];
  
  if (len - MinLen < (int)Distribution.size())
    return MinDist[len-MinLen];

  return Distribution.back()+(len-MinLen-Distribution.size()+1)*FinalSlope;
}
// ----------------------------------------------------------------
// Reads a penalty from file Current format: Max len cost pairs with
// increasing len and which are (linearly interpolated)
// ----------------------------------------------------------------
void PenaltyDist :: LoadPenaltyDist(char *FileName)
{
  FILE* Handle;
  int FirstLen,Len,LastLen;
  double FirstPen,Pen,LastPen;
  char *DirName = new char[FILENAME_MAX+1];

  strcpy(DirName, PAR.getC("eugene_dir"));
  strcat(DirName, MODELS_DIR);
  Handle = FileOpen(DirName, FileName,"r");

  Distribution.clear();
  
  //read the 1st pair which indicates MinLen
  fscanf(Handle,"%d %lf\n",&FirstLen,&FirstPen);
  MinLen = FirstLen;
  Distribution.push_back(FirstPen);

  //read the next pair
  if (fscanf(Handle,"%d %lf\n",&Len,&Pen) != 2)
    {
      fprintf(stderr, "Insufficient distribution data. Cannot compute final slope\n");
      exit(2);
    }

  while (fscanf(Handle,"%d %lf\n",&LastLen,&LastPen) == 2) {

    if (LastLen <= Len)  {
      fprintf(stderr, "Wrong distribution. Non increasing length: %d %lf\n",LastLen,LastPen);
      exit(2);
    }

    FinalSlope = (Pen-FirstPen)/(Len-FirstLen);
    for (int i = 1; i <= Len-FirstLen; i++) {
      Distribution.push_back(FirstPen+i*FinalSlope);
    }
    FirstLen = Len; FirstPen = Pen;
    Pen = LastPen; Len = LastLen;
  }
  
  FinalSlope = (Pen-FirstPen)/(Len-FirstLen);
  MaxLen = FirstLen;
  UpdateMinDist();
  UpdateMaxDelta();

  Normalize();

  fclose(Handle);

  delete [] DirName;
}
// ----------------------------------------------------------------
//  Default destructor.
// ----------------------------------------------------------------
PenaltyDist :: ~PenaltyDist  ()
{
  Distribution.clear();
}
// ----------------------------------------------------------------
// Plot the exp of a len dist 
// ----------------------------------------------------------------
void PenaltyDist :: PlotExp(FILE* h,int a,int b)
{
  for (int i = a; i<= b; i++)
    fprintf(h,"%d %lf\n",i,exp(-(*this)[i]));
}
// ----------------------------------------------------------------
// Dumps a len dist 
// ----------------------------------------------------------------
void PenaltyDist :: Dump(FILE* h)
{
  for (unsigned int i = 0; i < Distribution.size(); i++)
    fprintf(h,"%d %lf\n",i+MinLen,Distribution[i]);
}
// ----------------------------------------------------------------
// Computes the integral of the distribution of the exponential of the
// penalty. On integre betement sur Distribution puis
// analytiquement sur la queue
// ----------------------------------------------------------------
double PenaltyDist :: Integrate()
{
  double Norm = 0.0;
  double a,b;
  unsigned int i;
  
  a = exp(-Distribution[0]);
  for (i = 1; i < Distribution.size(); i++) {
    b = exp(-Distribution[i]);
    Norm += (a+b)/2;
    a = b;
  }

  // Now handle the tail. We have a = exp(-lastpenalty) and i =
  // DistSize
  b = (MaxLen <0 ? -NINFINITY : MaxLen);
  Norm += a*(1-exp(-FinalSlope*(b-i+1)))/FinalSlope;
  
  return Norm;
}
// ----------------------------------------------------------------
// Normalizes a PenaltyDistribution considered as a log probability 
// ----------------------------------------------------------------
void PenaltyDist :: Normalize()
{
  if ((Distribution.size() == 1) && (FinalSlope == 0.0)){
    Distribution[0] = 0.0;
    return;
  }

  double Norm = this->Integrate();
  for (unsigned int i = 0; i < Distribution.size(); i++)
    Distribution[i] += log(Norm);
}
