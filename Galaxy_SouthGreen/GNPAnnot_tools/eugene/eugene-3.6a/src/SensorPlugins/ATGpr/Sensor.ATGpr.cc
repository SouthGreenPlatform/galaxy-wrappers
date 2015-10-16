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
// $Id: Sensor.ATGpr.cc,v 1.8 2009-01-12 14:14:21 sallet Exp $
// ------------------------------------------------------------------
// File:     Sensor.ATGpr.cc
// Contents: Sensor ATGpr 
// ------------------------------------------------------------------

#include "Sensor.ATGpr.h"

extern Parameters PAR;

#define NORM(x,n) (((n)+(Max(-(n),x)))/(n))

/*************************************************************
 **                        SensorATGpr                      **
 *************************************************************/

// ----------------------
// Default constructor.
// ----------------------
SensorATGpr :: SensorATGpr (int n, DNASeq *X) : Sensor(n)
{
  char tempname[FILENAME_MAX+1];

  type = Type_Start;
  
  vPosF.clear();
  vValF.clear();
  vPosR.clear();
  vValR.clear();
  
  fprintf(stderr, "Reading start file (ATGpr)....................");
  fflush(stderr);
  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".atgpr");
  ReadATGprF(tempname, X->SeqLen);
  fprintf(stderr,"forward,");
  fflush(stderr);
  
  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".atgprR");
  ReadATGprR(tempname, X->SeqLen);
  fprintf(stderr," reverse done\n");
  
  CheckStart(X,vPosF, vPosR);

  // vectors for reverse are put in the increasing order
  reverse(vPosR.begin(), vPosR.end()); 
  reverse(vValR.begin(), vValR.end()); 
}

// ----------------------
//  Default destructor.
// ----------------------
SensorATGpr :: ~SensorATGpr ()
{
  vPosF.clear();
  vValF.clear();
  vPosR.clear();
  vValR.clear();
}

// ----------------------
//  Init start.
// ----------------------
void SensorATGpr :: Init (DNASeq *X)
{
  startP = PAR.getD("ATGpr.startP*",GetNumber());
  startB = PAR.getD("ATGpr.startB*",GetNumber());
  
  indexF = indexR = 0;
  PositionGiveInfo = -1;

  if (PAR.getI("Output.graph")) Plot(X);
}

// --------------------------
//  Read start forward file.
// --------------------------
void SensorATGpr :: ReadATGprF (char name[FILENAME_MAX+1], int Len)
{
  FILE *fp;
  int i,j = -1;
  double force;
  
  fp = FileOpen(NULL, tempname, "r", PAR.getI("EuGene.sloppy"));
  if (!fp) {
    exit(2);
  }
  
  while (1) {
    i = fscanf(fp,"%*s %lf %*s %*s %d %*s %*s %*s\n", &force, &j); 
    if (i == EOF) break;
    if (i < 2) {
      fprintf(stderr, "Error in ATGpr file %s, position %d\n", name, j);
      exit(2);
    }
    vPosF.push_back( j-1 );
    vValF.push_back( force );
  }
  if (j == -1) fprintf(stderr,"WARNING: empty ATGpr file !\n");
  fclose(fp);
}

// --------------------------
//  Read start reverse file.
// --------------------------
void SensorATGpr :: ReadATGprR (char name[FILENAME_MAX+1], int Len)
{
  FILE *fp;
  int i,j = -1;
  double force;
  
  if (!(fp = fopen(name, "r"))) {
    fprintf(stderr, "cannot open start file %s\n", name);
    exit(2);
  }

  while (1) {
    i = fscanf(fp,"%*s %lf %*s %*s %d %*s %*s %*s\n", &force, &j);
    if (i == EOF) break;
    if (i < 2) {
      fprintf(stderr, "Error in ATGpr file %s, position %d\n", name, j);
      exit(2);
    }

    //if (force > 0.1) {
    j = Len-j+2;
    vPosR.push_back( j-1 );
    vValR.push_back( force );
    //}
  }
  if (j == -1) fprintf(stderr,"WARNING: empty ATGpr file !\n");
  fclose(fp);
}

// ------------------------
//  GiveInfo signal start.
// ------------------------
void SensorATGpr :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  bool update = false;
  double f;
  
  // update indexes on vectors
  if ( (PositionGiveInfo == -1) || (pos != PositionGiveInfo+1) ) update = true;
  PositionGiveInfo = pos;
  
  // Start Forward
  if(!vPosF.empty()) {
    if (update) 
      indexF = lower_bound(vPosF.begin(), vPosF.end(), pos)-vPosF.begin();
    
    if((indexF<(int)vPosF.size()) && (vPosF[indexF] == pos)) {
      // correct read score  to be in ]0 inf[
      vValF[indexF] = fabs(vValF[indexF]); 
      if (vValF[indexF]==0) vValF[indexF]=0.001;

      f = pow(vValF[indexF], startB) * exp(-startP);
      d->sig[DATA::Start].weight[Signal::Forward] += log(f);
      d->sig[DATA::Start].weight[Signal::ForwardNo] += log(1.0-f);
      indexF++;
    }
  }
  
  // Start Reverse
  if (!vPosR.empty()) {
    if (update) 
      indexR = lower_bound(vPosR.begin(),vPosR.end(),pos)-vPosR.begin();
    
    if((indexR<(int)vPosR.size()) && (vPosR[indexR] == pos)) {
      // correct read score  to be in ]0 inf[
      if (vValR[indexR]<0) vValR[indexR] = -vValR[indexR];
      if (vValR[indexR]==0) vValR[indexR]=0.001;

      f = pow(vValR[indexR], startB) * exp(-startP);
      d->sig[DATA::Start].weight[Signal::Reverse]   += log(f);
      d->sig[DATA::Start].weight[Signal::ReverseNo] += log(1.0-f);
      indexR++;
    }
  }
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorATGpr :: Plot(DNASeq *X)
{
  for (int i =0; i < (int)vPosF.size(); i++)
    PlotStart(vPosF[i],(vPosF[i]%3)+1,NORM(log(vValF[i]),4.0));

  for (int i =0; i < (int)vPosR.size(); i++)
    PlotStart(vPosR[i],-((X->SeqLen-vPosR[i])%3)-1,NORM(log(vValR[i]),4.0));
}

// ------------------
//  Post analyse
// ------------------
void SensorATGpr :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}
