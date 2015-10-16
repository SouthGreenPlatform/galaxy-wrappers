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
// $Id: Sensor.User.cc,v 1.16 2005-08-30 12:40:07 bardou Exp $
// ------------------------------------------------------------------
// File:     Sensor.User.cc
// Contents: Sensor User
// ------------------------------------------------------------------

#include "Sensor.User.h"
#include "yacc.h"

extern Parameters PAR;

/*************************************************************
 **                        SensorUser                       **
 *************************************************************/

// ----------------------
//  Default constructor.
// ----------------------
SensorUser :: SensorUser (int n, DNASeq *X) : Sensor(n)
{
  // all types are possible.
  type = Type_Any;
}

// ----------------------
//  Default destructor.
// ----------------------
SensorUser :: ~SensorUser ()
{
}

// ----------------------
//  Init user.
// ----------------------
void SensorUser :: Init (DNASeq *X)
{
  char tempname[FILENAME_MAX+1];
  int errflag;

  fprintf(stderr,"Reading user data............");
  
  strcpy(tempname, PAR.getC("fstname"));
  strcat(tempname, ".user");
  errflag = Utilisateur(tempname, &Signals, &Contents);   //prise en compte de donnees utilisateur
    
  if (errflag) {
    fprintf(stderr,"none found\n");
  }
  else fprintf(stderr,"done\n");
}

// -----------------------
//  GiveInfo signal user.
// -----------------------
void SensorUser :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  Util(pos, Signals, d);
  Util(pos, Contents, d);
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorUser :: Plot(DNASeq *X)
{
}

// ------------------
//  Post analyse
// ------------------
void SensorUser :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}
