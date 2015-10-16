/* ------------------------------------------------------------------
   Copyright (C) 2004 INRA <eugene@ossau.toulouse.inra.fr>
  
   This program is open source; you can redistribute it and/or modify
   it under the terms of the Artistic License (see LICENSE file).
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
  
   You should have received a copy of Artistic License along with
   this program; if not, please see http://www.opensource.org
  
   $Id: gdIF.h,v 1.10 2005-08-30 12:31:38 bardou Exp $
   ------------------------------------------------------------------
   File:     gdIF.h
   Contents: Functions for interfacing with libGD
   ------------------------------------------------------------------*/

#ifndef  GDIF_H_INCLUDED
#define  GDIF_H_INCLUDED

void InitPNG   (int x, int y, int offset, int From, int To,
		int olap, int ILen,char *name);
void PlotBarF  (unsigned int nuc, signed char phase, double pos,
		double width, int col);
void PlotBarI  (unsigned int nuc, signed char phase, double pos,
		int width, int col);
void PlotLine  (unsigned int nuc1, unsigned int nuc2, signed char phase1,
		signed char phase2, double pos1, double pos2, int col);
void PlotString(unsigned int nuc, signed char phase, double pos,
		char st[], int col);
void ClosePNG();
void OutputHTMLFileNames(int firstImg, char* html_dir, FILE* OUT=stdout);

#endif
