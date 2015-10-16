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
// $Id: Sensor.Repeat.cc,v 1.33 2009-01-12 14:38:55 sallet Exp $
// ------------------------------------------------------------------
// File:     Sensor.Repeat.cc
// Contents: Sensor Repeat
// ------------------------------------------------------------------

#include "Sensor.Repeat.h"

extern Parameters PAR;

/*************************************************************
 **                        SensorRepeat                     **
 *************************************************************/

// ----------------------
// Default constructor.
// ----------------------
SensorRepeat :: SensorRepeat (int n, DNASeq *X) : Sensor(n)
{
  char tempname[FILENAME_MAX+1];
  FILE* ncfile;
  int deb, fin;
  char line[MAX_LINE];
  int read;

  type = Type_Content;

  fprintf(stderr,"Reading Intergenic regions....");
  fflush(stderr);
  
  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".ig");
  
  inputFormat_ = to_string(PAR.getC("Repeat.format", GetNumber(),1));

  if ( inputFormat_ == "GFF3" )
  {
    strcat(tempname,".gff3");
    ReadRepeatGff3(tempname, X->SeqLen);
  }
  else
  {
    ReadRepeat(tempname,X->SeqLen);
  }
}


// -----------------------------
//  Read Repeat region file
// -----------------------------
void SensorRepeat :: ReadRepeat (char *name, int SeqLen) 
{
  FILE* ncfile;
  int deb, fin;
  char line[MAX_LINE];
  int read;

  ncfile = FileOpen(NULL,name, "r",PAR.getI("EuGene.sloppy"));

  if (ncfile) { 

    while (fgets(line, MAX_LINE, ncfile)!= NULL) {

      read = sscanf(line, "%d %d %*s\n", &deb, &fin);

      int s = (int)vDeb.size();
      deb   = Max(1,deb)-1;
      fin   = Min(SeqLen,fin)-1;
      if(deb > fin  || (read <2) ||
	 (s != 0  &&  vFin[s-1] >= deb)) {
	fprintf(stderr,"\nError in ig file %s, line %d\n", name, s+1);
	exit(2);
	 }
	 vDeb.push_back( deb );
	 vFin.push_back( fin );
    }
    fprintf(stderr,"done\n");
  } else
  {
    fprintf(stderr,"no file\n");
  }
}

// -----------------------------
//  Read Repeat region gff3 file
// -----------------------------
void SensorRepeat :: ReadRepeatGff3 (char *name, int SeqLen) 
{
  int deb, fin;
  GeneFeatureSet * geneFeatureSet = new GeneFeatureSet (name);
  vector< GeneFeature *>::iterator it = geneFeatureSet->getIterator();
  int nbElement=geneFeatureSet->getNbFeature();
  //geneFeatureSet->printFeature();
  int i=0;
  while ( i<nbElement )
  {
    //(*it)->second();
    GeneFeature * tmpFeature = *it;
    string idSo=tmpFeature->getType();
    if ( idSo.find("SO:") == string::npos )
    {
      string tmp=GeneFeatureSet::soTerms_->getIdFromName(idSo);
      idSo=tmp;
    }
    int s = (int)vDeb.size();
    deb   = Max(1,tmpFeature->getLocus()->getStart())-1;
    fin   = Min(SeqLen,tmpFeature->getLocus()->getEnd())-1;
    if(deb > fin || (s != 0  &&  vFin[s-1] >= deb) || strcmp (idSo.c_str(),"SO:0000657")!=0) 
    {
      fprintf(stderr,"\nError in ig file %s, line %d\n", name, s+1);
      exit(2);
    }
    vDeb.push_back( deb );
    vFin.push_back( fin );
    it++;
    i++;
  }
  fprintf(stderr,"done\n");
  delete geneFeatureSet;
}

// ----------------------
//  Default destructor.
// ----------------------
SensorRepeat :: ~SensorRepeat ()
{
  vDeb.clear();
  vFin.clear();
}

// ----------------------
//  Init Repeat.
// ----------------------
void SensorRepeat :: Init (DNASeq *X)
{
  UTRPenalty = PAR.getD("Repeat.UTRPenalty*");
  intronPenalty = PAR.getD("Repeat.IntronPenalty*");
  exonPenalty = PAR.getD("Repeat.ExonPenalty*");
  
  index = 0;
  PositionGiveInfo = -1;

  if (PAR.getI("Output.graph")) Plot(X);
}

// --------------------------
//  GiveInfo Content Repeat.
// --------------------------
void SensorRepeat :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  bool update = false;
  if ( (PositionGiveInfo == -1) || (pos != PositionGiveInfo+1) ) update = true; // update indexes on vectors
  PositionGiveInfo = pos;

  if(!vFin.empty()) {
    if (update) 
      index = lower_bound(vFin.begin(), vFin.end(), pos) - vFin.begin();

    if (index < (int) vFin.size())
      if (pos >= vDeb[index]) { 
	for(int i=DATA::ExonF1; i<=DATA::ExonR3; i++)	 // Exon(6)
	  d->contents[i] -= exonPenalty;
	for(int i=DATA::IntronF; i<=DATA::IntronR; i++)  // Intron(2)
	  d->contents[i] -= intronPenalty; 
	for(int i=DATA::UTR5F; i<=DATA::IntronUTRR; i++) // UTR(4)+IntronUTR(2)
	  d->contents[i] -= UTRPenalty;
 	if (pos == vFin[index]) index++;
      }
  }
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorRepeat :: Plot(DNASeq *TheSeq)
{
  int i;
  
  for (i =0; i < (int)vDeb.size(); i++)
    PlotRepeat(vDeb[i],vFin[i]);
}

// ------------------
//  Post analyse
// ------------------
void SensorRepeat :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}
