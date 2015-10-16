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
// $Id: Sensor.GFF.cc,v 1.13 2009-01-12 14:32:06 sallet Exp $
// ------------------------------------------------------------------
// File:     Sensor.GFF.cc
// Contents: Sensor GFF  
// ------------------------------------------------------------------

#include "Sensor.GFF.h"

/*************************************************************
 **                        GFFObject                        **
 *************************************************************/
// -------------------------
//  Default constructor.
// -------------------------
GFFObject :: GFFObject ()
{
  name[0]        = '\000';
  feature[0]     = '\000';
  start  = end   = 0;
  strand = frame = '?';
  A = T = C = G = N = -1;
  GC = -1.0;
}

// -------------------------
//  Default constructor.
// -------------------------
GFFObject :: GFFObject (char* nam, char* fea, int s, int e, char st, char fr,
			int a, int t, int c, int g, int n, float gc)
{
  strcpy(name,nam);
  strcpy(feature,fea);
  start = s;
  end   = e;
  strand = st;
  frame  = fr;
  A = a;
  T = t;
  C = c;
  G = g;
  N = n;
  GC = gc;
}

// -------------------------
//  Print GFFObject
// -------------------------
void GFFObject :: Print ()
{
  printf("%s\t\t%s\t%d\t%d\t\t%c\t%c\n",
	 name, feature, start, end, strand, frame);
}

// -------------------------
//  Print GFF Header
// -------------------------
void GFFObject :: PrintHeader ()
{
  printf("#S.Name\tSource\tFeature\tStart\tEnd\tScore\tStrand\tFrame\n");
}

// -------------------------
//  Default destructor.
// -------------------------
GFFObject :: ~GFFObject () {}


/*************************************************************
 **                       SensorGFF                         **
 *************************************************************/

extern Parameters PAR;

// ----------------------
// Default constructor.
// ----------------------
SensorGFF :: SensorGFF (int n, DNASeq *X) : Sensor(n)
{
  type = Type_None;
}

// ----------------------
//  Default destructor.
// ----------------------
SensorGFF :: ~SensorGFF ()
{
  // Clear the data structures
  gffList.clear();

  if (PAR.getI("GFF.PostProcess")) fclose(ppfile);
}

// ----------------------
//  Init.
// ----------------------
void SensorGFF :: Init (DNASeq *X)
{
  char tempname[FILENAME_MAX+1];

  // Clear the data structure
  gffList.clear();

  fprintf(stderr, "Reading GFF file.............");
  fflush(stderr);

  inputFormat_ = to_string(PAR.getC("GFF.format", GetNumber(),1));
  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".gff");

   if ( inputFormat_ == "GFF3" )
  {
    strcat(tempname,".gff3");

    GeneFeatureSet * geneFeatureSet = new GeneFeatureSet (tempname) ;
    ReadGFF3(*geneFeatureSet, X);
    delete geneFeatureSet;
  }
  else
  {
    ReadGFF(tempname, X->SeqLen);
  }
  fprintf(stderr, "done\n");
  fflush(stderr);

  if (!IsInitialized) {
    if (PAR.getI("GFF.PostProcess")  &&  gffList[0]->GC != -1.0) {
      if (!(ppfile = FileOpen(NULL, "statSeqs.txt", "w"))) {
	fprintf(stderr, "cannot open GFF file %s\n", "statSeqs.txt");
	exit(2);
      }
      if(strlen(gffList[0]->name) > 7)
	fprintf(ppfile,"SeqName\t\tE/I\tNb\tLmin\tLmax\tGCmin\tGCmax\tV_GC\n");
      else
	fprintf(ppfile,"SeqName\tE/I\tNb\tLmin\tLmax\tGCmin\tGCmax\tV_GC\n");
    }
    IsInitialized = true;
  }

  if (PAR.getI("Output.graph")) Plot(X);

  //for (int i=0; i<(int)gffList.size(); i++) {
  //if(i==0)
  //gffList[i]->PrintHeader();
  //gffList[i]->Print();
  //fflush(stdout);
  //}
}

// --------------------------
//  Read start forward file.
// --------------------------
void SensorGFF :: ReadGFF (char name[FILENAME_MAX+1], int seqlen)
{
  FILE  *fp;
  char  line[MAX_LINE];
  int   i;
  char  *seqname = (char *)malloc(FILENAME_MAX*sizeof(char));
  char  *feature = (char *)malloc(FILENAME_MAX*sizeof(char));
  int   start, end;
  char  strand, frame, phase[2];
  int   a  = -1, t = -1, c = -1, g = -1, n = -1;
  float gc = -1.0;

  if (!(fp = FileOpen(NULL, name, "r"))) {
    fprintf(stderr, "cannot open GFF file %s\n", name);
    exit(2);
  }
  
  int j=0;
  while(fgets (line, MAX_LINE, fp) != NULL) {
    j++;
    if (line[0] != '#') {
      i = sscanf(line,"%s %*s %s %d %d %*s %c %c %d %d %d %d %d %f",
		 seqname, feature, &start, &end, &strand, phase,
		 &a, &t, &c, &g, &n, &gc);
	
      if (i < 6) {
	if (i==-1) {
	  if(j==1)
	    fprintf(stderr,"WARNING: empty GFF file !...");
	}
	else {
	  fprintf(stderr, "Error in GFF file %s, line %d.\n", name, j);
	  exit(2);
	}
      }
      else {
	frame = '.';
	if (phase[0] != '.') {
	  int f = -1;
	  if (strand == '+')      { f = (start  - 1)   % 3; }
	  else if (strand == '-') { f = (seqlen - end) % 3; }
	  char v[2];
	  if (f != -1) { sprintf(v, "%d" , ((f + atoi(phase)) % 3)); }
	  frame = v[0];
	}
	gffList.push_back(new GFFObject(seqname, feature, start, end,
					strand,  frame,
					a, t, c, g, n, gc));

      }
    }
  }
  fclose(fp);
}

// --------------------------
//  Read gff3 file.
// --------------------------

void SensorGFF :: ReadGFF3 (GeneFeatureSet & geneFeatureSet , DNASeq *X )
{

  char  *seqname = (char *)malloc(FILENAME_MAX*sizeof(char));
  char  *feature = (char *)malloc(FILENAME_MAX*sizeof(char));
  int   start, end;
  char  strand, frame, phase[2];
  int   a  = -1, t = -1, c = -1, g = -1, n = -1;
  float gc = -1.0;
  string idSo;
  int seqlen=X->SeqLen;

  vector<GeneFeature *>::iterator it = geneFeatureSet.getIterator();
  int nbGeneFeature=geneFeatureSet.getNbFeature();
  for ( int j=0 ; j < nbGeneFeature ; j++, it++ )
  {
     strcpy ( seqname, (*it)->getSeqId().c_str() );
     strcpy (feature, (*it)->getType().c_str());
     start = (*it)->getLocus()->getStart();
     end = (*it)->getLocus()->getEnd();
     strand = (*it)->getLocus()->getStrand();
     // En natif lecture de : a, t, c, g, n, gc ???
     
     idSo=(*it)->getType();
     if ( idSo.find("SO:") == string::npos )
     {
      string tmp=GeneFeatureSet::soTerms_->getIdFromName(idSo);
      idSo=tmp;
     }
     if ( idSo!="SO:0000316" && idSo!="SO:0000204" && idSo!="SO:0000205" )
     {
	continue;
     }

     phase[0] = (*it)->getPhase() ;
     frame = '.';
     if (phase[0]  != '.') //selection des CDS == exon
     {
	 int f = -1;
	  if (strand == '+')      { f = (start  - 1)   % 3; }
	  else if (strand == '-') { f = (seqlen - end) % 3; }
	  char v[2];
	  if (f != -1) { sprintf(v, "%d", ((f + atoi(phase)) % 3)); }
	  frame = v[0];
     	  
     }
     gffList.push_back(new GFFObject(seqname, feature, start, end,
					strand,  frame,
					a, t, c, g, n, gc));
  }
}

// ------------------------
//  GiveInfo.
// ------------------------
void SensorGFF :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
}

// ----------------------------
//  Plot Sensor information.
// ----------------------------
void SensorGFF :: Plot(DNASeq *X)
{
  const int PredWidth = 2;
  char stGC[10];

  for(int i=0; i<(int)gffList.size(); i++) {
    if (gffList[i]->frame != '.') {
      int f = ((int)gffList[i]->frame - 48) + 1;
      if (gffList[i]->strand == '-')
	f *= -1;

      for (int j=gffList[i]->start; j<=gffList[i]->end; j++)
	PlotBarI(j, f, 0.27, PredWidth, 9);

      if (gffList[0]->GC != -1) {
	sprintf(stGC, "%.2f", gffList[i]->GC);
	//strcat(stGC, " %");
	//PlotString((gffList[i]->end+gffList[i]->start)/2, f, 0.21, stGC, 9);
	PlotString((gffList[i]->end+gffList[i]->start)/2, f, -0.02, stGC, 9);
      }
      
      if (gffList[i]->strand == '.') {
	f *= -1;
	for (int j=gffList[i]->start; j<=gffList[i]->end; j++)
	  PlotBarI(j, f, 0.27, PredWidth, 9);
	if (gffList[0]->GC != -1) {
	  PlotString((gffList[i]->end+gffList[i]->start)/2, f, -0.02, stGC, 9);
	}
      }
    }
    else {
      for (int j=gffList[i]->start; j<=gffList[i]->end; j++)
	PlotBarI(j, 0, 0.27, PredWidth, 9);
    }
  }
}

// ------------------
//  Post analyse.
// ------------------
void SensorGFF :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
  if (!PAR.getI("GFF.PostProcess")) return;
  
  if (gffList[0]->A  != -1  &&  gffList[0]->T != -1  &&
      gffList[0]->C  != -1  &&  gffList[0]->G != -1  &&
      gffList[0]->GC != -1) {
    int   LminI  = 100000, LmaxI  = 0;
    int   LminE  = 100000, LmaxE  = 0;
    float GCminI = 101.0 , GCmaxI = 0.0;
    float GCminE = 101.0 , GCmaxE = 0.0;
    int   nbI = 0, nbE = 0;
  
    for(int i=0; i<(int)gffList.size(); i++) {
      int len = gffList[i]->A + gffList[i]->T + gffList[i]->C + gffList[i]->G;
            
      if (gffList[i]->feature[0] == 'E') {
	nbE++;
        if (len > LmaxE) LmaxE = len;
	if (len < LminE) LminE = len;
	if (len > 70 ) {
	  if (gffList[i]->GC > GCmaxE) GCmaxE = gffList[i]->GC;
	  if (gffList[i]->GC < GCminE) GCminE = gffList[i]->GC;
	}
      }
      else if (gffList[i]->feature[0] == 'I') {
	nbI++;
	if (len > LmaxI) LmaxI = len;
	if (len < LminI) LminI = len;
	if (len > 70) {
	  if (gffList[i]->GC > GCmaxI) GCmaxI = gffList[i]->GC;
	  if (gffList[i]->GC < GCminI) GCminI = gffList[i]->GC;
	}
      }
    }
    
    if (nbE == 0)
      fprintf(ppfile,"%s\tExon\t%d\t-\t-\t-\t-\t-\n", gffList[0]->name,nbE);
    else if (GCminE == 101.0  &&  GCmaxE == 0.0)
      fprintf(ppfile,"%s\tExon\t%d\t%d\t%d\tLen<70\tLen<70\tLen<70\n",
	      gffList[0]->name,nbE,LminE,LmaxE);
    else fprintf(ppfile,"%s\tExon\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\n",
		 gffList[0]->name,nbE,LminE,LmaxE,GCminE,GCmaxE,GCmaxE-GCminE);
    
    if (nbI == 0)
      fprintf(ppfile,"%s\tIntron\t%d\t-\t-\t-\t-\t-\n\n",gffList[0]->name,nbI);
    else if (GCminI == 101.0  &&  GCmaxI == 0.0)
      fprintf(ppfile,"%s\tIntron\t%d\t%d\t%d\tLen<70\tLen<70\tLen<70\n\n",
	      gffList[0]->name,nbI,LminI,LmaxI);
    else fprintf(ppfile,"%s\tIntron\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\n\n",
		 gffList[0]->name,nbI,LminI,LmaxI,GCminI,GCmaxI,GCmaxI-GCminI);
  }
}
