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
// $Id: Sensor.Riken.cc,v 1.30 2010-01-25 17:01:23 sallet Exp $
// ------------------------------------------------------------------
// File:     Sensor.Riken.cc
// Contents: Sensor Riken
// RAFL: Riken Arabidopsis Full Length cDNA
// ------------------------------------------------------------------

#include "Sensor.Riken.h"

extern Parameters PAR;


bool Before(const RAFLgene A, const RAFLgene B)
{ return (A.deb < B.deb); };


// ************
// * RAFLgene *     // RAFL: Riken Arabidopsis Full Length cDNA 
// ************
// ----------------------
// Default constructor.
// ----------------------
RAFLgene :: RAFLgene () { deb = fin = sens = 0; ID[0] = '0'; }

// ----------------------
// Default destructor.
// ----------------------
RAFLgene :: ~RAFLgene () {}



/*************************************************************
 **                        SensorRiken                      **
 *************************************************************/
// ----------------------
// Default constructor.
// ----------------------
// ------------->> WARNING: VIRER CE RAFL_A_TRAITER ! << -------------------
// (pour le -pd entre autres)
SensorRiken :: SensorRiken (int n, DNASeq *X) : Sensor(n)
{
  FILE *fRAFL;
  int i, j;
  int beg5, end5;
  int beg3, end3;
  char name[FILENAME_MAX+1];
  char tempname[FILENAME_MAX+1];
  RAFLgene tmp;
  std::vector <RAFLgene> RAFLtmp;

  type = Type_Content;

  StrandRespect = PAR.getI("Riken.StrandRespect");
  MAX_RIKEN_LENGTH = PAR.getI("Riken.Max_riken_length"); // default = 60000;
  MAX_RIKEN_EST_LENGTH = PAR.getI("Riken.Max_riken_est_length"); //default  = 3000;
  MIN_RIKEN_LENGTH = PAR.getI("Riken.Min_riken_length"); //default  = 120; // 2* riken overlap (60)
  MIN_RIKEN_EST_LENGTH = PAR.getI("Riken.Min_riken_est_length"); //default  = 10;
  MAX_OVERLAP = PAR.getI("Riken.Max_overlap"); //default  = 60
  MIN_EST_DIFF = PAR.getI("Riken.Min_est_diff") ; //def=100

  fprintf(stderr, "Reading RAFL gene............");
  fflush(stderr);
  strcpy(tempname, PAR.getC("fstname"));
  strcat(tempname, ".riken");
  inputFormat_ = to_string(PAR.getC("Riken.format", GetNumber(),1));
  if ( inputFormat_ == "GFF3" )
  {

    strcat(tempname,".gff3");
    GeneFeatureSet * geneFeatureSet = new GeneFeatureSet (tempname);
    readRikenGff3 (*geneFeatureSet, X, RAFLtmp);
    delete geneFeatureSet;
  }
  else
  {
    readRiken (tempname, X, RAFLtmp);
  }
  
  
  
  
  fprintf(stderr, "%d RAFL EST pairs read, ", RAFLtmp.size());
  
  sort(RAFLtmp.begin(), RAFLtmp.end(), Before);
  
  int RAFLindice=0;
  int RAFLtmpindice=0;
  
  //    for (RAFLtmpindice=0; RAFLtmpindice< (int)RAFLtmp.size(); RAFLtmpindice++) {
  //     fprintf(stderr, "\nRAFLtmp[%d]: deb=%d fin=%d sens=%d ID=%s\n",RAFLtmpindice,RAFLtmp[RAFLtmpindice].deb,RAFLtmp[RAFLtmpindice].fin,RAFLtmp[RAFLtmpindice].sens,RAFLtmp[RAFLtmpindice].ID);
  //  }
  
  if ((int)RAFLtmp.size()>0) RAFL.push_back(RAFLtmp[RAFLtmpindice]); // on retient le premier
  
  for (RAFLtmpindice=1; RAFLtmpindice<(int)RAFLtmp.size(); RAFLtmpindice++) { 
      // on analyse chaque tmp que l'on compare avec le dernier rafl retenu
    
    if (RAFL[RAFLindice].fin - RAFLtmp[RAFLtmpindice].deb < -3) { // pas d'overlap
      RAFL.push_back(RAFLtmp[RAFLtmpindice]);
      RAFLindice++;
      continue; // tt va bien, on passe au tmp suivant
    }
    if (RAFL[RAFLindice].fin - RAFLtmp[RAFLtmpindice].deb >= MAX_OVERLAP) { // cas de grand overlap
      fprintf(stderr,"fusion %s-%s...",RAFL[RAFLindice].ID,RAFLtmp[RAFLtmpindice].ID);
      if ((RAFL[RAFLindice].sens!=0) && 
	  (RAFLtmp[RAFLtmpindice].sens!=0) && 
	  (RAFL[RAFLindice].sens!=RAFLtmp[RAFLtmpindice].sens) &&
	  (RAFL[RAFLindice].sens!=2) ) { // orientations contraires
	fprintf(stderr, "\nWARNING: RAFL data inconsistent : contig containing rikens in different orientations (%s and %s)\n",RAFL[RAFLindice].ID,RAFLtmp[RAFLtmpindice].ID);
	RAFL[RAFLindice].sens = 2;
      }
      // grand overlap, fusion des deux riken (meme gene)
      RAFL[RAFLindice].deb = Min( RAFL[RAFLindice].deb , RAFLtmp[RAFLtmpindice].deb );
      RAFL[RAFLindice].fin = Max( RAFL[RAFLindice].fin , RAFLtmp[RAFLtmpindice].fin );
      if (RAFL[RAFLindice].sens==0) 
	RAFL[RAFLindice].sens = RAFLtmp[RAFLtmpindice].sens;
      continue;
    
    } else { 
      // little overlap, real overlapping genes
      int const HALF_HOLE_SIZE = 3;
      
      fprintf(stderr,"%s-%s overlapping...",RAFL[RAFLindice].ID,RAFLtmp[RAFLtmpindice].ID);
      RAFL[RAFLindice].fin = ((RAFLtmp[RAFLtmpindice].deb+RAFL[RAFLindice].fin)/2)-HALF_HOLE_SIZE;
      i=RAFLtmpindice; // on raccourcit le debut de tous les petits chevauchant:
      
      while ( (i<(int)RAFLtmp.size()) && (RAFL[RAFLindice].fin >= RAFLtmp[i].deb-2*HALF_HOLE_SIZE) ) {
	RAFLtmp[i].deb= RAFL[RAFLindice].fin+2*HALF_HOLE_SIZE;
	i++;
      }
      
      RAFL.push_back(RAFLtmp[RAFLtmpindice]);
      RAFLindice++;
      continue;
    }
  }
  
  fprintf(stderr,"resulting %d\n",RAFL.size());
  fflush(stderr);

  //    for (RAFLtmpindice=0; RAFLtmpindice< (int)RAFL.size(); RAFLtmpindice++) {
  //       fprintf(stderr, "\nRAFL[%d]: deb=%d fin=%d sens=%d ID=%s\n",RAFLtmpindice,RAFL[RAFLtmpindice].deb,RAFL[RAFLtmpindice].fin,RAFL[RAFLtmpindice].sens,RAFL[RAFLtmpindice].ID);
  //     }
}

void SensorRiken :: readRiken (char name[FILENAME_MAX+1],DNASeq *X, std::vector <RAFLgene> & RAFLtmp )
{
  FILE *fRAFL;
  int  j;
  int beg5, end5;
  int beg3, end3;
  RAFLgene tmp;
  fRAFL = FileOpen(NULL, name, "r");
  while ((j=fscanf (fRAFL,"%d %d %*s %d %d %*s %s",
	  &beg5, &end5, &beg3, &end3, name)) == 5) {
	    tmp.deb = Min(beg3, Min(end3, Min(beg5, end5)));
	    tmp.fin = Max(beg3, Max(end3, Max(beg5, end5)));
	    strcpy(tmp.ID, name);
	    tmp.sens = (((beg5+end5) <= (beg3+end3)) ? 1 : -1);
	    if (StrandRespect==0) tmp.sens=0;
    
 
	    if ( checkRAFL(tmp,beg5,end5,beg3,end3,X->SeqLen) == 0 )
	      continue;
	    RAFLtmp.push_back(tmp);
	  }
  
	  if (j != EOF) {
	    fprintf(stderr, "Incorrect RAFL file\n");
	    exit(2);
	  }
  
	  fclose(fRAFL);
}


void SensorRiken :: readRikenGff3 (GeneFeatureSet & geneFeatureSet ,DNASeq *X, std::vector <RAFLgene> & RAFLtmp)
{
  int  j;
  int beg5, end5;
  int beg3, end3;
  RAFLgene tmp;
  map < string, vector <GeneFeature *> >::iterator it= geneFeatureSet.getIteratorParentToChildren();
  //vector<GeneFeature *>::iterator it = geneFeatureSet.getIterator();
  int nbGeneFeature=geneFeatureSet.getNbParentFeature();
  int i=0;

  for (i = 0; i < nbGeneFeature ; i++, it++ )
  {
    GeneFeature * parent = geneFeatureSet.getGeneFeature( (*it).first );
    if ((*it).second.size() !=2 )
    {
      fprintf(stderr, "\nWARNING: Parent in gff3 hasn't 2 children (5' and 3' exon)\n");
    }
    GeneFeature * child5 = ((*it).second[0]);
    GeneFeature * child3 = ((*it).second[1]);
    beg5 = child5->getLocus()->getStart();
    end5 = child5->getLocus()->getEnd();
    if ( child5->getLocus()->getStrand() == '-' )
    {
      int swap=end5 ;
      end5=beg5;
      beg5 = swap;
    }
    beg3 = child3->getLocus()->getStart();
    end3 = child3->getLocus()->getEnd();
    if ( child3->getLocus()->getStrand() == '-' )
    {
      int swap=end3 ;
      end3=beg3;
      beg3 = swap;
    }
    tmp.deb = Min(beg3, Min(end3, Min(beg5, end5)));
    tmp.fin = Max(beg3, Max(end3, Max(beg5, end5)));
    strcpy(tmp.ID, parent->getId().c_str());
    tmp.sens = (((beg5+end5) <= (beg3+end3)) ? 1 : -1);

    if ( checkRAFL(tmp,beg5,end5,beg3,end3, X->SeqLen) == 0 )
      continue;
    
    RAFLtmp.push_back(tmp);
  }
}

int SensorRiken :: checkRAFL ( RAFLgene & tmp, int beg5, int end5, int beg3, int end3, int len )
{
  int res=1;
  char name[FILENAME_MAX+1];
  strcpy ( name, tmp.ID);
  if (StrandRespect==0) tmp.sens=0;
    
    // si le gene est trop court, on ne peut pas connaitre le sens !
  if ( (StrandRespect) && (abs((beg5+end5) - (beg3+end3)) < MIN_EST_DIFF) ) {
    tmp.sens = 0;
    fprintf(stderr, "\nWARNING: Check RAFL data: Riken %s has no significative orientation\n", name);
  }
    
  if (abs(tmp.deb-tmp.fin) > MAX_RIKEN_LENGTH) {
    fprintf(stderr, "\nWARNING: Check RAFL data: Riken %s rejected, too long transcript (%d kb)\n", name, abs(tmp.deb-tmp.fin)/1000);
    res=0;
  }
  if (abs(beg5-end5) > MAX_RIKEN_EST_LENGTH) {
    fprintf(stderr, "\nWARNING: Check RAFL data: Riken %s rejected, 5' EST mapped to %d bp\n",
	    name, abs(beg5-end5));
    res=0;
  }
  if (abs(beg3-end3) > MAX_RIKEN_EST_LENGTH) {
    fprintf(stderr, "\nWARNING: Check RAFL data: Riken %s rejected, 3'' EST mapped to %d bp\n",
	    name, abs(beg3-end3));
    res=0;
  }
  if (abs(tmp.deb-tmp.fin) < MIN_RIKEN_LENGTH) {
    fprintf(stderr, "\nWARNING: Check RAFL data: Riken %s rejected, too short transcript (%d bp)\n",
	    name, abs(tmp.deb-tmp.fin));
    res=0;
  }
  if (abs(beg5-end5) < MIN_RIKEN_EST_LENGTH) {
    fprintf(stderr, "\nWARNING: Check RAFL data: Riken %s rejected, 5' EST mapped to %d bp\n",
	    name, abs(beg5-end5));
    res=0;
  }
  if (abs(beg3-end3) < MIN_RIKEN_EST_LENGTH) {
    fprintf(stderr, "\nWARNING: Check RAFL data: Riken %s rejected, 3'' EST mapped to %d bp\n",
	    name, abs(beg3-end3));
    res=0;
  }
  if (tmp.deb <= 0  || tmp.fin > len) {
    fprintf(stderr, "\nWARNING: Check RAFL data: Riken %s rejected, coordinate(s) %d-%d out of range\n",
	    name, tmp.deb, tmp.fin);
    res=0;
  }
  return res;
}


// ----------------------
//  Default destructor.
// ----------------------
SensorRiken :: ~SensorRiken ()
{
  RAFL.clear();
}

// ----------------------
//  Init Riken.
// ----------------------
void SensorRiken :: Init (DNASeq *X)
{
  RAFLPenalty = PAR.getD("Riken.RAFLPenalty*"); //default  = -120;

  RAFLpos = RAFLindex = 0;

  if (PAR.getI("Output.graph")) Plot(X);
}

// -------------------------
//  GiveInfo Content Riken.
// -------------------------
void SensorRiken :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  // Calcul de la position par rapport aux genes RAFL (Riken Ara.Full-Length)
  // valeurs de RAFLpos:
  // 0-> en dehors,
  // 1-> frontiere 5prime (intergenique ou UTR5  obligatoire),
  // 2-> frontiere 3prime (intergenique ou UTR3  obligatoire),
  // 3-> dedans(penalisation IG)
  // 4-> bordure min (StrandRespect==0)
  // 5-> bordure max  (StrandRespect==0)

  // Out of seq range, do nothing
  if (pos<0 || pos>=X->SeqLen) return;
  
  if(!RAFL.empty()) {
    RAFLpos = 0;
    
    while (pos < RAFL[RAFLindex].deb-2) {
      if (RAFLindex-1 >= 0) RAFLindex--; 
      else return;
    }
    
    while (pos > RAFL[RAFLindex].fin) {
      // si on depasse le RAFL, on prend l'eventuel prochain
      if (RAFLindex+1 < (int)RAFL.size()) RAFLindex++; 
      else return;
    }
    
    // dans le Riken:
    if (pos >= RAFL[RAFLindex].deb-2) {
      // bordure min:
      if (pos == RAFL[RAFLindex].deb-2) {
	if ((RAFL[RAFLindex].sens== 1) || (RAFL[RAFLindex].sens== -1))
	  RAFLpos = ( (RAFL[RAFLindex].sens == 1) ? 1 : 2);
	else 
	  RAFLpos=4;
      }
      else {
	if (pos < RAFL[RAFLindex].fin) RAFLpos = 3;
	else {
	  if (pos == RAFL[RAFLindex].fin) { // frontiere max:
	    if ((RAFL[RAFLindex].sens== 1) || (RAFL[RAFLindex].sens== -1))
	      RAFLpos = ( (RAFL[RAFLindex].sens == -1) ? 1 : 2);
	    else 
	      RAFLpos=5;
	  }
	}
      }
    }

    //    fprintf(stdout,"i=%d, raflpos= %d\n",pos,RAFLpos);

    // Exons forward
    for(int i = DATA::ExonF1; i<= DATA::ExonF3; i++)
      if ((RAFLpos == 1) || (RAFLpos == 2) || (RAFLpos >= 4) ||
	  ((RAFLpos == 3) && (RAFL[RAFLindex].sens == -1)) )
	d->contents[i] += RAFLPenalty;  

    // Exons Reverse
    for(int i=DATA::ExonR1; i<=DATA::ExonR3; i++)
      if ((RAFLpos == 1) || (RAFLpos == 2) || (RAFLpos >= 4) ||
	  ((RAFLpos == 3) && (RAFL[RAFLindex].sens == 1)) )
	d->contents[i] += RAFLPenalty; 

    // IntronF
    if ( (RAFLpos == 1) || (RAFLpos == 2) || (RAFLpos >= 4) ||
	 ( (RAFLpos == 3) && (RAFL[RAFLindex].sens == -1)) )
      d->contents[DATA::IntronF] += RAFLPenalty;

    //IntronR
    if ( (RAFLpos == 1) || (RAFLpos == 2) || (RAFLpos >= 4) ||
	 ( (RAFLpos == 3) && (RAFL[RAFLindex].sens == 1)) )
      d->contents[DATA::IntronR] += RAFLPenalty;     

    //InterG
    if (RAFLpos == 3) {
      d->contents[DATA::InterG] += RAFLPenalty;
      d->contents[DATA::RNAF] += RAFLPenalty; 
      d->contents[DATA::RNAR] += RAFLPenalty; 
      
      
      // Warning SIGNAL START MODIFICATION //
      d->sig[DATA::Start].weight[Signal::ForwardNo] = 0.0;
      d->sig[DATA::Start].weight[Signal::ReverseNo] = 0.0;
    }

    //UTR5'F
    if ( (RAFLpos == 2) || (RAFLpos == 5) || 
	 ((RAFL[RAFLindex].sens == -1) && (RAFLpos >= 1) ) )
      d->contents[DATA::UTR5F] += RAFLPenalty;     
    
    //UTR5'R
    if ( (RAFLpos == 2) || (RAFLpos == 4) || 
	 ((RAFL[RAFLindex].sens == 1) && (RAFLpos >= 1) ) )
      d->contents[DATA::UTR5R]+= RAFLPenalty;     
    
    //UTR3'F
    if ( (RAFLpos == 1) || (RAFLpos == 4) || 
	 ((RAFL[RAFLindex].sens == -1) && (RAFLpos >= 1) ) )
      d->contents[DATA::UTR3F]+= RAFLPenalty; 
    
    //UTR3'R
    if ( (RAFLpos == 1) || (RAFLpos == 5) || 
	 ((RAFL[RAFLindex].sens == 1) && (RAFLpos >= 1) ) )
      d->contents[DATA::UTR3R]+= RAFLPenalty;
  }
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorRiken :: Plot(DNASeq *TheSeq)
{
  int plotcolor=0;

  for (int j =0; j<(int)RAFL.size() ;j++) {
    plotcolor= ( (RAFL[j].sens==2) ? 2 : 2-RAFL[j].sens);
    PlotStartReg(RAFL[j].deb, 0, plotcolor);
    PlotEndReg(RAFL[j].fin,0,plotcolor);
  }
}

// ------------------
//  Post analyse
// ------------------
void SensorRiken :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}
