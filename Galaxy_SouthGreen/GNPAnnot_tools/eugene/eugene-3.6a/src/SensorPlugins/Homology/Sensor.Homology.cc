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
// $Id: Sensor.Homology.cc,v 1.25 2009-01-20 09:40:47 sallet Exp $
// ------------------------------------------------------------------
// File:     Sensor.Homology.cc
// Contents: Sensor Homology  
// Comment:  Travail sur la prise en compte des informations d'homologies.
// A partir de hits tblastx, on tente d'aider la prediction.
// Nouveau parametre TblastxP dans le fichier .par
// FORMATS necessaire:
// fichier matrice proteique: fichier BLOSUM ou PAM classique
// fichier .tblastx : comme les .blast, mais derniere colonne= seq prot. du
// hit subject.
// Ex:
// 802 864 104 0.0 +1 AB015498_354_416 354 416 PQGQTPLFPRIFGHEAGG*EKVGLWRVLEKV*LILHQEI
// ------------------------------------------------------------------

#include "Sensor.Homology.h"

extern Parameters PAR;

/*************************************************************
 **                      SensorHomology                     **
 *************************************************************/

// ----------------------
// fonction temporaire
// etude de la prise en compte possible des tblastx
// ----------------------
double SensorHomology :: tblastxupdate (int hitnb, double hitscore, double pen, double base) {
  // return (pen); // test avec score constant pour chq nt d'un hit
  // test d'integration de la donnee nbre de hits:
  //	return (hitnb);    
  //	return (hitnb+pen);
  //	return (hitnb*base +pen);
  // test d'integration de la donnee score du hit:
  //	return (hitscore);
  //	return (hitscore +pen);
  return (hitscore*base +pen);
  //	return (hitscore*base);
}

// --------------------------------------------------
// return the maximum number of hits for this contig 
// (overlaping hits in the same frame) 
// --------------------------------------------------
int SensorHomology :: hitsmaxnumber (int position, int frame, int maxlen) {
  int maxnumber=0;
  for (int i=position; i<=maxlen ; i++) {
    if (TblastxNumber[frame][i] == 0) return maxnumber;
    maxnumber = Max(maxnumber,(int)TblastxNumber[frame][i]);
  }
  return maxnumber; // case of hit reaching the end of the seq
}

// ----------------------
//  Default constructor.
// ----------------------
SensorHomology :: SensorHomology(int n, DNASeq *X) : Sensor(n)
{
  type = Type_Content;

  int i, j;
  int Len = X->SeqLen; 
  FILE *protmatfile;
  ProtMat* PROTMAT;
  char tempname[FILENAME_MAX+1];
  int nmax;
  int contigbeg;
  char *tmpdir = new char[FILENAME_MAX+1];
  
  fileExt       = PAR.getC("Homology.FileExtension", GetNumber(),1);
  int MaxHitLen = PAR.getD("Homology.MaxHitLen", GetNumber(),1);
  inputFormat_  = to_string(PAR.getC("Homology.format", GetNumber(),1));

  TblastxNumber = new short int*[6];
  TblastxScore = new float*[6];
  for (j = 0; j < 6; j++) {
    TblastxNumber[j]= new short int[Len+1];
    TblastxScore[j] = new float[Len+1];
    for (i = 0; i<= Len; i++) {
      TblastxNumber[j][i]=0;
      TblastxScore[j][i]=0.0;
    }
  }

  strcpy(tmpdir, PAR.getC("eugene_dir"));
  strcat(tmpdir, MODELS_DIR);
  strcpy(tempname,PAR.getC("Homology.protmatname", GetNumber()));
  // Lecture de la matrice proteique ("BLOSUM62" par defaut)
  protmatfile=FileOpen(tmpdir, tempname, "rt");
  if (protmatfile == NULL) {
    fprintf (stderr,"\ncannot open protein matrix file %s\n",tempname);
    exit (2);
  }
  // chargement classe protmat
  fprintf(stderr,"Reading protmat file.........");
   fflush(stderr);
  if (fichier2protmat(protmatfile , PROTMAT)) {
    fprintf(stderr,"error when reading protmat file %s\n",tempname);
    exit(2);
  }
  fclose(protmatfile);
  fprintf(stderr,"done\n");
  fflush(stderr);
  
  fprintf(stderr,"Reading tblastx data.........");
  fflush(stderr);
  strcpy(tempname, PAR.getC("fstname"));
    // concatenation of file extension
  strcat(tempname, fileExt);


  if ( inputFormat_ == "GFF3" )
  {
    strcat(tempname,".gff3");
    GeneFeatureSet * geneFeatureSet = new GeneFeatureSet (tempname);
    ReadHomologyGff3(*geneFeatureSet, X, MaxHitLen, PROTMAT);
    delete geneFeatureSet;
  }
  else
  {
    ReadHomology(tempname, X, MaxHitLen, PROTMAT);
  }
  fprintf(stderr, "done\n");
  fflush(stderr);

  // weighting by the ratio of hits number.
  for (j=0; j<6 ;j++) {
    contigbeg=1;
    nmax=0;
    for (i=0 ; i<=Len ; i++) {
      if (TblastxNumber[j][i] > 0) {
	if (contigbeg == 1) {
	  nmax = hitsmaxnumber(i,j,Len);
	}
	contigbeg=0;
	TblastxScore[j][i] /= nmax;
      }
      else {
	nmax=0;
	contigbeg=1;
      }
    }
  }
  delete [] tmpdir;
}

// ----------------------
//  Default destructor.
// ----------------------
SensorHomology :: ~SensorHomology ()
{
  for(int i=0;i<6;i++){
    delete TblastxNumber[i];
    delete TblastxScore[i];
  }
  delete [] TblastxNumber;
  delete [] TblastxScore;
}



void SensorHomology:: ReadHomology (char name[FILENAME_MAX+1],DNASeq *X,  const int MaxHitLen,ProtMat* PROTMAT )
{
  int i, j;
  
  int Len = X->SeqLen; 
  FILE *ftblastx;
  FILE *protmatfile;
  int deb,fin,phase,ProtDeb,ProtFin,sens;
  int score = 0;
  double bits;

  char tampon;
  char* paire= new char[3];
  paire[0]=paire[1]= '0';
  paire[2]='\0';
  char tempname[FILENAME_MAX+1];
  int nmax;
  int contigbeg;
  char *tmpdir = new char[FILENAME_MAX+1];

  ftblastx = FileOpen(NULL,name, "r", PAR.getI("EuGene.sloppy"));
  if (ftblastx) {
  
    while (!feof(ftblastx)) {
      if (fscanf(ftblastx,"%d %d %lf %*s %d %*s %d %d",
	  &deb, &fin, &bits, &phase, &ProtDeb, &ProtFin) != 6) {
	    if (feof(ftblastx)) break;
	    fprintf(stderr,"error(1) in tblastxfile format %s\n",name);
	    exit(2);
	  }
	  if (abs(fin-deb) > MaxHitLen) {
	    fprintf(stderr,"Similarity of extreme length rejected. Check tblastx file %s\n",name);
	    fscanf(ftblastx,"%*s");
	    continue;
	  }
	  sens= ( (phase < 0) ? -1 : 1 );
	  phase = ph06(phase);
    
	  j=deb;
	  deb= Min (deb,fin);
	  fin= Max (j,fin);
    
          // lecture de la suite (sequence hit subject)
	  tampon=fgetc(ftblastx);
	  
          // en cas de plsrs separateurs avant la derniere colonne
	  while (isspace(tampon)) tampon=fgetc(ftblastx);
	  for (i = deb-1; i < fin; i++)  {
	    if ( (feof(ftblastx)) ||
			 ( (isspace(tampon)) && (i != (fin)) ) ) {
	      fprintf(stderr,"error(2) in tblastx format, file %s\n",name);
	      exit(2);
			 }
			 if ( (i-deb+1)%3 == 0 ) {
			   if (i >deb-1) tampon=fgetc(ftblastx);
			   paire[1]=tampon;
			   paire[0]= ( (sens>0) ? X->AA(i,0) : X->AA(deb+fin-i-1,1) );
			   score= PROTMAT->VAL[PROTMAT->mot2indice(paire)];
			 }
			 TblastxNumber[phase][i] ++;
      
			 TblastxScore[phase][i]  += score;
	  }
    }
    fclose(ftblastx);
    }
    
}

void SensorHomology ::ReadHomologyGff3(GeneFeatureSet & geneFeatureSet ,DNASeq *X,  const int MaxHitLen, ProtMat* PROTMAT)
{

  int Len = X->SeqLen; 
  int deb,fin,phase,ProtDeb,ProtFin,sens;
  int score = 0;
  double bits= 0;
  string feature;
  char tampon,strand;
  char* paire= new char[3];
  paire[0]=paire[1]= '0';
  paire[2]='\0';

  vector<GeneFeature *>::iterator it = geneFeatureSet.getIterator();
  int nbGeneFeature=geneFeatureSet.getNbFeature();
  int i=0;

  for ( i=0 ; i < nbGeneFeature ; i++, it++ )
  {
	  feature  = (*it)->getType();
	  deb      = (*it)->getLocus()->getStart();
	  fin      = (*it)->getLocus()->getEnd();
	  strand   = (*it)->getLocus()->getStrand();
	  if ( ! (*it)->hasTarget() )
	  {
	    fprintf( stderr, "Warn : skipped : %s \n", ((*it)->getId()).c_str() );
	    fflush(stderr);
	    continue;
	  }
	  
	  bits     = (*it)->getAttributes()->getTarget()->getScoreHit();
	  phase    = (*it)->getAttributes()->getTarget()->getFrameHit();
	  if (phase == 0) //not specified
	  { 
	     phase  = X->Pos2Frame(deb,strand); 
	  }
          else 
	  {
	    int computedFrame=X->Pos2Frame(deb,strand);
	    if (phase != computedFrame) //check frame computrd and read are the same.
	    {
	      fprintf( stderr, "Warn skipped : %s,  computed frame (%d) and input frame (%d) are different \n",(*it)->getId().c_str(), computedFrame, phase);
	      fflush(stderr);
	      continue;
	    }
	  }
	  
	 
	  ProtDeb=(*it)->getAttributes()->getTarget()->getLocus()->getStart();
	  ProtFin=(*it)->getAttributes()->getTarget()->getLocus()->getEnd();
    
	  if (abs(fin-deb) > MaxHitLen) {
	    fprintf(stderr,"Similarity of extreme length rejected. Check tblastx file \n");
	    fflush(stderr);
	    continue;
	  }
	  sens   = ( (strand == '-') ? -1 : 1 );
          //converti la frame blast en piste eugene
	  phase  = ph06(phase);

	  if ( (*it)->getAttributes()->getTarget()->getSequenceData() == "")
	  {
	    fprintf( stderr, "Warn Feature hasn't target sequence, skipped : %s \n", ((*it)->getString()).c_str() );
	    fflush(stderr);
	    continue;
	  }
	  
	  char subjectHit [MAX_GFF_LINE] ;
	  strcpy (subjectHit, ((*it)->getAttributes()->getTarget()->getSequenceData()).c_str() );
	  int j=0;
	   
	  tampon=subjectHit[j];
	  for (int l = deb-1; l < fin; l++)  {
		 if ( (l-deb+1)%3 == 0 ) {
		   if (l >deb-1) 
		       tampon=subjectHit[j];
			   
		   paire[1]=tampon;
		   paire[0]= ( (sens>0) ? X->AA(l,0) : X->AA(deb+fin-l-1,1) );
		   score= PROTMAT->VAL[PROTMAT->mot2indice(paire)];
		   j++;
		 }
		 TblastxNumber[phase][l] ++;
		 TblastxScore[phase][l]  += score;
	  }
    }
}

// ----------------------
//  Init Homol.
// ----------------------
void SensorHomology :: Init (DNASeq *X)
{
  TblastxP= PAR.getD("Homology.TblastxP*", GetNumber());
  TblastxB= PAR.getD("Homology.TblastxB*", GetNumber());

  if (PAR.getI("Output.graph")) Plot(X);
}

// --------------------------
//  GiveInfo Content TBlastX.
// --------------------------
void SensorHomology :: GiveInfo(DNASeq *X, int pos, DATA *d)
{
  int k;
  for (k = 0; k<6; k++) {
    if (TblastxNumber[k][pos]>0)
      d->contents[k] += tblastxupdate(TblastxNumber[k][pos],TblastxScore[k][pos],TblastxP,TblastxB);
  }
}

// -----------------------------------------------
//  Convertit les phases 0-6 en 1 2 3 -1 -2 -3 0.
// -----------------------------------------------
int SensorHomology :: PhaseAdapt(char p)
{
  if (p >= 12) return 0;
  else if (p < 3) return (1+p);
  else if (p < 6) return (2-p);
  else if (p < 9) return (p-2);
  else return (5-p);
}

// -----------------------------------------------
//  Convertit les phases 1 2 3 -1 -2 -3 0 en 0-6.
// -----------------------------------------------
char SensorHomology :: ph06(char p)
{
  if (p == 0) return 6;
  else if (p > 0) return (p-1);
  else return 2-p;   
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorHomology :: Plot(DNASeq *X)
{
  for (int pos =0; pos < X->SeqLen ; pos++){
    for(int j=0;j<6;j++) {
      if (TblastxNumber[j][pos]>0) {
	PlotBarI (pos,PhaseAdapt(j),0.6, TblastxNumber[j][pos], ((TblastxScore[j][pos]>0)? 8- (int) (Min(2,1+(int)TblastxScore[j][pos]/4)) : 8 ) );
      }
    }
  }
}

// --------------------------
//  Affichage des TBlastX.
// --------------------------
void SensorHomology :: PlotAt(int pos)
{
  //  if (PAR.getI("Output.graph")) {
  for(int j=0;j<6;j++) {
    if (TblastxNumber[j][pos]>0) {
      PlotBarI (pos,PhaseAdapt(j),0.6, TblastxNumber[j][pos], ((TblastxScore[j][pos]>0)? 8- (int) (Min(2,1+(int)TblastxScore[j][pos]/4)) : 8 ) );
    }
  }
}

// ------------------
//  Post analyse
// ------------------
void SensorHomology :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}
