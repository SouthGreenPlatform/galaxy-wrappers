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
// $Id: Sensor.MarkovProt.cc,v 1.15 2010-01-25 17:01:23 sallet Exp $
// ------------------------------------------------------------------
// File:     Sensor.MarkovProt.cc
// Contents: Sensor MarkovProt
// ------------------------------------------------------------------

#include "Sensor.MarkovProt.h"

extern Parameters PAR;

double SCOREMIN = -5; // for graph option only (plot)

/*************************************************************
 **                        SensorMarkovProt                     **
 *************************************************************/
double Recadre (double val){
  if (val < SCOREMIN) return SCOREMIN;
  return val;
}

// ---------------------------------------------------- 
// Normalisation sur les 6 phases + non codant
// ---------------------------------------------------- 
void AmplifyScore(double Score [], unsigned int normopt)
{
  double Sum = 0.0;
  double Min = -log(0.0);
  int i;

  for  (i = 0;  i < 9;  i ++) 
    if (Min > Score[i]) Min = Score[i];

  switch (normopt) {
  case 0:
    for  (i = 0;  i < 9;  i ++) 
      Score[i] = exp(Min-Score[i]);
    break;

  case 1:
    for  (i = 0;  i < 9;  i ++) {
      Score[i] = exp(Score[i]-Min);
      Sum += Score[i];
    }
    for  (i = 0;  i < 9;  i ++) 
      Score [i] /= Sum;
    break;

  case 2:
    for (i = 0; i < 9; i++) Score[i] = exp(Score[i]-Min);
    for (i = 0; i < 6; i++) {
      Sum += Score[i];
      Score [i] /= (Score[6]+Score[7]+Score[8]+Score[i]);
   }
    Score[6] /= (Sum+Score[6]+Score[7]+Score[8]);
    Score[7] /= (Sum+Score[6]+Score[7]+Score[8]);
    Score[8] /= (Sum+Score[6]+Score[7]+Score[8]);
    break;
   }
  return;
}

// ----------------------
//  Default constructor.
// ----------------------
SensorMarkovProt :: SensorMarkovProt (int n, DNASeq *X) : Sensor(n)
{
  FILE *fp;
  char *tmpdir = new char[FILENAME_MAX+1];

  type = Type_Content;

  if (!IsInitialized) {
    minGC = PAR.getD("MarkovProt.minGC",GetNumber())/100;
    maxGC = PAR.getD("MarkovProt.maxGC",GetNumber())/100;
    
    maxorder = PAR.getI("MarkovProt.maxorder");
    order    = PAR.getI("MarkovProt.order");
    
    strcpy(tmpdir, PAR.getC("eugene_dir"));
    strcat(tmpdir, MODELS_DIR);
    if (! (fp = FileOpen(tmpdir, PAR.getC("MarkovProt.matname",GetNumber()), "rb"))) {
      fprintf(stderr, "cannot open matrix file %s\n", PAR.getC("matname"));
      exit(2);
    }
    
    fprintf(stderr,"Reading MarkovProt model...");
    fflush(stderr);
    ModeleProt= new TabChaine<ChainePROT21,unsigned short int> (maxorder,new ChainePROT21);
    //  ModeleProt= new TabChaine<ChainePROT21,double> (maxorder,new ChainePROT21);
    // load the coding model
    if ( ModeleProt->chargefichier(fp) ) {
      fprintf(stderr,"Proteic Model unreadable in %s. Aborting (be sure that maxorder and the matrix max. order correspond).\n",PAR.getC("matname"));
      exit(1);
    }
    fclose(fp);
    
    ProbacodonGeneral = new TabChaine<ChaineADN,double> (2,new ChaineADN);
    // initialisation is in Init, because DNASeq X is needed
    
    fprintf(stderr,"done\n");
    fflush(stderr);

    IsInitialized = true;
  }

  // for the moment, the GC rate is computed on the entire sequence
  GCrate = (X->Markov0[BitG] + X->Markov0[BitC]);
  Probacodon = ProbacodonGeneral;
  Probacodon->initialisation(GCrate);

  delete [] tmpdir;
}

// ----------------------
//  Default destructor.
// ----------------------
SensorMarkovProt :: ~SensorMarkovProt ()
{
  delete ModeleProt;
  delete Probacodon;
}

// ----------------------
//  Init MarkovProt.
// ----------------------
void SensorMarkovProt :: Init (DNASeq *X)
{
  if(PAR.getI("Output.graph")) Plot(X);
}

// --------------------------
//  GiveInfo Content MarkovProt.
// --------------------------
void SensorMarkovProt :: GiveInfo(DNASeq *X, int pos, DATA *d)
{
  int  i,j,k,n,mode=0;
  double NonCodingF,NonCodingR=0.0; // non coding forward/reverse
  char* peptid= new char[order+2];
  char* codon= new char[4];
  char tampon='\0';
  int position=pos;

  codon[3]='\0';
  peptid[order+1]='\0';

  // If the model is not in its GC% area, simply do nothing
  if ((X->Markov0[BitG] + X->Markov0[BitC]) <= minGC ||
      (X->Markov0[BitG] + X->Markov0[BitC]) > maxGC)
    return;

  // else if the current nuc. is unknow, again do nothing.
  if ((*X)[pos] == 'n') 
    return;

  // out of merges (sequence extremity), again do nothing (to be ameliorated)
  if ( (pos <= 3*(order+1)) || (pos >= (X->SeqLen - 3*(order+1))))
    return;

// ExonsF
  for(j=0;j<6;j++) {
    mode= ( (j<3)? 0 : 2 ); // mode for the DNASeq functions
    // position= relative position (strand-dependent), 
    // pos = absolute position (from the left extremity of the sequence)
    position= ( (j<3)? pos : X->SeqLen -pos -1 );
    // n is an indice corrector, frame-dependent, 
    // wich allows to obtain the begin of the "position" containing codon (for j).
    n= (j%3) - (position%3);
    if ( n>0 ) n-=3;
    k=-1;
    for (i=position +n -3*(order);i<=position;i+=3) {
      // for each codon in this frame j (i= codon first letter position) 
      k++;
      tampon= X->AA(i,mode);
      peptid[k]=tampon;
      if (tampon=='X') break; 
      //      printf("position=%d,mod=%d,i=%d,j=%d,mode=%d,n=%d,nt(i)=%c,tampon=%c\n",position,(position%3),i,j,mode,n,X->nt(i,mode),tampon);
    }
    // last codon (last AA of the peptid)
    codon[0]=X->nt(i-3,mode);
    codon[1]=X->nt(i-2,mode);
    codon[2]=X->nt(i-1,mode);
    //    printf("position=%d,frame=%d,pep= %s ,P= %f ,lP(pep)= %f ,codon=%s, nbrecodons=%d, lP/N=%f\n",position,j,peptid,ModeleProt->proba (peptid,order),ModeleProt->logproba (peptid,order),codon, NBRECODONS[Probacodon->mot2indice(codon, 3,0) -21],log( ModeleProt->proba(peptid,order)/NBRECODONS[Probacodon->mot2indice(codon, 3,0) -21] ));

    // special degenerated caracter "X";
    if (tampon == 'X')
      d->contents[j] += log(1/21) / 3;
    else {
      d->contents[j] +=
	( log(  ModeleProt->usi2real(ModeleProt->proba(peptid,order))/NBRECODONS[Probacodon->mot2indice(codon, 3,0) -21] ) ) / 3 ;
    }
  }

  // neutral model = nucleotidic probability depending on the GC rate
  NonCodingF = log (X->Markov(pos));
  NonCodingR = log (X->MarkovR(pos));

  d->contents[DATA::InterG] += (NonCodingF+NonCodingR)/2;
  d->contents[DATA::IntronF] += NonCodingF;
  d->contents[DATA::IntronR] += NonCodingR;
  d->contents[DATA::IntronUTRF] += NonCodingF;
  d->contents[DATA::IntronUTRR] += NonCodingR;
  d->contents[DATA::UTR5F] += NonCodingF;
  d->contents[DATA::UTR5R] += NonCodingR;
  d->contents[DATA::UTR3F] += NonCodingF;
  d->contents[DATA::UTR3R] += NonCodingR;
  d->contents[DATA::RNAF] += NonCodingF;
  d->contents[DATA::RNAR] += NonCodingF;

  delete peptid;
  delete codon;
  return;
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorMarkovProt :: Plot(DNASeq *TheSeq)
{
  int window, normopt;
  DATA data;
  double *Score = data.contents;
  double NScore[9], LScore[9] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,-1.0,-1.0};  
  int i,j,p;
//  fprintf(stderr,"PLOT1...");
  window = PAR.getI("Output.window")*2+1;
  normopt = PAR.getI("Output.normopt");

  for (j = 0 ; j < 9 ; j++)
    NScore[j] = Score[j] = 0.0;
  for (i = 0; i < window/2; i++) {
    GiveInfo(TheSeq,i, &data);
    for (j = 0 ; j < 9 ; j++) Score[j] = Recadre(Score[j]);
    for (j = 0 ; j < 9 ; j++) {
      NScore[j] += Score[j];
      Score[j] = 0;
    }
  }

  for (i=0; i<TheSeq->SeqLen; i++)  {
    if (i-window/2 >= 0) {
      GiveInfo(TheSeq,i-window/2, &data);
      for (j = 0 ; j < 9 ; j++)	Score[j] = Recadre(Score[j]);
      for (j = 0 ; j < 9 ; j++) {
	NScore[j] -= Score[j];
	Score[j] = 0;
      }
    }
    if (i+window/2 < TheSeq->SeqLen) {
      GiveInfo(TheSeq,i+window/2, &data);
      for (j = 0 ; j < 9 ; j++) Score[j] = Recadre(Score[j]);
      for (j = 0 ; j < 9 ; j++) 
	NScore[j] += Score[j];
    }
    for (j = 0 ; j < 9 ; j++) Score[j] = NScore[j];

    AmplifyScore(Score,normopt);
    
    if (LScore[0] < 0) for (j = 0 ; j < 9 ; j++) LScore[j] = Score[j];

    p = ((i == 0) ? 0 : i-1);
    
    PlotLine(p,i, 1, 1,LScore[0],Score[0],6);
    PlotLine(p,i, 2, 2,LScore[1],Score[1],6);
    PlotLine(p,i, 3, 3,LScore[2],Score[2],6);
    PlotLine(p,i,-1,-1,LScore[3],Score[3],6);
    PlotLine(p,i,-2,-2,LScore[4],Score[4],6);
    PlotLine(p,i,-3,-3,LScore[5],Score[5],6);
    PlotLine(p,i, 4, 4,LScore[6],Score[6],6);
    PlotLine(p,i,-4,-4,LScore[7],Score[7],6);
    PlotLine(p,i, 0, 0,LScore[8],Score[8],6);

    for (j = 0 ; j < 9 ; j++) {
      LScore[j] = Score[j];
      Score[j] = 0;
    }
  }
}
// ------------------
//  Post analyse
// ------------------
void SensorMarkovProt :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}
