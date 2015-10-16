
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
// $Id: Sensor.MarkovIMM.cc,v 1.17 2010-01-25 17:01:23 sallet Exp $
// ------------------------------------------------------------------
// File:     Sensor.MarkovIMM.cc
// Contents: Sensor MarkovIMM
// ------------------------------------------------------------------

#include "Sensor.MarkovIMM.h"

extern Parameters PAR;

/*************************************************************
 **                     SensorMarkovIMM                     **
 *************************************************************/

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
    // pas de normalisation
    
    break;
    
  case 1:
    // on suppose que une phase au plus peut coder
    for  (i = 0;  i < 9;  i ++) {
      Score[i] = exp(Score[i]-Min);
      Sum += Score[i];
    }
    
    for  (i = 0;  i < 9;  i ++) 
      Score [i] /= Sum;
    break;
    
  case 2:
    // chaque phase peut coder independamment
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
SensorMarkovIMM :: SensorMarkovIMM (int n, DNASeq *X) : Sensor(n)
{
  FILE *fp;
  int i;
  std::vector<BString_Array*> IMMatrix;
  std::string matrixName;
  bool is_initialised = false;
  char *tmpdir = new char[FILENAME_MAX+1];

  type = Type_Content;


  IntergenicModel = (PAR.getI("MarkovIMM.IntergenicModel",GetNumber()));
  
  maxOrder = PAR.getI("MarkovIMM.maxOrder",GetNumber());  
  minGC = PAR.getD("MarkovIMM.minGC",GetNumber())/100;
  maxGC = PAR.getD("MarkovIMM.maxGC",GetNumber())/100;

  matrixName = PAR.getC("MarkovIMM.matname",GetNumber());

  // Look if the same matrice has already been loaded
  // The name of the matrix is considered as enough to define it.
  for (unsigned int k=0; k<IMMatrixList.size(); k++) 
    if (matrixName == matrixNameList[k]) {
      refCount[k]++;
      IMMatrix_index = k;
      is_initialised = true;
      k = IMMatrixList.size();
    }

  if (!is_initialised) { 
    IMMatrix_index = IMMatrixList.size();

    // update a local variable, copied in the list when initialized
    for (i=0; i<7; i++) IMMatrix.push_back(NULL); 

    strcpy(tmpdir, PAR.getC("eugene_dir"));
    strcat(tmpdir, MODELS_DIR);
    if (! (fp = FileOpen(tmpdir , PAR.getC("MarkovIMM.matname",GetNumber()), "rb"))) {
      fprintf(stderr, "cannot open matrix file %s\n", PAR.getC("MarkovIMM.matname"));
      exit(2);
    }
    
    fprintf(stderr,"Reading IMM...");
    fflush(stderr);
    
    // On essaie d'abord de charger les 5 modeles fondamentaux (3 ex/int/interG)
    for  (i = 0;  i < 5;  i ++) {
      IMMatrix[i] = new BString_Array(MODEL_LEN, ALPHABET_SIZE);
      if (IMMatrix[i]->Read(fp)) {
	fprintf(stderr,"Model %d unreadable in %s. Aborting.\n",i+1,PAR.getC("MarkovIMM.matname"));
	exit(1);
      } 
      fprintf(stderr," %d",i+1);
      fflush(stderr);
    }
    
    // On essaie ensuite de lire un 6eme modele. Si cela echoue,
    // le modele intronique est utilise pour les UTR
    IMMatrix[6] = new BString_Array(MODEL_LEN, ALPHABET_SIZE);
    if (IMMatrix[6]->Read(fp)) {
      fprintf(stderr,"- No UTR model found, using intronic model. ");
      delete IMMatrix[6];
      IMMatrix[6] = IMMatrix[3];
      IMMatrix[5] = IMMatrix[3];
    } else {
      fprintf(stderr," 6");
      IMMatrix[5] = new BString_Array(MODEL_LEN, ALPHABET_SIZE);
      if (IMMatrix[5]->Read(fp)) {
	fprintf(stderr,"- No second UTR model found, using intronic model. ");
	delete IMMatrix[5];
	IMMatrix[5] = IMMatrix[3];
      } else fprintf(stderr," 7");
    }  
    fprintf(stderr," ...done\n");
    fclose(fp);

    IMMatrixList.push_back( IMMatrix );
    matrixNameList.push_back(matrixName);
    refCount.push_back(1);
  }

  delete [] tmpdir;
}

// ----------------------
//  Default destructor.
// ----------------------
SensorMarkovIMM :: ~SensorMarkovIMM ()
{

  // If nobody uses the Matrix model, delete
  if (--refCount[IMMatrix_index] == 0){
    
    if (IMMatrixList[IMMatrix_index][5] != IMMatrixList[IMMatrix_index][3]) 
      delete  IMMatrixList[IMMatrix_index][5];

    if (IMMatrixList[IMMatrix_index][6] != IMMatrixList[IMMatrix_index][3]) 
      delete  IMMatrixList[IMMatrix_index][6];
  
    for  (int i = 0;  i < 5;  i ++)
      delete  IMMatrixList[IMMatrix_index][i];
  }
}

// ----------------------
//  Init MarkovIMM.
// ----------------------
void SensorMarkovIMM :: Init (DNASeq *X)
{
  if(PAR.getI("Output.graph")) Plot(X);
}

// ----------------------------
//  GiveInfo Content MarkovIMM.
// ----------------------------
void SensorMarkovIMM :: GiveInfo(DNASeq *X, int pos, DATA *d)
{
  int  Rev,rModelLen,fModelLen;
  int  indexF, indexR;
  std::vector<BString_Array*> IMMatrix; // to locally store the used IMMatrix

  IMMatrix = IMMatrixList[IMMatrix_index];

  // If the model is not in its GC% area, simply do nothing
  if ((X->Markov0[BitG] + X->Markov0[BitC]) <= minGC ||
      (X->Markov0[BitG] + X->Markov0[BitC]) > maxGC)
    return;

  // else if the current nuc. is unknown, again do nothing.
  if ((*X)[pos] == 'n') 
    return;

  // Compute possible maximum order
  Rev = X->SeqLen - pos;
  fModelLen = Min(maxOrder+1,Min(Rev,MODEL_LEN));
  rModelLen = Min(maxOrder+1,Min(pos+1,MODEL_LEN));
  
  // and indexes in the IMMatrix models
  indexF = (*IMMatrix[0]).AntiString_To_Sub(X,pos,fModelLen);
  indexR = (*IMMatrix[0]).String_To_Sub(X,pos-rModelLen+1,rModelLen);

  // Exons F
  d->contents[DATA::ExonF1] += log((double)(*IMMatrix[2-((pos+2)%3)])[indexF]/65535.0);
  d->contents[DATA::ExonF2] += log((double)(*IMMatrix[2-((pos+1)%3)])[indexF]/65535.0);
  d->contents[DATA::ExonF3] += log((double)(*IMMatrix[2-(pos%3)])[indexF]/65535.0);
  
  // Exons R
  d->contents[DATA::ExonR1] += log((double)(*IMMatrix[2-((Rev+1)%3)])[indexR]/65535.0);
  d->contents[DATA::ExonR2] += log((double)(*IMMatrix[2-(Rev%3)])[indexR]/65535.0);
  d->contents[DATA::ExonR3] += log((double)(*IMMatrix[2-((Rev+2)%3)])[indexR]/65535.0);
  
  // Introns F/R
  d->contents[DATA::IntronF] += log((double)(*IMMatrix[3])[indexF]/65535.0);
  d->contents[DATA::IntronR] += log((double)(*IMMatrix[3])[indexR]/65535.0);
  d->contents[DATA::IntronUTRF] += log((double)(*IMMatrix[3])[indexF]/65535.0);
  d->contents[DATA::IntronUTRR] += log((double)(*IMMatrix[3])[indexR]/65535.0);
  
  // InterG and RNA
  switch (IntergenicModel)
    {
      // Before 3.4, known as useM0asIG=true: uses a Oth order markov
      // model with seq. GC%
    case 0:  
      d->contents[DATA::InterG] += log(X->GC_AT(pos));
      d->contents[DATA::RNAF]   += log(X->GC_AT(pos));
      d->contents[DATA::RNAR]   += log(X->GC_AT(pos));
      break;
     
      // new assymetric model (3.4). This breaks the syymetry between
      // reverse and forward but seems top work better on C. Elegans
    case 1:
      d->contents[DATA::InterG] += log(((double)(*IMMatrix[4])[indexF]/65535.0));
      d->contents[DATA::RNAF]   += log(((double)(*IMMatrix[4])[indexF]/65535.0));
      d->contents[DATA::RNAR]   += log(((double)(*IMMatrix[4])[indexF]/65535.0));
      break;

      // old symetric model (default before 3.4). Mixes forward and reverse probabilities.
    case 2:
      d->contents[DATA::InterG] += log(((double)(*IMMatrix[4])[indexF] + (double)(*IMMatrix[4])[indexR])/131071.0);
      d->contents[DATA::RNAF]   += log(((double)(*IMMatrix[4])[indexF] + (double)(*IMMatrix[4])[indexR])/131071.0);
      d->contents[DATA::RNAR]   += log(((double)(*IMMatrix[4])[indexF] + (double)(*IMMatrix[4])[indexR])/131071.0);
      break;

    default:
      fprintf(stderr,"Error: Incorrect value for parameter IntergenicModel\n");
      exit(1);
    }
  
  // UTR 5' F/R
  d->contents[DATA::UTR5F] += log((double)(*IMMatrix[6])[indexF]/65535.0);
  d->contents[DATA::UTR5R] += log((double)(*IMMatrix[6])[indexR]/65535.0);
  
  // UTR 3' F/R
  d->contents[DATA::UTR3F] += log((double)(*IMMatrix[5])[indexF]/65535.0);
  d->contents[DATA::UTR3R] += log((double)(*IMMatrix[5])[indexR]/65535.0);
  
  return;
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorMarkovIMM :: Plot(DNASeq *TheSeq)
{
  int window, normopt;
  DATA data;
  double *Score = data.contents;
  double NScore[9], LScore[9] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,-1.0,-1.0};  
  int i,j,p;

  window = PAR.getI("Output.window")*2+1;
  normopt = PAR.getI("Output.normopt");
  
  for (j = 0 ; j < 9 ; j++)
    NScore[j] = Score[j] = 0.0;
  
  for (i = 0; i < window/2; i++) {
    GiveInfo(TheSeq,i, &data);
    for (j = 0 ; j < 9 ; j++) {
      NScore[j] += Score[j];
      Score[j] = 0;
    }
  }

  for (i=0; i<TheSeq->SeqLen; i++)  {

    if (i-window/2 >= 0) {
      GiveInfo(TheSeq,i-window/2, &data);
      for (j = 0 ; j < 9 ; j++) {
	NScore[j] -= Score[j];
	Score[j] = 0;
      }
    }

    if (i+window/2 < TheSeq->SeqLen) {
      GiveInfo(TheSeq,i+window/2, &data);
      for (j = 0 ; j < 9 ; j++) 
	NScore[j] += Score[j];
    }
    
    for (j = 0 ; j < 9 ; j++) Score[j] = NScore[j];
    AmplifyScore(Score,normopt);
    
    if (LScore[0] < 0) for (j = 0 ; j < 9 ; j++) LScore[j] = Score[j];

    p = ((i == 0) ? 0 : i-1);
    
    PlotLine(p,i, 1, 1,LScore[0],Score[0],3);
    PlotLine(p,i, 2, 2,LScore[1],Score[1],3);
    PlotLine(p,i, 3, 3,LScore[2],Score[2],3);
    PlotLine(p,i,-1,-1,LScore[3],Score[3],3);
    PlotLine(p,i,-2,-2,LScore[4],Score[4],3);
    PlotLine(p,i,-3,-3,LScore[5],Score[5],3);
    PlotLine(p,i, 4, 4,LScore[6],Score[6],3);
    PlotLine(p,i,-4,-4,LScore[7],Score[7],3);
    PlotLine(p,i, 0, 0,LScore[8],Score[8],3);
      
    for (j = 0 ; j < 9 ; j++) {
      LScore[j] = Score[j];
      Score[j] = 0;
    }
  }
}
// ------------------
//  Post analyse
// ------------------
void SensorMarkovIMM :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}
