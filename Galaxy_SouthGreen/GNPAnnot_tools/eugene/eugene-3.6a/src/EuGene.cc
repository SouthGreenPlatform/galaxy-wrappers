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
// $Id: EuGene.cc,v 1.105 2009-04-22 13:54:57 tschiex Exp $
// ------------------------------------------------------------------
// File:     EuGene.cc
// Contents: This program finds exons/introns and intergenic regions 
//           (including UTRs)
// ------------------------------------------------------------------


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstdio>
#include <cstdlib>

#ifdef __APPLE__
// MacOS-X kludge. cmath undefines these macros. Turn them into inlines 
#include <math.h>
inline int (isinf)(double r) { return isinf(r); }
inline int (isnan)(double r) { return isnan(r); }
#endif

#include <cmath>
#include <cctype>
#include <climits>
#include <cfloat>
#include <ctime>
#include <cassert>
#include <cerrno>
#include <unistd.h>
#include <string.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif

#include "Const.h"
#include "System.h"
#include "DAG.h"
#include "AltEst.h"
#include "SensorIF.h"
#include "MSensor.h"
#include "Param.h"
#include "DNASeq.h"
#include "BackP.h"
#include "PenaltyDist.h"
#include "Prediction.h"
#include "Parametrization/ParaOptimization.h"

#ifndef FILENAME_MAX
#define FILENAME_MAX        1024
#endif

// ------------------ Globals --------------------
MasterSensor*    MS;
Parameters       PAR;
ParaOptimization OPTIM;

// -------------------------------------------------------------------------
// Compute Optimal Prediction
// -------------------------------------------------------------------------
Prediction* Predict (DNASeq* TheSeq, int From, int To, MasterSensor* MSensor)
{

  int   j, k;
  int   Data_Len = TheSeq->SeqLen;
  DATA	Data;
  int   Forward =  1;//PAR.getI("Sense"); 
  int   Dir = (Forward ? 1 : -1);

  int   GCVerbose = PAR.getI("EuGene.VerboseGC");
  int   GCLatency = PAR.getI("EuGene.GCLatency");
  DAG   *Dag;
  Prediction *prediction;
	
  // DynaProg end at the lastpos + 1 to account for final signals.
  int	FirstNuc = (Forward ? From : To+1);
  int   LastNuc  = (Forward ? To+1 : From);

  // Initialize the central DAG
  // one more position here for prediction of partial elements at extremities.
  Dag = new DAG(FirstNuc-Dir, LastNuc+Dir, PAR,TheSeq);
  
  Dag->LoadDistLength();
  Dag->WeightThePrior();

  // --------------------------------------------------------------------------
  // Demarrage de la programmation dynamique
  // --------------------------------------------------------------------------
  for (int nuc = FirstNuc; nuc != LastNuc+Dir; nuc += Dir) {

    // recuperation des infos
    MSensor->GetInfoAt(TheSeq, nuc, &Data);
    if (Forward) 
      Dag->ShortestPathAlgoForward(nuc,Data);
    else
      Dag->ShortestPathAlgoBackward(nuc,Data);

    if (nuc && (nuc % GCLatency == 0)) Dag->MarkAndSweep(nuc,GCVerbose,GCLatency);
  }

  Dag->WeightThePrior();
  Dag->BuildPrediction(From,To,Forward);

  Dag->pred->TrimAndUpdate(TheSeq);
  prediction = Dag->pred;

  delete Dag;

  return prediction;
}

// Same with whole sequence by default
Prediction* Predict (DNASeq* TheSeq, MasterSensor* MSensor)
{
	return Predict(TheSeq,0,TheSeq->SeqLen-1,MSensor);
}
// -------------------------------------------------------------------------
// Compute alternative Predictions based on EST
// -------------------------------------------------------------------------
Prediction* AltPredict (DNASeq* TheSeq, int From, int To, MasterSensor* MSensor, 
		AltEst *AltEstDB, Prediction *optpred, int idx)
{
  int   j, k;
  int   Data_Len = TheSeq->SeqLen;
  DATA	Data;
  int   Forward =  1;//PAR.getI("Sense"); 
  int   Dir = (Forward ? 1 : -1);
  int   GCVerbose = PAR.getI("EuGene.VerboseGC");
  int   GCLatency = PAR.getI("EuGene.GCLatency");
  DAG* Dag;
  Prediction *prediction;

   // DynaProg end at the lastpos + 1 to account for final signals.
  int	FirstNuc = (Forward ? From : To+1);
  int   LastNuc  = (Forward ? To+1 : From);

  if (!AltEstDB->voae_AltEst[idx].CompatibleWith(optpred))
    {

      Dag = new DAG(FirstNuc-Dir, LastNuc+Dir, PAR,TheSeq);
	  
      Dag->LoadDistLength();
      Dag->WeightThePrior();

      for (int nuc = FirstNuc; nuc != LastNuc+Dir; nuc += Dir) {
	
	// recuperation des infos
	MSensor->GetInfoAt(TheSeq, nuc, &Data);
	AltEstDB->Penalize(idx,nuc,&Data);

	if (Forward) 
	  Dag->ShortestPathAlgoForward(nuc,Data);
	else
	  Dag->ShortestPathAlgoBackward(nuc,Data);
	
	if (nuc && (nuc % GCLatency == 0)) Dag->MarkAndSweep(nuc,GCVerbose,GCLatency);
      }
      Dag->WeightThePrior();
      Dag->BuildPrediction(From, To, Forward);

      Dag->pred->TrimAndUpdate(TheSeq);
      prediction = Dag->pred;

      delete Dag;
      return prediction;
    }
  return NULL;

}

// -------------------------------------------------------------------------
// Read a fasta file
// -------------------------------------------------------------------------
DNASeq* ReadSequence (char* sequence_name)
{
  DNASeq *TheSeq;
  
  fprintf(stderr,"-------------------------------------");
  fprintf(stderr,"--------------------------------\nLoading sequence...");
  fflush(stderr);
  
  TheSeq = new DNASeq(sequence_name);
  
  fprintf(stderr,"%s, %d bases read, ",TheSeq->Name, TheSeq->SeqLen);
  
  fprintf (stderr,"GC Proportion = %.1f%%\n", 
	   (TheSeq->Markov0[BitG]+TheSeq->Markov0[BitC])*100.0);
  
  return TheSeq;
}

// -------------------------------------------------------------------------
//            MAIN
// -------------------------------------------------------------------------
int main  (int argc, char * argv [])
{
  DNASeq     *TheSeq;
  int        Data_Len;
  Prediction *pred;
  FILE       *MISC_INFO;
  char       prefixName[FILENAME_MAX+1];
  char       grname[FILENAME_MAX+1];
  char       miname[FILENAME_MAX+1];
  int        graph;

  fprintf(stderr,"-------------------------------------"
	  "--------------------------------\n");

  // Lecture de la ligne d'arg et du fichier .par
  PAR.initParam(argc, argv);

  if (PAR.getI("ParaOptimization.Use")) 
    OPTIM.ParaOptimize(argc, argv);
  else 
    {
      // Objectif : limiter les appels à la MAP
      graph = PAR.getI("Output.graph");

      int sequence;
      for (sequence = optind; sequence < argc ; sequence++) {
	
	PAR.set("fstname", argv[sequence]);

	// --------------------------------------------------------------------
	// Lecture de la sequence    
	// --------------------------------------------------------------------
	TheSeq = ReadSequence( PAR.getC("fstname") );
	Data_Len = TheSeq->SeqLen;
    
	// --------------------------------------------------------------------
	// Calcul des positions de prédiction
	// --------------------------------------------------------------------
 	int fromPos = PAR.getI("EuGene.from",0,true);
  	int toPos   = PAR.getI("EuGene.to",0,true);

  	if (toPos) 
		toPos = Min(Data_Len-1,toPos);
	else 
		toPos = Data_Len-1;

	// --------------------------------------------------------------------
	// Prefix output file name
	// --------------------------------------------------------------------
	strcpy(prefixName, PAR.getC("Output.Prefix"));
	strcat(prefixName, BaseName(PAR.getC("fstname")));
	if ( rindex(prefixName, '.') != NULL ) {
	  if (!strcmp(rindex(prefixName, '.'), ".fasta") ||
	      !strcmp(rindex(prefixName, '.'), ".fsa")   ||
	      !strcmp(rindex(prefixName, '.'), ".tfa")   ||
	      !strcmp(rindex(prefixName, '.'), ".txt"))
	    *rindex(prefixName, '.') = 0;     // on enleve l'extension
	}
	PAR.set("prefixName", prefixName);

	// --------------------------------------------------------------------
	// Preparation sortie graphique + Scores
	// --------------------------------------------------------------------
	if (graph) {
	  int gto       = PAR.getI("Output.gto");
	  int gfrom     = PAR.getI("Output.gfrom");
	  int glen      = PAR.getI("Output.glen");
      		  
	  // Construction du nom de sortie (*.png)
	  strcpy(grname, prefixName);
	  
	  if ((gfrom <= 0)|| (gfrom >= Data_Len))
	    gfrom = 1;
	  if ((gto <= 0)  || (gto <= gfrom) || (gto > Data_Len))
	    gto = Data_Len;
	  
	  if ((PAR.getI("Output.gfrom")!=-1) || PAR.getI("Output.gto")!=-1) {
	    sprintf(grname+strlen(grname), ".%d", gfrom);
	    sprintf(grname+strlen(grname), "-%d", gto);
	  }
	  
	  gfrom--;
	  gto--;
	  
	  if (glen < 0)
	    glen = ((gto-gfrom+1 <= 6000) ? gto-gfrom+1 : 6000);
	  
	  InitPNG(PAR.getI("Output.resx"),   PAR.getI("Output.resy"),
		  PAR.getI("Output.offset"), gfrom, gto,
		  PAR.getI("Output.golap"), glen, grname);
	}
	
	// --------------------------------------------------------------------
	// Init MasterSensor
	// --------------------------------------------------------------------
	MS = new MasterSensor();
	MS->InitMaster(TheSeq);

	// --------------------------------------------------------------------
	// Predict: 1st main prediction
	// --------------------------------------------------------------------
	pred = Predict(TheSeq, fromPos, toPos, MS);
	fprintf(stderr,"Optimal path length = %.4f\n",- pred->optimalPath);
	
	// --------------------------------------------------------------------
	// Textual and graphical output
	// --------------------------------------------------------------------
	if (graph)
	  pred->PlotPred();

	if ( ! PAR.getI("AltEst.use") )
		pred->Print(TheSeq, MS);

	strcpy(miname, prefixName);
	MISC_INFO = FileOpen(NULL, strcat(miname, ".misc_info"), "wb");
	pred->PrintGeneInfo(MISC_INFO);
	MS->PostAnalyse(pred, MISC_INFO);
	fclose(MISC_INFO);

    	if (graph) {
	  fprintf(stderr,"Dumping images (\"%s.---.png\")...", grname);
	  fflush(stderr);
	  ClosePNG();
	  fprintf(stderr, "done\n");
	}

	// --------------------------------------------------------------------
	// Load Alternative EST data (if any)
	// --------------------------------------------------------------------
	if (PAR.getI("AltEst.use")) {

	  int ExonBorderMatchThreshold = PAR.getI("AltEst.ExonBorderMatchThreshold");
	  int RepredictMargin = PAR.getI("AltEst.RepredictMargin");
	  int newGene = 0; // if a splice variant has no base gene, it is a "new" gene. counter needed for gene number
	  AltEst *AltEstDB = new AltEst(TheSeq);

	  std::vector <Prediction*> vPred;
	  Prediction* AltPred;
	  Gene *baseGene;
	  for (int altidx = 0; altidx < AltEstDB->totalAltEstNumber; altidx++)
	    {
	      int localFrom,localTo;

	      localFrom = Max(fromPos,AltEstDB->voae_AltEst[altidx].GetStart()-RepredictMargin);
	      localTo = Min(toPos, AltEstDB->voae_AltEst[altidx].GetEnd()+RepredictMargin);

	      AltPred = AltPredict(TheSeq,localFrom,localTo,MS,AltEstDB,pred,altidx);

	      if (AltPred) {
            if ( (AltPred->vGene[0]->cdsStart == -1) || (AltPred->vGene[0]->cdsEnd == -1))
            {
              delete AltPred;
              continue;
          	}
			//fprintf(stderr,"start=%d, end=%d\n",AltEstDB->voae_AltEst[altidx].GetStart(),AltEstDB->voae_AltEst[altidx].GetEnd());
		AltPred->DeleteOutOfRange(AltEstDB->voae_AltEst[altidx].GetStart(),AltEstDB->voae_AltEst[altidx].GetEnd());
		if (AltPred->IsOriginal(pred,vPred,ExonBorderMatchThreshold))
		  {
		    fprintf(stderr,"Optimal path length = %.4f\n",- AltPred->optimalPath);
		    baseGene = pred->FindGene(AltEstDB->voae_AltEst[altidx].GetStart(),AltEstDB->voae_AltEst[altidx].GetEnd());
		    if (baseGene) {
		    	baseGene->hasvariant++;
		    	AltPred->vGene[0]->isvariant = true;
		    	AltPred->vGene[0]->hasvariant = baseGene->hasvariant;
		    	AltPred->vGene[0]->geneNumber = baseGene->geneNumber;
			//fprintf(stderr,"Alt: %d,hasv=%d,isv=%d Opt gn=%d,hasv=%d,isv=%d\n",AltPred->vGene[0]->geneNumber,AltPred->vGene[0]->hasvariant,AltPred->vGene[0]->isvariant,baseGene->geneNumber,baseGene->hasvariant,AltPred->vGene[0]->isvariant);
			baseGene->tuStart = ( baseGene->tuStart ) ? Min(baseGene->tuStart,AltPred->vGene[0]->trStart) 
													  : Min(baseGene->trStart,AltPred->vGene[0]->trStart);
			baseGene->tuEnd   = ( baseGene->tuEnd )   ? Max(baseGene->tuEnd,AltPred->vGene[0]->trEnd)     
								  : Max(baseGene->trEnd,AltPred->vGene[0]->trEnd);
		    //AltPred ->Print(TheSeq, MS,NULL,1);
			}
			else 
			{
				fprintf(stderr,"New gene predicted by alternative spliced gene prediction.\n");
				AltPred->vGene[0]->geneNumber = pred->nbGene + newGene++;
			}
		    vPred.push_back(AltPred);
		  }
		else delete AltPred;	
	      }
	    }
		pred->Print(TheSeq, MS);
		for (int idx = 0; idx < vPred.size(); idx++) 
		{
			vPred[idx]->Print(TheSeq, MS,NULL,1);
		}
	  //delete AltPred;	
	}

	// Free used memory
	delete TheSeq; 
	delete MS;
	delete pred;
	
	fflush(stderr);
	fflush(stdout);
      } // fin de traitement de chaque séquence....

      fprintf(stderr,"-------------------------------------"	      "--------------------------------\n");

      return  0;
    }
}
