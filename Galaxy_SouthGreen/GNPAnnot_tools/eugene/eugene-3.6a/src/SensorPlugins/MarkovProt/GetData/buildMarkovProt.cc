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
// $Id: buildMarkovProt.cc,v 1.8 2008-07-02 11:20:09 tschiex Exp $
// ------------------------------------------------------------------
// File:     buildMarkovProt.cc
// Contents: This program can build a proteic markov model file for eugene
// from a proteic data base file in multifasta format (typicaly SwissProt).
// A proteic markov model is generaly less accurate than a nucleic one,
// but is not species-specific, so it can be used as a general coding model,
// without requiring a species-specific genes training data set.
// A set of species specific proteins can also be used.
// For more explanation, see Foissac's PhD thesis (available in 2004 I hope...).

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "../../0_SensorTk/markov.h"
#include "../../0_SensorTk/markov.cc"
#include "../../../System.h"
#include "../../../System.cc"
#include "../../../DNASeq.h"
#include "../../../DNASeq.cc"

int ORDREMAX = 2;  // maximum of 2 recommended (to avoid overfitting)
const char* DBfile     = "/Annotation/Z_Sylvain/sprot_rel41_noFragm.fasta";  // DB infile
char* MATRIXfile = "SwissProt41.noFragm.mininfo1000.order2.usi.bin";  // outfile matrix
int OCCURENCE_MINI =1000;

int main()
{
  int i,j,k;
  FILE* file;
  TabChaine<ChainePROT21,int>* COUNT = new TabChaine<ChainePROT21,int>(ORDREMAX,new ChainePROT21);
  COUNT->initialisation();
  TabChaine<ChainePROT21,unsigned short int> MOD(ORDREMAX,new ChainePROT21);
  MOD.initialisation();

  //----------------------------------------------------------------//
  // READ the PROTEIC DBfile file (multifasta, e.g. SwissProt)
  fprintf(stderr," - opening proteic DB file...");
  file= fopen (DBfile,"rt");
  if (file==NULL) {
    printf("cannot open swissprot file %s\n",DBfile);
    exit(1);
  }
  fprintf(stderr,"done\n");

  //----------------------------------------------------------------//
  // COUNT words occurences
  fprintf(stderr," - counting words occurences...");
  COUNT->fichier2compte(file);
  fclose(file);
  fprintf(stderr,"done\n");

  //----------------------------------------------------------------//
  // ADD 1 pseudocount to each word occurence (to avoid null probs)
  //fprintf(stderr," - adding pseudocounts...");
  //  COUNT->pseudocount(1);
  // fprintf(stderr,"done\n");
  COUNT->affichage(0);  //(verbose, displays the counts)

  //----------------------------------------------------------------//
  // COMPUTE PROBS (store results in matrix MOD)
  fprintf(stderr," - computing probs...");
  MOD.compte2probas(COUNT,OCCURENCE_MINI);
  fprintf(stderr,"done\n");
  MOD.affichage(0);  //(verbose, displays the probs)

  //----------------------------------------------------------------//
  // SAVE data in the matrix file
  fprintf(stderr," - saving data in matrix file...");
  file= fopen (MATRIXfile,"wb");
  if (file==NULL) {
    printf("cannot open sp matrix file %s\n",MATRIXfile);
    exit(1);
  }
  MOD.sauve2fichier(file);
  fclose(file);
  fprintf(stderr,"done\n");

  //----------------------------------------------------------------//
  // LOAD data from a matrix file
  //  MOD.initialisation();
  //file= fopen ("../Markov/sp.modele.dou.o2.bin","rb");
  //if (file==NULL) printf("cannot open matrix file\n");
  //MOD.chargefichier(file);
  //fclose (file);
  //MOD.affichage(0);
  // g++ buildMarkovProt.cc -o buildMarkovProt
}
