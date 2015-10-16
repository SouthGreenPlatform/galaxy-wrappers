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
// $Id: CEM.cc,v 1.1 2004-10-20 12:24:21 tschiex Exp $
// ------------------------------------------------------------------
// File:     CEM.cc
// Contents: Classification of sequences according to their Markov model
// ------------------------------------------------------------------


#ifdef HAVE_CONFIG_H
#include "../../../../config.h"
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#else
#ifndef HAVE_GETOPT
#include "../../../getopt.h"
#endif
#endif

#include "../../../System.cc"
#include "../../../Const.h"
#include "../../0_SensorTk/EndianConv.h"
#include "strarray.h"
#include <vector>

// Constantes
const int  MODEL_LEN = 9;                         
const int  SIMPLE_MODEL_LEN = 6;
const int  ALPHABET_SIZE = 4;
const int  MAX_NAME_LEN = 256;
const int  MAX_STRING_LEN = 6000;
const double  MINIMUM_DELTA_VALUE = 0.1;
const int  WINDOW_SIZE = 48;

void  Find_Stop_Codons  (char [], int, int []);
int  Read_String  (FILE *, char * &, int &, char []);
void SeedInit (int seed);

int main  (int argc, char * argv []) {
  
  int ModelCount = 2; 
  int carg, errflag = 0;
  FILE  *Infile;
  char  *Line, File_Name [MAX_NAME_LEN], Name [MAX_NAME_LEN];
  int  i, j, k, Input_Size, T,NumIter = 5;
  

  SeedInit(0);

  while ((carg = getopt(argc, argv, "n:k:")) != -1) {
    switch (carg) {
    case 'k':
      ModelCount = atoi(optarg);
      break;

    case 'n':
      NumIter = atoi(optarg);
      break;
      
    default:
      errflag++;
    }
  }
  
  if (errflag) {
    fprintf (stderr, "USAGE: %s -k n Exons.fasta\n", argv [0]);
    exit (-1);
  }
  
  Line = (char *) Safe_malloc (INIT_SIZE);
  Line [0] = ' ';
  Input_Size = INIT_SIZE;
  
  // Lecture fichier de toutes les CDS
  std::vector<char *> AllCDS;
  std::vector<char *> AllNames;
  
  Infile = FileOpen (NULL,argv [optind], "r");
  while  (Read_String (Infile, Line, Input_Size, Name)) {
    int  Stop [7];
    T = strlen (Line + 1);
    assert (T < Input_Size - 1);
    // Verification de la presence de STOP en phase
    Find_Stop_Codons  (Line, T, Stop);
    
    if  (Stop [0]){
      printf ("*** CDS %s has stop codon and is ignored.\n", Name);
      break;
    }
    
    char * Tmp = new char[T+2];
    strcpy(Tmp+1,Line+1);
    Tmp[0]='\0';
    AllCDS.push_back(Tmp);
    Tmp = new char[strlen(Name)+1];
    strcpy(Tmp,Name);
    AllNames.push_back(Tmp);
  }
  printf ("Number of CDS = %ld\n", AllCDS.size());
  std::vector<int>  CDSClass(AllCDS.size(),0);
  std::vector<int>  CDSBestClass(AllCDS.size(),0);
  double OverallBestLL = -1e30;

  // the models for all classes
  String_Array  *IMM[3*ModelCount], *Lambda[3*ModelCount];

  for  (i = 0;  i < 3*ModelCount;  i ++) {
    IMM [i] = new String_Array (MODEL_LEN, ALPHABET_SIZE);
    Lambda [i] = new String_Array (MODEL_LEN - 1, ALPHABET_SIZE);
  }

  while (NumIter>0) {

    printf("Iteration %d\n",NumIter);

    // split the CDS randomly
    for (i=0; i < AllCDS.size(); i++)
      CDSClass[i] = rand() % ModelCount;

    double OverallLL = -1e30;

    while (1){

      for (k=0; k<AllCDS.size(); k++) {
	T = strlen(AllCDS[k]+1);
	Reverse_Complement(AllCDS[k],T);
      }

      // initialize  the models for each class
      for  (i = 0;  i < 3*ModelCount;  i ++) 
	IMM [i] -> Set (0.0);

      // estimate the models
      for (k=0; k<AllCDS.size(); k++) {
	
	T = strlen(AllCDS[k]+1);
	for  (j = T;  j > 0;  j --)
	  for  (i = 0;  i < MODEL_LEN && i < j;  i ++)
	    IMM [CDSClass[k]*3+(j%3)] ->Incr(AllCDS[k] + j - i, i + 1, 1.0);
      }
      
      // Normalize and condense models
      for  (i = 0;  i < 3*ModelCount;  i ++) {
	Lambda [i] -> Set_Lambda (* IMM [i]);
	IMM [i] -> Normalize ();
	IMM [i] -> Condense (* Lambda [i]);
      }
      
      // for each CDS
      int changed =0;
      double NewOverallLL = 0.0;

      for (k=0; k<AllCDS.size(); k++) {
	int Bestclass = -1;
	double BestLL = log(0.0); // -infinity
	double LL = 0.0;

	T = strlen(AllCDS[k]+1);
	Reverse_Complement(AllCDS[k],T);

	for (i=0; i< ModelCount; i++) {
	  LL=0;
	  for (j=1;j<=T;j++) 
	    LL += log(IMM[i*3+2-((j+1)%3)]->Anti(AllCDS[k]+j,Min((T-j),MODEL_LEN)));
	  //	  printf("CDS %s Model %d LL %lf\n",AllNames[k],i,LL);

	  if (LL > BestLL) {
	    BestLL=LL;
	    Bestclass = i;
	  }
	}

	//	printf("CDS %s Classe %d->%d\n",AllNames[k],CDSClass[k],Bestclass);
	if (Bestclass != CDSClass[k]) changed++;
	CDSClass[k] = Bestclass;
	NewOverallLL += BestLL;
      }

      printf("%d changes, ",changed);
      printf("Loglike = %lf\n",NewOverallLL,OverallLL);

      if (NewOverallLL > OverallBestLL) {
	OverallBestLL = NewOverallLL;
	for (i=0; i<CDSClass.size(); i++)
	  CDSBestClass[i] = CDSClass[i];
      }

      if (!changed || ((NewOverallLL-OverallLL) <= 1e-6))
	break;

	OverallLL = NewOverallLL;
    }
    NumIter--;
  }
  
  printf("Final LL = %lf\n",OverallBestLL);

  for (i=0; i<ModelCount; i++) 
    for (k=0; k < AllCDS.size(); k++) 
      if (CDSClass[k] == i) {
	printf("%s:%d %s\n",AllNames[k],CDSClass[k],AllCDS[k]+1);
      }

  return  0;
}

//  Set  Stop [0 .. 6]  TRUE  or  FALSE   according to whether
//  X [1 .. T] has a stop codon in the corresponding reading frame.
//  Stop [6]  is always set  FALSE .
void  Find_Stop_Codons  (char X [], int T, int Stop [])  {
  static int  Next [10] [4] =
    {{ 0,  2,  0,  1},     //  0  a, g
     { 3,  4,  5,  6},     //  1  t
     { 0,  2,  0,  6},     //  2  c
     { 7,  2,  7,  1},     //  3  ta
     { 8,  2,  0,  6},     //  4  tc
     { 7,  2,  0,  1},     //  5  tg
     { 9,  4,  5,  6},     //  6  ct, tt
     { 0,  2,  0,  1},     //  7  taa, tag, tga    Forward stop
     { 0,  2,  0,  1},     //  8  tca              Reverse stop
     { 7,  2,  7,  1}};    //  9  cta, tta         Reverse stop
  int  i, State;
  
  for  (i = 0;  i < 7;  i ++)
    Stop [i] = 0;
  
  State = 0;
  for  (i = 1;  i <= T;  i ++) {
    switch  (tolower (X [i])) {
    case  'a' :
      State = Next [State] [0];
      break;
    case  'c' :
      State = Next [State] [1];
      break;
    case  'g' :
      State = Next [State] [2];
      break;
    case  't' :
      State = Next [State] [3];
      break;
    default :
      fprintf (stderr, "Unexpected character %c\n", X [i]);
      State = 0;
    }
    if  (State == 7)
      Stop [i % 3] = TRUE;
    else if  (State > 7)
      Stop [3 + i % 3] = TRUE;
  }
  
  return;
}

/* Read next string from fp (assuming one string per line with each
*  string preceded by a name tag) into T [0 ..]  which has Size
*  characters.  Allocate extra memory if needed and adjust Size
*  accordingly.  Returns TRUE if successful, FALSE otherwise (e.g.,
*  EOF). */
int  Read_String  (FILE * fp, char * & T, int & Size, char Name []) {
  
  long int  Len;
  int  Ch;
  
  if  (fscanf (fp, "%s", Name) == EOF)
    return  FALSE;
  
  while  ((Ch = fgetc (fp)) != EOF && (Ch == ' ' || Ch == '\t'))
    ;
  if  (Ch == EOF)
    return  FALSE;

  ungetc (Ch, fp);
  T [0] = '\0';
  Len = 1;
  while  ((Ch = fgetc (fp)) != EOF && Ch != '\n') {
    if  (isspace (Ch))
      continue;
    
    if  (Len >= Size) {
      Size += INCR_SIZE;
      T = (char *) Safe_realloc (T, Size);
    }
    Ch = tolower (Ch);
    switch  (Ch)
      {
      case  'a' :
      case  'c' :
      case  'g' :
      case  't' :
	break;
      case  's' :
	Ch = 'c';
	break;
      case  'w' :
	Ch = 'a';
	break;
      case  'r' :
	Ch = 'a';
	break;
      case  'y' :
	Ch = 'c';
	break;
      case  'm' :
	Ch = 'a';
	break;
      case  'k' :
	Ch = 'g';
	break;
      case  'b' :
	Ch = 'c';
	break;
      case  'd' :
	Ch = 'a';
	break;
      case  'h' :
	Ch = 'a';
	break;
      case  'v' :
	Ch = 'a';
	break;
      case  'n' :
	Ch = 'a';
	break;
      default :
	fprintf (stderr, "Unexpected character `%c\' in string %s\n",
		 Ch, Name);
	Ch = 'a';
      }
    T [Len ++] = Ch;
  }
  
  T [Len] = '\0';
  
  return  TRUE;
}
/* --------------------------------------------------------------------
// Random generator: seed initialisation.  If seed = 0, generates a
// random seed from pid and time else use the given seed
// -------------------------------------------------------------------- */
void SeedInit (int seed)
{
  unsigned int heure;
  
  if (!seed)
  {
    heure = (unsigned int) time(NULL);
    srand(heure);
  }
  else
    srand(seed);
}

