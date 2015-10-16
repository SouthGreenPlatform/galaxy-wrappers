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
// $Id: TrainIMM.cc,v 1.5 2006-11-17 08:14:02 gouzy Exp $
// ------------------------------------------------------------------
// File:     TrainIMM.cc
// Contents: Estimation of parameters (inspired of Glimmer (Salzberg))
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


// Constantes
const int  MODEL_LEN = 9;
const int  SIMPLE_MODEL_LEN = 6;
const int  ALPHABET_SIZE = 4;
const int  MAX_NAME_LEN = 256;
const int  MAX_STRING_LEN = 6000;
const double  MINIMUM_DELTA_VALUE = 0.1;
const int  WINDOW_SIZE = 48;



void  Find_Stop_Codons  (char [], int, int []);
int  Read_String  (FILE *, char * &, long int &, char []);

int main  (int argc, char * argv [])
{
  int ModelCount=5; // 3 exons, introns, intergenique
  int carg, errflag = 0;
  FILE  *Infile;
  String_Array  *IMM[7], *Lambda[7];
  char  *Line, File_Name [MAX_NAME_LEN], Name [MAX_NAME_LEN];
  char *UTR5File,*UTR3File,*IgFile;
  long int  i, j, Input_Size, Num_Strings, T;
  double IgWeight = 1.0;

  UTR3File = UTR5File = IgFile = NULL;

  while ((carg = getopt(argc, argv, "h5:3:I:")) != -1)
  {
    switch (carg)
    {
    case 'h':
      IgWeight = 0.5;
      break;

    case '5':   // a 5'UTR file is available
      UTR5File = optarg;
      ModelCount++;
      break;

    case '3':
      UTR3File = optarg;
      ModelCount++;
      break;

    case 'I':
      IgFile = optarg;
      break;
    default:
      errflag++;
    }
  }

  if (((argc - optind) != 2) || (UTR5File && UTR3File == NULL))
    errflag++;

  if (errflag)
  {
    fprintf (stderr, "USAGE: %s Exons.fasta Introns.fasta [-3 UTR3.fasta [-5 UTR5.fasta]]\n", argv [0]);
    fprintf(stderr,"            [-I Intergenic.fasta].\n");
    fprintf(stderr,"NB: if an UTR5 file is provided, the UTR3 must also be provided.\n");
    exit (-1);
  }

  for  (i = 0;  i < 7;  i ++)
  {
    IMM [i] = new String_Array (MODEL_LEN, ALPHABET_SIZE);
    IMM [i] -> Set (0.0);
  }
  Num_Strings = 0;
  Line = (char *) Safe_malloc (INIT_SIZE);
  Line [0] = ' ';
  Input_Size = INIT_SIZE;

  // Lecture fichier exon
  Infile = FileOpen (NULL,argv [optind], "r");
  while  (Read_String (Infile, Line, Input_Size, Name))
  {
    int  Stop [7];
    T = strlen (Line + 1);
    assert (T < Input_Size - 1);
    // Verification de la presence de STOP en phase
    Find_Stop_Codons  (Line, T, Stop);
    if  (Stop [0])
      printf ("*** String %ld has stop codon\n", Num_Strings + 1);

    // On utilise des modeles avec contexte en avant.
    Reverse_Complement (Line, T);
    for  (j = T;  j > 0;  j --)
      for  (i = 0;  i < MODEL_LEN && i < j;  i ++)
        IMM [j % 3] ->Incr(Line + j - i, i + 1, 1.0);
    Num_Strings ++;
  }
  printf ("Number of exon strings = %ld\n", Num_Strings);
  Num_Strings = 0;

  // Les introns et l'intergenique s'il n'y en a pas
  if (!IgFile) printf ("Intergenic model built by intron autocomplementation.\n");
  Infile = FileOpen (NULL,argv[optind+1], "r");
  while  (Read_String (Infile, Line, Input_Size, Name))
  {
    T = strlen (Line + 1);
    assert (T < Input_Size - 1);

    // sans intergenique on mixe les introns dans les deux sens
    if (!IgFile)
      for  (j = T;  j > 0;  j --)
        for  (i = 0;  i < MODEL_LEN && i < j;  i ++)
          IMM [4] ->Incr(Line + j - i, i + 1, IgWeight);
    Reverse_Complement (Line, T);
    // intron et mixage intergenic
    for  (j = T;  j > 0;  j --)
      for  (i = 0;  i < MODEL_LEN && i < j;  i ++)
      {
        IMM [3] ->Incr(Line + j - i, i + 1, 1.0);
        if (!IgFile) IMM [4]->Incr(Line + j - i, i + 1, IgWeight);
      }
    Num_Strings ++;
  }
  printf ("Number of introns strings = %ld\n", Num_Strings);
  Num_Strings = 0;

  // Intergenic (if provided)
  if (IgFile)
  {
    Infile = FileOpen(NULL,IgFile, "r");
    while  (Read_String (Infile, Line, Input_Size, Name))
    {
      T = strlen (Line + 1);
      assert (T < Input_Size - 1);
      Reverse_Complement (Line, T);
      for  (j = T;  j > 0;  j --)
        for  (i = 0;  i < MODEL_LEN && i < j;  i ++)
          IMM [4] ->Incr(Line + j - i, i + 1, 1.0);
      Num_Strings ++;
    }
    printf ("Number of Intergenic strings = %ld\n", Num_Strings);
    Num_Strings = 0;
  }


  // UTR 5
  if (UTR5File)
  {
    Infile = FileOpen(NULL,UTR5File, "r");
    while  (Read_String (Infile, Line, Input_Size, Name))
    {
      T = strlen (Line + 1);
      assert (T < Input_Size - 1);
      Reverse_Complement (Line, T);
      for  (j = T;  j > 0;  j --)
        for  (i = 0;  i < MODEL_LEN && i < j;  i ++)
          IMM [5] ->Incr(Line + j - i, i + 1, 1.0);
      Num_Strings ++;
    }
    printf ("Number of 5UTR strings = %ld\n", Num_Strings);
    Num_Strings = 0;
  }

  // UTR3
  if (UTR3File)
  {
    Infile = FileOpen(NULL,UTR3File, "r");
    while  (Read_String (Infile, Line, Input_Size, Name))
    {
      T = strlen (Line + 1);
      assert (T < Input_Size - 1);
      Reverse_Complement (Line, T);
      for  (j = T;  j > 0;  j --)
        for  (i = 0;  i < MODEL_LEN && i < j;  i ++)
          IMM [6] ->Incr(Line + j - i, i + 1, 1.0);
      Num_Strings ++;
    }
    printf ("Number of 3UTR strings = %ld\n", Num_Strings);
    Num_Strings = 0;
  }

  for  (i = 0;  i < ModelCount;  i ++)
  {
    Lambda [i] = new String_Array (MODEL_LEN - 1, ALPHABET_SIZE);
    Lambda [i] -> Set_Lambda (* IMM [i]);
    IMM [i] -> Normalize ();
    IMM [i] -> Condense (* Lambda [i]);
  }

  strcat(strcpy (File_Name, argv[optind]), ".bin");
  FILE *fp;
  if (! (fp = FileOpen(NULL,File_Name, "wb")))
  {
    fprintf(stderr, "cannot open output file %s\n",  File_Name);
    exit(2);
  }

  printf("Dumping models in %s\n",File_Name);
  for  (i = 0;  i < ModelCount;  i ++)
    IMM [i] -> BinWrite (fp);

  fclose (fp);

  return  0;
}

//  Set  Stop [0 .. 6]  1  or  0   according to whether
//  X [1 .. T] has a stop codon in the corresponding reading frame.
//  Stop [6]  is always set  0 .
void  Find_Stop_Codons  (char X [], int T, int Stop [])
{
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
  for  (i = 1;  i <= T;  i ++)
  {
    switch  (tolower (X [i]))
    {
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
      Stop [i % 3] = 1;
    else if  (State > 7)
      Stop [3 + i % 3] = 1;
  }

  return;
}

/* Read next string from  fp  (assuming one string per line with each string
*  preceded by a name tag) into  T [1 ..]
*  which has  Size  characters.  Allocate extra memory if needed
*  and adjust  Size  accordingly.  Return  1  if successful,  0
*  otherwise (e.g., EOF). */
int  Read_String  (FILE * fp, char * & T, long int & Size, char Name [])
{

  long int  Len;
  int  Ch;

  if  (fscanf (fp, "%s", Name) == EOF)
    return  0;

  while  ((Ch = fgetc (fp)) != EOF && (Ch == ' ' || Ch == '\t'))
    ;
  if  (Ch == EOF)
    return  0;

  ungetc (Ch, fp);
  T [0] = '\0';
  Len = 1;
  while  ((Ch = fgetc (fp)) != EOF && Ch != '\n')
  {
    if  (isspace (Ch))
      continue;

    if  (Len >= Size)
    {
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

  return  1;
}
