//
//   A. L. Delcher
//
//     File:  ~delcher/TIGR/strarray.h
//  Version:  1.01a  4 Dec 97
//            Fixed divide-by-zero bug
//  Version:  1.02  25 Feb 98
//            Removed usused variables.
//
//    Copyright (c) 1997 by Arthur Delcher, Steven Salzberg, Simon
//    Kasif, and Owen White.  All rights reserved.  Redistribution
//    is not permitted without the express written permission of
//    the authors.
//
//  Definitions for an array class whose index is strings.
//
// Hacked by T. Schiex

#ifndef  STRARRAY_H_INCLUDED
#define  STRARRAY_H_INCLUDED

#include <math.h>
#include <string.h>
#include <ctype.h>

const long int  INCR_SIZE = 10000;
const long int  INIT_SIZE = 10000;
const double  MINIMUM_DELTA_VAL = 1e-5;
const double  MINIMUM_COUNT = 0.5;
const double  MIN_CORRECT_PREDICTS = 15.0;
const int  SAMPLE_SIZE_BOUND = 400;
const double  SAMPLE_WEIGHT = 400.0;
const double  SIGNIFICANCE_THRESHOLD = 0.01;
const double  Z_SCORE_THRESHOLD = 2.33;
const double  CHI2_THRESHOLD = 11.3;    // 0.01 significance  v=3
const double  BINARY_PRECISION = 1.0/65535.0;

float  New_Analyze  (float [], float [], int, float []);
char  Complement  (char);
void  Reverse_Complement  (char [], long int);

class  String_Array {
 private:
  float  * Val;
  int  Max_Str_Len;
  int  Alphabet_Size;
  int  Num_Entries;
  
  int  String_To_Sub  (char [], int);
  int  AntiString_To_Sub  (char [], int);
  
 public:
  String_Array  ();
  String_Array  (int, int);
  ~ String_Array  ();
  void  Condense  (String_Array &);
  void  Incr  (char [], int, double);
  double  Min_Val  ();
  void  Normalize  ();
  void  Set  (double = 0.0);
  void  Set_Lambda  (String_Array &);
  int  Value  (char);
  int   AntiValue  (char);    
  void  BinWrite  (FILE *);   
  double  operator []  (int i)    {  return  Val [i];    }
  float &  operator ()  (char [], int);
  float & Anti (char [], int);
};

//  Default constructor.
String_Array :: String_Array  () {
  Val = NULL;
  Num_Entries = Max_Str_Len = Alphabet_Size = 0;
}

//  Construct a  String_Array  that can be indexed by strings of length
//  up to  L  over an alphabet of size  S .
String_Array :: String_Array  (int L, int S) {
  int  i;
  
  myassert (L > 0 && S > 0);
  
  Max_Str_Len = L;
  Alphabet_Size = S;
  
  Num_Entries = int (pow (S, L + 1) - 1) / (S - 1);
  Val = (float *) Safe_malloc (Num_Entries * sizeof (float));
  
  Val [0] = 1.0;
  for  (i = 1;  i < Num_Entries;  i ++)
    Val [i] = 0.0;
}

//  Destroy this  String_Array  by freeing its memory.
String_Array :: ~ String_Array  ()  {
  free (Val);
}

//  Reset the delta values for  this  string array to be
//  used as a simple model that is equivalent to the context model
//  using the lambda values in  L .
void  String_Array :: Condense  (String_Array & L)  {
  long int  i, Ct, Prev_i, Prev_i_Start, Prev_Size;
  long int  Size, Sub;
  
  myassert (Max_Str_Len == L . Max_Str_Len + 1);
  
  Val [0] = 1.0;
  
  Ct = Prev_i = Prev_i_Start = Sub = 0;
  Prev_Size = 1;
  Size = Alphabet_Size;
  
  for  (i = 1;  i < Num_Entries;  i ++) {
    Val [i] = L . Val [Sub] * Val [i]
      + (1.0 - L . Val [Sub]) * Val [Prev_i_Start + Prev_i];
    
    if  (++ Ct == Size) {
      Ct = 0;
      Prev_i_Start += Prev_Size;
      Prev_Size  = Size;
      Size *= Alphabet_Size;
      Prev_i = 0;
    }
    else
      Prev_i = (Prev_i + 1) % Prev_Size;
    
    if  (i % Alphabet_Size == 0)
      Sub ++;
  }
   
  return;
}

//  Add  D  to the entry at the subscript indicated by string
//  S [0 .. L-1]  in this  String_Array .
void  String_Array :: Incr  (char S [], int L, double D) {
  int  i;
  
   i = String_To_Sub (S, L);
   Val [i] += D;
   
   return;
}

//  Return the minimum value in this  String_Array .
double  String_Array :: Min_Val  () {
  double  Min;
  int  i;
  
  Min = Val [0];
  
  for  (i = 1;  i < Num_Entries;  i ++) {
    if  (Val [i] < Min)
      Min = Val [i];
  }
  
  return  Min;
}

//  Convert the values in this  String_Array  to probabilities.
void  String_Array :: Normalize  () {
  double  Sum;
  long int  j, Ct, Prev_Size, Prev_Start, Prev_Sub, Size, Start;
  
  Val [0] = 1.0;
  
  Ct = 0;
  Prev_Size = 1;
  Prev_Start = 0;
  Prev_Sub = 0;
  Size = Alphabet_Size;
  
  for  (Start = 1;  Start < Num_Entries;  Start += Alphabet_Size)  {
    if  (Prev_Size > 1)
      for  (j = 0;  j < Alphabet_Size;  j++) {
	//	 printf("P %d %f %f\n", Start+j,Val [Start + j],Val [Prev_Start + Prev_Sub + j]);
	Val [Start + j] += Val [Prev_Start + Prev_Sub + j];
      }
    
    Sum = 0.0;
    for  (j = 0;  j < Alphabet_Size;  j++) {
      Sum += Val [Start + j];
      //       printf("A %d %d %f\n", Start,j,Val [Start + j]);
    }
    
    if  (Sum == 0.0)  {
      for  (j = 0;  j < Alphabet_Size;  j++)
	Val [Start + j] = 0.0;
    }
    else  {
      for  (j = 0;  j < Alphabet_Size;  j++) {
	Val [Start + j] /= Sum;
	if  (Val [Start + j] < MINIMUM_DELTA_VAL)
	  Val [Start + j] = MINIMUM_DELTA_VAL;
      }
    }
    // Bug ?
    //    if  ((Ct += Alphabet_Size) == Size) {
    if  (++ Ct == Size) {
      Ct = 0;
      Prev_Start += Prev_Size;
      Prev_Size = Size;
      Size *= Alphabet_Size;
      Prev_Sub = 0;
    }
    else
      Prev_Sub = (Prev_Sub + Alphabet_Size) % Prev_Size;
  }
  return;
}


//  Set the value of all elements in this  String_Array  to  X .
void  String_Array :: Set  (double X) {
  int  i;
  
  Val [0] = 1.0;
  for  (i = 1;  i < Num_Entries;  i ++)
    Val [i] = X;
  
  return;
}

//  Estimate lambda values for  this  string array from
//  delta counts in  D .
void  String_Array :: Set_Lambda  (String_Array & D) {
  float  * Prob, Sum;
  int  Len;
  long int  Prev_Start, Prev_Sub, Size, Sub;
  long int  i, Ct, Prev_Size, Used;
  
  myassert (Max_Str_Len == D . Max_Str_Len - 1);
  
  //   printf ("Context Length Usage\n");
  // printf (" Len      Used\n");
  Len = 1;
  Used = 0;
  Val [0] = 1.0;
  
  Prob = (float *) Safe_malloc (D . Num_Entries * sizeof (float));
  Sum = 0.0;
  for  (i = 1;  i <= Alphabet_Size;  i ++)
    Sum += D . Val [i];
  for  (i = 1;  i <= Alphabet_Size;  i ++)
    Prob [i] = D . Val [i] / Sum;
  
  Sub = Alphabet_Size + 1;
  Prev_Start = 1;
  Prev_Sub = 0;
  Ct = 0;
  Prev_Size = 1;
  Size = Alphabet_Size;
  for  (i = 1;  i < Num_Entries;  i ++)  {
    Val [i] = New_Analyze (D . Val + Sub, Prob + Prev_Start + Prev_Sub,
			   Alphabet_Size, Prob + Sub);
    if  (Val [i] > 0.0)
      Used ++;
    
    if  (++ Ct == Size)  {
      //           printf ("%3d %9ld (%3d%%)\n", Len ++, Used,
      //                           int (0.5 + (100.0 * Used) / Size));
      Used = 0;
      Ct = 0;
      Prev_Start += Size;
      Prev_Size *= Alphabet_Size;
      Size *= Alphabet_Size;
      Prev_Sub = 0;
    }  else  {
      Prev_Sub = (Prev_Sub + Alphabet_Size) % Size;
    }
    
    Sub += Alphabet_Size;
  }
  
  free (Prob);
  
  return;
}

//  Convert string  S [0 .. L-1]  to a subscript in the  Val  array.
int  String_Array :: String_To_Sub  (char S [], int L)  {
  int  i, Sub;
  
  myassert (L <= Max_Str_Len);
  
  if  (L == 0)
    return  0;
  
  Sub = 0;
  for  (i = 0;  i < L;  i ++)  
    Sub = Sub * Alphabet_Size + Value (S [i]);
  
  Sub += int (pow (Alphabet_Size, L) - 1) / (Alphabet_Size - 1);
  
  return  Sub;
}

//  Convert the antiparallel of string S [0 .. L-1] to a subscript in
//  the Val array.
int  String_Array :: AntiString_To_Sub  (char S [], int L)
{
  int  i, Sub;
  
  myassert (L <= Max_Str_Len);
  
  if  (L == 0)
    return  0;
  
  Sub = 0;
  for  (i = L-1;  i >=0;  i --)
    Sub = Sub * Alphabet_Size + AntiValue (S [i]);
  
  Sub += int (pow (Alphabet_Size, L) - 1) / (Alphabet_Size - 1);
  
  return  Sub;
}


//  Return the numeric value of character  Ch .
int  String_Array :: Value  (char Ch)  {
  switch  (tolower (Ch)) {
  case  'a' :
    return  0;
  case  'c' :
    return  1;
  case  'g' :
    return  2;
  case  't' :
  case  'u' :   
    return  3;
  default :
    fprintf (stderr, "ERROR:  Unexpected character \'%c\' (ASCII %d)\n",
	     Ch, int (Ch));
    exit (-1);
  }
}

//  Return the numeric value of character  Ch .
//  the complement char is used
int  String_Array :: AntiValue  (char Ch) {
  switch  (tolower (Ch))
    {
    case  'u' :   
    case  't' :
      return  0;
    case  'g' :
      return  1;
    case  'c' :
      return  2;
    case  'a' :
      return  3;
    default :
      fprintf (stderr, "ERROR:  Unexpected character \'%c\' (ASCII %d)\n",
	       Ch, int (Ch));
      exit (-1);
    }
}

// Return a reference to the value associated with the antiparallel
// of S [0 .. L-1]   
float & String_Array :: Anti (char S [], int L)  {
  int  i;
  
  i = AntiString_To_Sub (S, L);
  return  Val [i];
}

//  Return a reference to the value associated with string  S [0 .. L-1] .
float &  String_Array :: operator ()  (char S [], int L) {
  int  i;
  
   i = String_To_Sub (S, L);
   return  Val [i];
}

// Writinga binary version on the model. Does not care of Endianness
// of the architecture (taken care when read)
void  String_Array :: BinWrite  (FILE *fp) {
  int    i;
  unsigned short  buffer;
  
  (void) fwrite(&Max_Str_Len,   sizeof(int), 1, fp);
  (void) fwrite(&Alphabet_Size, sizeof(int), 1, fp);
  (void) fwrite(&Num_Entries,   sizeof(int), 1, fp);
  
  for (i = 0;  i < Num_Entries;  i ++) {
    buffer = (unsigned short) floor((double) (Val[i]/BINARY_PRECISION + 0.5));
    (void) fwrite(&buffer,  sizeof(unsigned short), 1, fp);
  }
  
  return;
}

//  Return the lambda value for the frequency counts in  A [0 .. N-1]
//  corresponding to the longer context string compared to the predicted
//  probabilities in  P [0 .. N-1]  for the next shorter suffix context string.
//  Also set  C [0 .. N-1]  to the probabilities for this longer context
//  string.
float  New_Analyze  (float A [], float P [], int N, float C []) {
  const int  CHI2_ENTRIES = 7;
  float  Chi2_Val [CHI2_ENTRIES] = {2.37, 4.11, 6.25, 7.81, 9.35, 11.3, 12.8};
  float  Significance [CHI2_ENTRIES]
    = {0.50, 0.75, 0.90, 0.95, 0.975, 0.99, 0.995};
  double  A_Sum, Chi2, E, Lambda;
  int  i;
  
  A_Sum = 0.0;
  for  (i = 0;  i < N;  i ++)
    A_Sum += A [i];
  if  (A_Sum >= SAMPLE_SIZE_BOUND)  {
    for  (i = 0;  i < N;  i ++)
      C [i] = (A [i] + P [i]) / (A_Sum + 1.0);
    return  1.0;
  }
  
  Chi2 = 0.0;
  for  (i = 0;  i < N;  i ++)  {
    E = A_Sum * P [i];
    if  (E != 0.0)
      Chi2 += pow (A [i] - E, 2.0) / E;
  }
  for  (i = 0;  i < CHI2_ENTRIES && Chi2_Val [i] < Chi2;  i ++)
    ;
  if  (i == 0)
    Lambda = 0.0;
  else if  (i == CHI2_ENTRIES)
    Lambda = 1.0;
  else
    Lambda = Significance [i - 1]
      + ((Chi2 - Chi2_Val [i - 1]) / (Chi2_Val [i] - Chi2_Val [i - 1]))
      * (Significance [i] - Significance [i - 1]);
  Lambda *= A_Sum / SAMPLE_WEIGHT;
  if  (Lambda > 1.0)
    Lambda = 1.0;
  
  for  (i = 0;  i < N;  i ++)
    C [i] = Lambda * ((A [i] + P [i]) / (A_Sum + 1.0))
      + (1.0 - Lambda) * P [i];
  
  return  Lambda;
}


/* Returns the DNA complement of  Ch . */
char  Complement  (char Ch)  {
  switch  (tolower (Ch)) {
  case  'a' :
    return  't';
  case  'c' :
    return  'g';
  case  'g' :
    return  'c';
  case  't' :
    return  'a';
  case  'r' :          // a or g
    return  'y';
  case  'y' :          // c or t
    return  'r';
  case  's' :          // c or g
    return  's';
  case  'w' :          // a or t
    return  'w';
  case  'm' :          // a or c
    return  'k';
  case  'k' :          // g or t
    return  'm';
  case  'b' :          // c, g or t
    return  'v';
  case  'd' :          // a, g or t
    return  'h';
  case  'h' :          // a, c or t
    return  'd';
  case  'v' :          // a, c or g
    return  'b';
  default :            // anything
    return  'n';
  }
}

//  Set  S [1 .. T]  to its DNA reverse complement.
void  Reverse_Complement  (char S [], long int T)  {
  char  Ch;
  long int  i, j;
  
  for  (i = 1, j = T;  i < j;  i ++, j --) {
    Ch = S [j];
    S [j] = Complement (S [i]);
    S [i] = Complement (Ch);
  }
  
  if  (i == j)
    S [i] = Complement (S [i]);
}

#endif
