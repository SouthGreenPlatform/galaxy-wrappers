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
// $Id: System.cc,v 1.15 2006-06-02 08:44:20 tschiex Exp $
// ------------------------------------------------------------------
// File:     System.cc
// Contents: utilitary functions
// ------------------------------------------------------------------


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif
#include <time.h>
#include <string>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/times.h>

#if defined (__sun)
#include <ieeefp.h>
int isinf(double x) { return !finite(x) && x==x; }
#endif

#include "Const.h"
#include "System.h"


// ------------------------------------------------------------------
// BASENAME: returns a pointer to the filename, w/o any
// leading prefix
// ------------------------------------------------------------------
char * BaseName(char *path)
{
  char *lead = rindex(path,'/');
  
  return ((lead != NULL)  ? lead+1 : path);
}
// ------------------------------------------------------------------
// ProbeFile: test if a file exists and has non empty size
// return 1 if Ok.
// ------------------------------------------------------------------
int ProbeFile(const char *defdir, const char *filename)
{
  struct stat FileStat;

  if ((stat(filename,&FileStat) == 0) && S_ISREG(FileStat.st_mode))
	return (FileStat.st_size != 0);

  if (defdir) {
    char buffer[FILENAME_MAX];
    strcat(strcat(strcpy(buffer, defdir), "/"), filename);
    if ((stat(buffer,&FileStat) == 0) && S_ISREG(FileStat.st_mode))
	return (FileStat.st_size != 0);
  }
  return 0;
}
// ------------------------------------------------------------------
// Function to open a file by using environment variables.  We first
// search in local directory and if not found in directory defdir
// ------------------------------------------------------------------

FILE *FileOpen(const char *defdir, const char *filename,
	       const char *mode,   int sloppy)
{
  FILE  *fp;
  char buffer[FILENAME_MAX];
  
  if ((fp = fopen(filename, mode)))
    return fp;
  
  if (defdir) {
    strcat(strcat(strcpy(buffer, defdir), "/"), filename);
    fp = fopen (buffer, mode);
  }
  
  if  (!fp) {
    if (sloppy)
      fprintf (stderr, "WARNING: Could not open file %s ", filename);
    else {
      fprintf (stderr, "ERROR: Could not open file %s \n", filename);
      exit (1);
    }
  }
  return  fp;
}
// ------------------------------------------------------------------
// Malloc sur
// ------------------------------------------------------------------

void *  Safe_malloc  (size_t Len)     
{
  void  * P;
  
  
  if  ((P = malloc (Len)) == NULL)
    {
      perror("malloc");
      exit (1);
    }
  return  P;
}

// ------------------------------------------------------------------
// Realloc sur
// ------------------------------------------------------------------

void *  Safe_realloc  (void * Q, size_t Len)
{
  void  * P;
  
  
  if  ((P = realloc (Q, Len)) == NULL)
    {
      perror("realloc");
      exit (1);
    }  
  return  P;
}


// ------------------------------------------------------------------
// Provides the current date on the form jjMMMaaaa 
// where jj, aaaa are numbers and MMM caracters
// ------------------------------------------------------------------
void GetStrDate (char* d)
{

  char *m,*j,*a;
  m = new char[4]; j = new char[3]; a = new char[5];
  time_t t=time(0);

  strcpy(d,ctime(&t));
  sscanf(d, "%*s %s %s %*s %s", m,j,a);
  strcpy(d,j);strcat(d,m);strcat(d,a);

  delete [] m; delete [] j; delete [] a;
}
// --------------------------------------------------------------------
// Timer management functions
// -------------------------------------------------------------------- 
double cpuTime()
{
  static struct tms buf;
  
  times(&buf);
  return ((double) (buf.tms_utime+buf.tms_stime+buf.tms_cutime+buf.tms_cstime))
    / ((double) sysconf(_SC_CLK_TCK));
}
