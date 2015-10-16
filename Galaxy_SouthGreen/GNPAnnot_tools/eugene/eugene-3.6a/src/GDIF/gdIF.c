/* ------------------------------------------------------------------
   Copyright (C) 2004 INRA <eugene@ossau.toulouse.inra.fr>
  
   This program is open source; you can redistribute it and/or modify
   it under the terms of the Artistic License (see LICENSE file).
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
  
   You should have received a copy of Artistic License along with
   this program; if not, please see http://www.opensource.org
  
   $Id: gdIF.c,v 1.19 2005-11-22 19:26:50 tschiex Exp $
   ------------------------------------------------------------------
   File:     gdIF.c
   Contents: Functions for interfacing with libGD
   ------------------------------------------------------------------*/

#include <stdlib.h>
#include <math.h>
#include <gd.h>
#include <gdfontt.h>
#include <gdfonts.h>
#include <gdfontmb.h>
#include <gdfontl.h>
#include <gdfontg.h>
#ifdef __sun__
#include <strings.h>
#endif
#ifdef __linux__
#include <string.h>
#endif

#define LMargin 20
#define RMargin 20
#define TMargin 20
#define BMargin 20

#define NLigne  4

struct Image
{
  gdImagePtr im;
  int Left;
  int Right;
};

struct Image *images;

int rx,ry,Step;
int Len,From,To,Offset;
int ImgLen;
int OvLap;
int NbIm;
int Col[16];

char* Ylab[9] = {"IR","-3","-2","-1","IG","+1","+2","+3","IF"};
char *FName;
char TName[1024];

#define Ystretch 1.2

#define ToX(nuc) ((int)(((nuc)*rx)/ImgLen+LMargin))
#define ToY(ph,pos) ((int)(((((NLigne-(double)(ph))*Ystretch) - (pos)+1.0)/((2*NLigne+1)*Ystretch))*ry+TMargin))


void InitImg(struct Image *image, int n)
{
  int i;
  char str[12];

  image->im = gdImageCreate(rx+LMargin+RMargin,ry+TMargin+BMargin);
  
  Col[0] = gdImageColorAllocate(image->im, 255, 255, 255); /* white      */
  Col[1] = gdImageColorAllocate(image->im, 255,   0,   0); /* red        */
  Col[2] = gdImageColorAllocate(image->im,   0,   0, 255); /* blue       */
  Col[3] = gdImageColorAllocate(image->im,   0,   0,   0); /* black      */
  Col[4] = gdImageColorAllocate(image->im, 255,   0, 255); /* violet     */
  Col[5] = gdImageColorAllocate(image->im, 120, 255, 210); /* vert pale  */
  Col[6] = gdImageColorAllocate(image->im, 180, 180, 180); /* gris 1     */
  Col[7] = gdImageColorAllocate(image->im, 200, 200, 200); /* gris 2     */
  Col[8] = gdImageColorAllocate(image->im, 220, 220, 220); /* gris 3     */
  Col[9] = gdImageColorAllocate(image->im, 255, 165,   0); /* orange     */
  Col[10]= gdImageColorAllocate(image->im, 244, 164,  96); /* sandybrown */
  Col[11]= gdImageColorAllocate(image->im,   0, 255,   0); /* vert       */
  
  gdImageColorTransparent(image->im, Col[0]);
  
  gdImageLine(image->im, LMargin-1,  TMargin-1, LMargin-1,
	      ry+TMargin+5, Col[3]);
  gdImageLine(image->im, rx+LMargin, TMargin-1, rx+LMargin,
	      ry+TMargin+5, Col[3]);
  
  for (i = -NLigne; i <= NLigne; i++) {
    gdImageString(image->im,gdFontTiny,5, ToY(i,0.6),(unsigned char *)Ylab[i+4],Col[3]);
    gdImageString(image->im,gdFontTiny,rx+LMargin+5, ToY(i,0.6),(unsigned char *)Ylab[i+4],Col[3]);
  }
  
  image->Left = n*(ImgLen-OvLap)+From;
  image->Right = image->Left+ImgLen;
  
  sprintf(str,"%d",image->Right+Offset);
  gdImageString(image->im,gdFontTiny,
		rx+LMargin - (strlen(str)*gdFontTiny->w/2), ry+TMargin+9,(unsigned char *)str,Col[3]);
  gdImageString(image->im,gdFontTiny,
		rx+LMargin - (strlen(str)*gdFontTiny->w/2), TMargin-15,(unsigned char *)str,Col[3]);
  
  sprintf(str,"%d",image->Left+Offset);
  gdImageString(image->im,gdFontTiny,
		LMargin - gdFontTiny->w/2, ry+TMargin+9,(unsigned char *)str,Col[3]);
  gdImageString(image->im,gdFontTiny,
		LMargin - gdFontTiny->w/2, TMargin-15,(unsigned char *)str,Col[3]);
  
  
  for (i = Step; i<= ImgLen-Step; i+= Step) {
    sprintf(str,"%i",i+image->Left+Offset);
    gdImageString(image->im,gdFontTiny,
		  ToX(i) - (strlen(str)*gdFontTiny->w/2), ry+TMargin+9,(unsigned char *)str,Col[3]);
    gdImageString(image->im,gdFontTiny,
		  ToX(i) - (strlen(str)*gdFontTiny->w/2), TMargin-15,(unsigned char *)str,Col[3]);
    gdImageLine(image->im,ToX(i),ry+TMargin,ToX(i),ry+TMargin+5,Col[3]);
    gdImageLine(image->im,ToX(i),TMargin,ToX(i),TMargin-5,Col[3]);
  }
}

/* Initilisation des donnees images : resolution... */
void InitPNG(int x,int y, int offset, int DFrom, int DTo, int olap, int ILen,char *name)
{
  int i,lstep;

  FName = name;
  
  rx = x;
  ry = y;
  Len = DTo-DFrom+1;
  From = DFrom;
  To = DTo;
  ImgLen = ILen;
  Offset = offset;
  
  Step = ImgLen/(rx/80);
  lstep = (int)pow(10,(int)floor(log10(Step)));
  Step = Step/lstep;

  if (Step > 7) Step = 10;
  if (Step == 6) Step = 5;
  if (Step == 3) Step = 2;
  if (Step == 4) Step = 5;
    
  Step *= lstep;

  OvLap = ((olap <0) ? Step : olap);
  NbIm = 1+(Len-OvLap-1)/(ImgLen-OvLap);
  
  images = (struct Image *)malloc(sizeof(struct Image)*NbIm);

  for (i = 0; i< NbIm; i++)
    InitImg(&images[i],i);
}

void SaveClosePNG(struct Image *image, int num)
{
  FILE *out;
  char str[12];

  gdImageInterlace(image->im, 1);
  strcpy(TName,FName);
  sprintf(str,".%03d",num);
  strcat(TName,str);
#if GIF
  strcat(TName,".gif");
  out = fopen(TName,"wb");
  gdImageGif(image->im, out);
#else
  strcat(TName,".png");
  out = fopen(TName,"wb");
  gdImagePng(image->im, out);
#endif
  fclose(out);
  gdImageDestroy(image->im);
  image->im = NULL;
  image->Left = -1;
  image->Right = -1;
}


int CheckIn(unsigned int nuc, struct Image *image)
{
  return ((nuc >= image->Left) && (nuc < image->Right));
}

void PlotBarI(unsigned int nuc, signed char phase, double pos, int width, int col)
{
  int i;

  for (i=0; i< NbIm; i++)
    if (CheckIn(nuc,&images[i]))
      gdImageLine(images[i].im,ToX(nuc-images[i].Left),ToY(phase, pos)-width,ToX(nuc-images[i].Left),ToY(phase, pos)+width, Col[col]);
}

void PlotBarF(unsigned int nuc, signed char phase, double pos, double width, int col)
{
  int i;

  for (i=0; i< NbIm; i++)
    if (CheckIn(nuc,&images[i]))
      gdImageLine(images[i].im,ToX(nuc-images[i].Left),ToY(phase, pos-(width/2)),ToX(nuc-images[i].Left),ToY(phase, pos+(width/2)), Col[col]);
}

void PlotLine(unsigned int nuc1,  unsigned int nuc2,
	      signed char phase1, signed char phase2,
	      double pos1,        double pos2, int col)
{
  int i;

  for (i=0; i< NbIm; i++)
    if ((CheckIn(nuc1,&images[i]))&& (CheckIn(nuc2,&images[i])))
      gdImageLine(images[i].im,ToX(nuc1-images[i].Left),ToY(phase1, pos1),ToX(nuc2-images[i].Left),ToY(phase2, pos2), Col[col]);
}

void PlotString(unsigned int nuc, signed char phase, double pos, char st[], int col)
{
  int i;
  for (i=0; i< NbIm; i++)
    if (CheckIn(nuc, &images[i]))
      gdImageString(images[i].im, gdFontTiny,
		    (ToX(nuc-images[i].Left)-(strlen(st)*gdFontTiny->w/2)),
		    ToY(phase, pos), (unsigned char *)st, Col[col]);
}

void ClosePNG()
{
  int i;

  for (i=0; i< NbIm; i++)
    SaveClosePNG(&images[i],i);
}

// Parametre = 1 -> toute 1ere image de la 1ere sequence traitée
// Parametre = 0 -> la suite....*.fasta
void OutputHTMLFileNames(int firstImg, char* html_dir, FILE* OUT)
{
  int i;
  char str[12];
  char *basename;

  basename = rindex(FName,'/');
  if (basename == NULL) basename = FName;
  else basename++;
  
  if (firstImg) {
    strcpy(TName,basename);
    sprintf(str,".%03d",0);
    strcat(TName,str);
#if GIF
    strcat(TName,".gif");
#else
    strcat(TName,".png");
#endif
    fprintf(OUT, "				  <img src=\"%s\" name=\"show\"\n",
	    TName);
    fprintf(OUT, "				       width=\"%d\" height=\"%d\">\n",
	    rx+40, ry+40);
    fprintf(OUT, "				</td>\n"
	    "			      </tr>\n"
	    "			      <tr>\n"
	    "				<td width=\"100\" align=\"center\"\n"
	    "				    bgcolor=\"#c0dbe2\">\n"
	    "				  <img src=\"%s/Images/first.jpg\"\n"
	    "				       onclick=\"first();\"\n"
	    "				       title=\"Jump to beginning\" "
	    "align=\"middle\">\n"
	    "				  &nbsp; &nbsp; &nbsp;\n"
	    "				  <img src=\"%s/Images/back.jpg\"\n"
	    "				       onclick=\"previous();\"\n"
	    "				       title=\"Back\" "
	    "align=\"middle\">\n"
	    "				</td>\n"
	    "				<td align=\"center\" "
	    "bgcolor=\"#9bc7d0\">\n"
	    "				  <select name=\"slide\" "
	    "onChange=\"change();\" size=\"1\">\n", html_dir, html_dir);
    fprintf(OUT, "				    <option value=\"%s\"\n",TName);
    fprintf(OUT, "					    selected>%s : %d -> %d"
	    "</option>\n",TName,From,ImgLen+From);
  }
 
  for (i=firstImg; i< NbIm; i++) {
    strcpy(TName,basename);
    sprintf(str,".%03d",i);
    strcat(TName,str);
#if GIF
    strcat(TName,".gif");
#else
    strcat(TName,".png");
#endif
    fprintf(OUT, "                                    <option value=\"%s\">",TName);
    if (i==0)
      fprintf(OUT, "%s : %d -> %d</option>\n",
	      TName, ImgLen*i+From,
	      (ImgLen*(i+1))+From < Len ? (ImgLen*(i+1))+From : Len);
    else
      fprintf(OUT, "%s : %d -> %d</option>\n",
	      TName, ImgLen*i+From-(OvLap*i),
	      (ImgLen*(i+1))+From-(OvLap*i) < Len ? (ImgLen*(i+1))+From-(OvLap*i) :Len);
  }
}
