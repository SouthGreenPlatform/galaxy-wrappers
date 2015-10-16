/* ------------------------------------------------------------------
   Copyright (C) 2004 INRA <eugene@ossau.toulouse.inra.fr>
  
   This program is open source; you can redistribute it and/or modify
   it under the terms of the Artistic License (see LICENSE file).
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
  
   You should have received a copy of Artistic License along with
   this program; if not, please see http://www.opensource.org
  
   $Id: yacc.y,v 1.7 2005-06-01 15:44:23 tschiex Exp $
   ------------------------------------------------------------------
   File:     yacc.y
   Contents: analyseur syntaxique du langage utilisateur eugene
   ------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
/*----------- DECLARATION et DEFINITION des unites lexicales -------------*/
/*------------------------------------------------------------------------*/
%{
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#include "../../Const.h"
#include "../../SensorIF.h"
#include "../../System.h"

#include "structure.h"


//---------------------------------------------------------------------------------
  // La prediction (pour verifier la coherence)
  char* Choice;

  // les listes d'info utilisateur
  ptUTIL SignalUser,ContentsUser;

  // pour afficher une raison (source)
  char* raison;

  // declarations externes
  extern int yylex();
  extern char *yytext;
  extern FILE *yyin;

  /* nouveauté Linux compatible Solaris */
  void yyerror(char* s) {
    fprintf (stderr, "%s : \"%s\"\n",s,yytext);
    exit(0);
  }

//---------------------------------------------------------------------------------
  char* ContentsText[13] = {
    "exon f1","exon f2","exon f3",
    "exon r1","exon r2","exon r3", 
    "intron f", "intron r",
    "intergenic",
    "utr5 f","utr5 r",
    "utr3 f","utr3 r"
  };

  char* SignalText[8] = {
    "start f", "start r",
    "stop f","stop r",
    "acceptor f","acceptor r",
    "donor f", "donor r"
  };

 
//---------------------------------------------------------------------------------
// Execution de la regle
// les aretes "contents" sont passees par une exponentielle
// les aretes "signal" sont bornees entre 0 et 1
//---------------------------------------------------------------------------------
void ExecuteRegle(ptUTIL ut, DATA *d)
{
  if (ut->tab <= 12) // contents edges
    d->contents[ut->tab] += ut->delta;
  else
    switch (ut->tab) {
    case 13:
      d->sig[DATA::Start].weight[Signal::Forward] = log(Max(0.0,Min(1.0,ut->delta)));
      d->sig[DATA::Start].weight[Signal::ForwardNo] = log(1.0-Max(0.0,Min(1.0,ut->delta)));
      break;
    case 14:
      d->sig[DATA::Start].weight[Signal::Reverse] = log(Max(0.0,Min(1.0,ut->delta)));
      d->sig[DATA::Start].weight[Signal::ReverseNo] = log(1.0-Max(0.0,Min(1.0,ut->delta)));
      break;
    case 15:
      d->sig[DATA::Stop].weight[Signal::Forward]  = log(Max(0.0,Min(1.0,ut->delta)));
      d->sig[DATA::Stop].weight[Signal::ForwardNo]  = log(1.0-Max(0.0,Min(1.0,ut->delta)));
      break;
    case 16:
      d->sig[DATA::Stop].weight[Signal::Reverse]  = log(Max(0.0,Min(1.0,ut->delta)));
      d->sig[DATA::Stop].weight[Signal::ReverseNo]  = log(1.0-Max(0.0,Min(1.0,ut->delta)));
      break;
    case 17:
      d->sig[DATA::Acc].weight[Signal::Forward]   = log(Max(0.0,Min(1.0,ut->delta)));
      d->sig[DATA::Acc].weight[Signal::ForwardNo]   = log(1.0-Max(0.0,Min(1.0,ut->delta)));
      break;
    case 18:
      d->sig[DATA::Acc].weight[Signal::Reverse]   = log(Max(0.0,Min(1.0,ut->delta)));
      d->sig[DATA::Acc].weight[Signal::ReverseNo]   = log(1.0-Max(0.0,Min(1.0,ut->delta)));
      break;
    case 19:
      d->sig[DATA::Don].weight[Signal::Forward]   = log(Max(0.0,Min(1.0,ut->delta)));
      d->sig[DATA::Don].weight[Signal::ForwardNo]   = log(1.0-Max(0.0,Min(1.0,ut->delta)));
      break;
    case 20:
      d->sig[DATA::Don].weight[Signal::Reverse]   = log(Max(0.0,Min(1.0,ut->delta)));
      d->sig[DATA::Don].weight[Signal::ReverseNo]   = log(1.0-Max(0.0,Min(1.0,ut->delta)));
      break;
    }
}
 
//---------------------------------------------------------------------------------
// Application des regles de ut valables a la position i
// Initialement, ut a toujours un Bloc d'info inutile en header
// c'est donc ut-> suiv le bloc a traiter.
//---------------------------------------------------------------------------------
void Util(int i, ptUTIL ut, DATA *d)
{
  // liste epuisee ou elements trop avances: on sort
  if ((ut == NULL) || (ut->suiv == NULL) ||  (ut->suiv->n1-1) > i) return;
  
  //element actif, on agit !
  if ((i >= (ut->suiv->n1-1)) && (i <= (ut->suiv->n2-1)))
    ExecuteRegle(ut->suiv, d); 
  
  // on avance
  Util(i, ut->suiv, d);
}
//---------------------------------------------------------------------------------
// Ecrire un util au format texte dans un flot
//---------------------------------------------------------------------------------
void WriteContents(ptUTIL ut,FILE *flot)
{
  fprintf(flot,"%s ",ContentsText[ut->tab]);
  fprintf(flot,"[%d..%d] ",ut->n1,ut->n2);
  if (isnan(ut->delta))  fprintf(flot,"infinity\n");
  else fprintf(flot,"%a\n",ut->delta);
}
//---------------------------------------------------------------------------------
void WriteSignal(ptUTIL ut,FILE *flot)
{
  fprintf(flot,"%s ",SignalText[ut->tab-13]);
  fprintf(flot,"%d ",ut->n1);
  fprintf(flot,"%a\n",ut->delta);
}
//---------------------------------------------------------------------------------
void WriteUtils(ptUTIL ut,FILE *flot)
{
  if (ut->suiv == NULL) return;
  
  if (ut->suiv->tab <= 12) // contents edge
    WriteContents(ut->suiv,flot);
  else WriteSignal(ut->suiv,flot);

  WriteUtils(ut->suiv,flot);
}

//---------------------------------------------------------------------------------
// Insere a a la bonne position dans r, retourne la position suivante
//---------------------------------------------------------------------------------
ptUTIL suivant(ptUTIL a,ptUTIL r)
{
  static ptUTIL s ;
  
  // au bout... on insere
  if (r->suiv == NULL) {
    r->suiv = a; 
    s = NULL;
  }
  else if (a->n1 >= r->suiv->n1) // pas assez loin.
    suivant(a,r->suiv); // on avance 
  else { // on est assez loin, on insere
    s = r->suiv;
    r->suiv = a;
  }
  return s;
}
//---------------------------------------------------------------------------------
// Interface malloc
//---------------------------------------------------------------------------------
ptUTIL mall()
{
  ptUTIL pu;
  pu = (ptUTIL)malloc(sizeof(UTIL)) ;
  if (pu!=NULL) {return pu;}
  else {fprintf (stderr, "Allocation failed\n"); exit(1);};
}
//---------------------------------------------------------------------------------
// Message de verification info (Signal)
//---------------------------------------------------------------------------------
void MessageSignal(ptUTIL ut,int cond,char* m)
{
  if (ut->rais!=raison) 
    {raison = ut->rais;printf("Source: %s\n",raison);}; 
  if (cond)       
    {printf("Info. %s at %d consistent\n",m,ut->n1);}
  else 
    {printf("Info. %s at %d inconsistent\n",m,ut->n1);}
}
//---------------------------------------------------------------------------------
// Message de verification info (Contents)
//---------------------------------------------------------------------------------
void MessageContenu(ptUTIL ut,char* m,int ph)
{
  int i;
  if (ut->rais!=raison) 
    raison = ut->rais;printf("Source: %s\n",raison);
  
  printf("Info %s on [%d..%d]:",m,ut->n1,ut->n2); 
  for (i=ut->n1-1;i<=ut->n2-1;i++)
    if (Choice[i]==ph)
      printf("confirmed at %d\n",i+1);
    else 
      printf("infirmed at %d\n",i+1);
}
//---------------------------------------------------------------------------------
// Analyse des infos et incoherences
//---------------------------------------------------------------------------------
void Retour(ptUTIL r)
{
  int i;
  ptUTIL ut=r;
  raison = "#";
  while (ut!=NULL)
    {
      if (ut->check==false) {ut=ut->suiv;continue;};
      if (ut->tab== 13) 
	{ 
	  i = ut->n1-1;
	  MessageSignal(ut,((Choice[i-1]==12)&&(Choice[i]>=0)&&(Choice[i]<=2)),"start forward");
	 }
      else if (ut->tab== 14)
	{
	  i = ut->n1-1;
	  MessageSignal(ut,((Choice[i]==12)&&(Choice[i]>=3)&&(Choice[i]<=5)),"start reverse");
	  }
      else if (ut->tab== 15)
	{
	  i = ut->n1-1;
	  MessageSignal(ut,((Choice[i-1]>=0)&&(Choice[i-1]<=2)&&(Choice[i]==12)),"stop forward");
	}
      else if (ut->tab==16)	 
	{
	  i = ut->n1-1;
	  MessageSignal(ut,((Choice[i]==12)&&(Choice[i+1]>=3)&&(Choice[i+1]<=5)),"Stop reverse");
	}
      else if (ut->tab==17)
	{
	  i = ut->n1-1;
	  MessageSignal(ut,(!((Choice[i-1]>=0)&&(Choice[i-1]<=2))&&(Choice[i-1]!=12)&&(Choice[i]>=0)&&(Choice[i]<=2)),"acceptor forward");
	}
      else if (ut->tab== 18)
	{
	  i = ut->n1-1;
	  MessageSignal(ut,(!((Choice[i+1]>=3)&&(Choice[i+1]<=5))&&(Choice[i+1]!=12)&&(Choice[i]>=3)&&(Choice[i]<=5)),"acceptor reverse");
	}
      else if (ut->tab== 19)
	{
	  i = ut->n1-1;
	  MessageSignal(ut,((Choice[i-1]>=0)&&(Choice[i-1]<=2)&&!((Choice[i]>=0)&&(Choice[i]<=2))&&(Choice[i]!=12)),"donnor forward");
	}
      else if (ut->tab== 20)	 
	{
	  i = ut->n1-1;
	  MessageSignal(ut,(!((Choice[i]>=3)&&(Choice[i]<=5))&&(Choice[i]!=12)&&(Choice[i+1]>=3)&&(Choice[i+1]<=5)),"donnor reverse");
	} 
      else if (ut->tab == 0)	  
	  MessageContenu(ut,"exonic in frame 1",0);
      else if (ut->tab == 1)
	  MessageContenu(ut,"exonic in frame 2",1);
      else if (ut->tab == 2)
	  MessageContenu(ut,"exonic in frame 3",2);
      else if (ut->tab == 3)
	  MessageContenu(ut,"exonic in frame -1",3);
      else if (ut->tab == 4)
	  MessageContenu(ut,"exonic in frame -2",4);
      else if (ut->tab == 5)
	  MessageContenu(ut,"exonic in frame -3",5);
      else if (ut->tab == 6)
	  MessageContenu(ut,"intronic forward",6);
      else if (ut->tab == 7)
	  MessageContenu(ut,"intronic reverse",7);
      else if (ut->tab== 8)
	  MessageContenu(ut,"intergenic",8);
      else if (ut->tab == 9)
	  MessageContenu(ut,"UTR5 forward",9);
      else if (ut->tab == 10)
	  MessageContenu(ut,"UTR5 reverse",10);
      else if (ut->tab == 11)
	  MessageContenu(ut,"UTR3 forward",11);
      else if (ut->tab == 12)
	  MessageContenu(ut,"UTR3 reverse",12);
      ut=ut->suiv;
    }
}


%}


%union {
  int entier;
  double reel;
  int tab;
  bool check;
  char* chaine;
};


/*---------- mots cles ----------*/
 
%token <chaine> RAISON
%token <entier> NUM BRIN
%token <reel> POIDS
%token <check> CH
%token CO CF PT NL
%token <tab> TYPE EXON INTRON INTERGENIQUE UTR3 UTR5

/*------ ordre de priorite ------*/

%start input


/*---------------------------------------------------------------------------*/
/*---------------------------  REGLES DE LA GRAMMAIRE -----------------------*/
/*---------------------------------------------------------------------------*/
%%

/*--- input : une entree quelconque ---*/
input   :
/*empty */ 
{
  //cellule INIT
  raison = "unknown";
  SignalUser = mall();
  SignalUser->tab = -1;
  SignalUser->n1 = 0;
  SignalUser->n2 = 0;
  SignalUser->delta = 0;
  SignalUser->rais = raison;
  SignalUser->check = false;
  SignalUser->suiv = NULL;

  ContentsUser = mall();
  ContentsUser->tab = -1;
  ContentsUser->n1 = 0;
  ContentsUser->n2 = 0;
  ContentsUser->delta = 0;
  ContentsUser->rais = raison;
  ContentsUser->check = false;
  ContentsUser->suiv = NULL;
}
| input NL
| input line NL
;

line :

RAISON 
{
  raison = $1;
}

|TYPE BRIN NUM POIDS 
{ 
  ptUTIL us = mall(); 
  us->tab = $1+$2;
  us->n1 = $3;
  us->n2 = $3;
  us->delta = $4;
  us->rais = raison;
  us->suiv = suivant(us,SignalUser);
}

|EXON BRIN NUM CO NUM PT NUM CF POIDS 
{ 
  ptUTIL us = mall(); 
  us->tab = 3*$2+$3-1;
  us->n1 = $5;
  us->n2 = $7;
  us->delta = $9;
  us->rais = raison;
  us->suiv = suivant(us,ContentsUser);
}

|INTRON BRIN CO NUM PT NUM CF POIDS 
{ 
  ptUTIL us = mall(); 
  us->tab = 6+$2;
  us->n1 = $4;
  us->n2 = $6;
  us->delta = $8;
  us->rais = raison; 
  us->suiv = suivant(us,ContentsUser);
}

|UTR3 BRIN CO NUM PT NUM CF POIDS 
{
  ptUTIL us = mall(); 
  us->tab = 11+$2;
  us->n1 = $4;
  us->n2 = $6;
  us->delta = $8;
  us->rais = raison;
  us->suiv = suivant(us,ContentsUser);
}

|UTR5 BRIN CO NUM PT NUM CF POIDS 
{ 
  ptUTIL us = mall(); 
  us->tab = 9+$2;
  us->n1 = $4;
  us->n2 = $6;
  us->delta = $8;
  us->rais = raison;
  us->suiv = suivant(us,ContentsUser);
}

|INTERGENIQUE CO NUM PT NUM CF POIDS 
{ 
  ptUTIL us = mall(); 
  us->tab = 8;
  us->n1 = $3;
  us->n2 = $5;
  us->delta = $7;
  us->rais = raison;
  us->suiv = suivant(us,ContentsUser);
}

; 


%%

/*------------------------------------------------------------------------*/
/*---------------------- DECLARATION DE PROCEDURES C ---------------------*/
/*------------------------------------------------------------------------*/


int Utilisateur(char *nom_fich, ptUTIL *SigInfo, ptUTIL *ConInfo) {

  yyin = fopen(nom_fich, "r" );
  if (yyin) {
    yyparse();
    fclose(yyin);
    *SigInfo = SignalUser;
    SignalUser = NULL;
    *ConInfo = ContentsUser;
    ContentsUser = NULL;
    return 0;
  } else  return (1);
}

