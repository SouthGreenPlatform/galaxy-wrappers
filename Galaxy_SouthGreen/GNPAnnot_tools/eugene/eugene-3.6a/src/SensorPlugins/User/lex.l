/* ------------------------------------------------------------------
   Copyright (C) 2004 INRA <eugene@ossau.toulouse.inra.fr>
  
   This program is open source; you can redistribute it and/or modify
   it under the terms of the Artistic License (see LICENSE file).
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
  
   You should have received a copy of Artistic License along with
   this program; if not, please see http://www.opensource.org
  
   $Id: lex.l,v 1.3 2004-09-16 13:36:14 cros Exp $
   ------------------------------------------------------------------
   File:     lex.l
   Contents: ANALYSEUR LEXICAL du langage utilisateur eugene           
     repere les elements de la grammaire :       
       - RAISON START STOP DONNEUR ACCEPTEUR     
       - EXON INTRON INTERGENIQUE NUM CERTITUDE  
   ------------------------------------------------------------------*/


%{
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "yacc.h" 
%}

WS [ \t]+
NL [\n]+  
RC [\r]+

%%

{WS} 
{NL} {return(NL);}
{RC}   


"#".*    {yylval.chaine = strdup(yytext);return(RAISON);}

"start"  {yylval.tab = 13;return(TYPE);}

"stop" {yylval.tab = 15; return(TYPE);}

"donor"  {yylval.tab = 19;return(TYPE);}

"acceptor"  {yylval.tab = 17;return(TYPE);}

"exon"  {yylval.tab = 0;return(EXON);}

"intron"  {yylval.tab = 0;return(INTRON);}

"utr3" {yylval.tab = 0;return(UTR3);}

"utr5" {yylval.tab = 0;return(UTR5);}

"intergenic"  {yylval.tab = 0;return(INTERGENIQUE);}

"f"  { yylval.entier = 0 ; return(BRIN);}
"r"  { yylval.entier = 1; return(BRIN);}

-?[0-9]+  {yylval.entier = atoi( yytext ); return(NUM);}

"[" {return(CO);}
"]" {return(CF);}
"."+ {return(PT);}

"check"  { yylval.check = true; return(CH);}
"nocheck"  { yylval.check = false; return(CH);}

"infinity" {yylval.reel = -log(0.0); return(POIDS);}

"-infinity" {yylval.reel = log(0.0); return(POIDS);}


%-?[0-9]|[0-9]+"."[0-9]+|-?0x[0-9a-h]"."?[0-9a-h]*p[+-][0-9]+  {yylval.reel = strtod(yytext,NULL); return(POIDS);}


-?[0-9]+"."[0-9]+  {yylval.reel = strtod(yytext,NULL); return(POIDS);}

.  printf( "Unknown character : %s\n", yytext );


%%

// necessary when compiling with g++
int yywrap(void) {return 1;} 
