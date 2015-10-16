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
// $Id: markov.cc,v 1.20 2008-07-02 11:19:44 tschiex Exp $
// ------------------------------------------------------------------
// File:     markov.cc
// Contents: class Chaine, template TabChaine
// ------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <typeinfo>
#include <ctype.h>
#include <assert.h>
#include <string>

#include "EndianConv.h"

#include "markov.h"


const char* CODEGENETIQUE = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
const int NBRECODONS[64]     = {2,2,2,2,4,4,4,4,6,6,6,6,3,3,1,3,2,2,2,2,4,4,4,4,6,6,6,6,6,6,6,6,2,2,2,2,4,4,4,4,4,4,4,4,4,4,4,4,3,2,3,2,6,6,6,6,3,2,1,2,6,2,6,2};

/////////////////////////////////////////////////////////
//////                CLASSE Chaine                //////
/////////////////////////////////////////////////////////
//////                  DEBUT                      //////
/////////////////////////////////////////////////////////

// constructeur par defaut
Chaine :: Chaine()
  : taille (0), lettre(NULL) {
  }

// destructeur par defaut
Chaine :: ~Chaine() {
  if (lettre != NULL) delete [] lettre; 
}

// constructeur
Chaine :: Chaine(const char* const data)
  : taille(strlen(data)),
    lettre(new unsigned char[taille+1])
{
  for (int i=0;i< (signed)taille;i++) {
    lettre[i]=data[i];
  }
  lettre[taille]=0;
}

// Constructeurs/Destructeurs de chaines nucleiques et proteiques
ChaineADN :: ChaineADN(void) : Chaine("ACGT") {}
ChaineADN :: ~ChaineADN(void) {}
ChainePROT :: ChainePROT(void) : Chaine("ACDEFGHIKLMNPQRSTVWY") {}
ChainePROT :: ~ChainePROT (void) {}
ChainePROT21 :: ChainePROT21(void) : Chaine("ACDEFGHIKLMNPQRSTVWY*") {}
ChainePROT21 :: ~ChainePROT21 (void) {}
ChainePROT24 :: ChainePROT24(void) : Chaine("ABCDEFGHIKLMNPQRSTVWXYZ*") {}
ChainePROT24 :: ~ChainePROT24 (void) {}

// Pour acceder aux valeurs du tableau
char Chaine :: operator [] (int i) const
{
  return lettre[i];
}

int Chaine :: operator [] (char c) const
{
  int i;
  for (i=0;i<(signed)taille;i++) {
    if(lettre[i]==c) break;
  }
  //  if(i!=taille)return i; // version securisee
  return i; // permet si ==taille de tester si la lettre existe
}

// Visualisation du tableau
void Chaine :: affichage() const {
  printf ("affichage de l'alphabet: ");
  printf ("taille=%d, lettres=%s\n",taille,lettre);
  //  for (int i=0;i<taille;i++) {
  //    printf("lettre[%d]:%c, ",i,(*this)[i]);
  //  }
  //  printf("\n");
}


// Indices des mots d'une certaine longueur
char* Chaine :: code2mot(int code, int lgr) const
{
  if (lgr==0) return NULL;
  int* clef;
  char* tampon;
  int i,reste;
  tampon = new char[lgr+1];
  for (i=0;i<lgr;i++) {
    tampon[i]='X';
  }
  clef= new int[lgr];
  reste= code;
  for (i=lgr-1; i>0; i--) {
    clef[i]=(int)pow((double)taille,(double)i);
    tampon[lgr-i-1]=lettre[(int)reste/clef[i]];
    reste = reste % clef[i];
  }
  tampon[lgr-1]=lettre[reste];
  tampon[lgr]=0;
  delete [] clef;
  //printf("fonction code2mot:\n lgr=%d,lgr du tampon=%d tampon=%s\n",
  //lgr,strlen(tampon),tampon);
  return tampon;
}

int Chaine :: mot2code(const char* const mot) const
{
  int i,code,lgr;
  lgr=strlen(mot);
  code=0;
  for(i=0;i<lgr;i++) {
    code+= (*this)[mot[i]]*(int)pow(taille,(lgr-1-i));
  }
  return code;
}

int Chaine :: mot2code(const char* const mot,int lgr) const
{
  int i,j,lgrmot,code;
  lgrmot=strlen(mot); // vraie longueur
  code=0;
  for(i=0;i<lgr;i++) {
    j= ( (i<lgrmot) ? (*this)[mot[i]] : 0 );
    code+= j *(int)pow(taille,(lgr-1-i));
  }
  return code;
}
int Chaine::mot2code(const char* const mot,int lgr, int deb) const
{
  int i,code;
  code=0;
  for(i=deb ; i< deb+lgr ; i++) {
    code+= (*this)[mot[i]]*(int)pow(taille,(lgr-1-i+deb));
  }
  return code;
}
// code du mot de longueur lgr commencant par seq[deb-fin],
// en le completant par la premiere lettre de l'alphabet
int Chaine::mot2code(const char* const mot,int lgr, int deb, int fin) const
{
  int i,j,code;
  code=0;
  for(i=deb ; i< deb+lgr ; i++) {
    j= ( (i<=fin) ? (*this)[mot[i]] : 0 );
    code+= j *(int)pow(taille,(lgr-1-i+deb));
  }
  return code;
}

/////////////////////////////////////////////////////////
//////                CLASSE Chaine                //////
/////////////////////////////////////////////////////////
//////                    FIN                      //////
/////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////
//////                CLASSE TabChaine             //////
/////////////////////////////////////////////////////////
//////                    DEBUT                    //////
/////////////////////////////////////////////////////////

// constructeurs/destructeur
template<class CHAINE, typename T> TabChaine<CHAINE,T> :: TabChaine() 
  : lgrmax(0), alphabet(new ChaineADN)
{
  indexlgrmots=NULL;
  nbrevaleurs=0;
  VAL=NULL;
}
template<class CHAINE, typename T> TabChaine<CHAINE,T> :: ~TabChaine()
{
  if(indexlgrmots!=NULL) delete [] indexlgrmots;
  if(VAL!=NULL) delete [] VAL;
}
template<class CHAINE, typename T> TabChaine<CHAINE,T> :: TabChaine(int ordre, CHAINE* _chaine)
  : lgrmax(ordre+1), alphabet(_chaine)
{
  int i;
  indexlgrmots= new int[lgrmax+1];
  for (i=0; i< lgrmax+1;i++) {
    indexlgrmots[i]= int (pow (alphabet->taille,i) - 1) / (alphabet->taille - 1);
  }
  nbrevaleurs= int (pow (alphabet->taille, lgrmax+1 ) - 1) / (alphabet->taille - 1);
  VAL= new T[nbrevaleurs];
  for (i=0; i<nbrevaleurs ;i++) { VAL[i] = 0; }
  /*
    if (typeid(VAL[0]) == typeid(i)) {
    for (i=0; i<nbrevaleurs ;i++) { VAL[i] = 0; }
    }
    else{
    VAL[0]=1.0;
    for (i=1; i<nbrevaleurs ;i++) VAL[i] = 0.0;
    }
  */
}

// POURQUOI NE COMPILE PAS AVEC "CONST"?!!!
template<class CHAINE, typename T> void TabChaine<CHAINE,T> :: affichage(int l)
{
  printf("Affichage de la classe template TabChaine:\n");
  if(alphabet != NULL) alphabet->affichage();
  printf("lgrmax=%d, taille alphabet=%d, nbrevaleurs=%d\n",lgrmax,alphabet->taille,nbrevaleurs);
  if(l>=0) affichagevaleurs(l);
}

template<class CHAINE, typename T> void TabChaine<CHAINE,T> :: affichagevaleurs(int l) 
{
  int i;
  std::string affiche;
  //  char* affiche = new char[20];

  printf("affichage des valeurs:\n");  
  if ( ( typeid(T) == typeid(int) ) || ( typeid(T) == typeid(unsigned)) ) affiche= "%3d %s : %3d\n";
  else affiche= "%3d %s : %5.4lf\n";

  for(i=1; i<nbrevaleurs;i++) {
    if ((l==0) || (indice2lgrmot(i) == l)) {
      if (typeid(T) == typeid(unsigned short int)) printf(affiche.c_str(),i,indice2mot(i),usi2real(VAL[i]));
      else printf(affiche.c_str(),i,indice2mot(i),VAL[i]);
    }
  }
}

template<class CHAINE, typename T> int TabChaine<CHAINE,T> :: indice2index(int indice) const
{
  int i,j = 0;
  if(indice<nbrevaleurs) {
    for (i=0;i<=lgrmax;i++) {
      if (indice>=indexlgrmots[i]) j=indexlgrmots[i];
      else break;
    }
    return j;
  }
  else return 0;
}

template<class CHAINE, typename T> int TabChaine<CHAINE,T> :: indice2lgrmot(int indice) const
{
  int i,j = 0;
  if(indice<nbrevaleurs) {
    for (i=0;i<=lgrmax;i++) {
      if (indice>=indexlgrmots[i]) j=i;
      else break;
    }
    return j;
  }
  else return 0;
}

template<class CHAINE, typename T> T TabChaine<CHAINE,T> :: indice2VAL(int indice) const
{ return VAL[indice]; }

template<class CHAINE, typename T> int TabChaine<CHAINE,T> :: mot2indice(char* mot) const
{
  int lgr=strlen(mot);
  if(lgr<=lgrmax) return(indexlgrmots[lgr]+alphabet->mot2code(mot));
  else return 0;
}
template<class CHAINE, typename T> int TabChaine<CHAINE,T> :: mot2indice(char* mot,int lgr) const
{
  if(lgr<=lgrmax) return(indexlgrmots[lgr]+alphabet->mot2code(mot,lgr));
  else return 0;
}
template<class CHAINE, typename T> int TabChaine<CHAINE,T> :: mot2indice(char* mot,int lgr, int deb) const
{
  if(lgr<=lgrmax) return(indexlgrmots[lgr]+alphabet->mot2code(mot,lgr,deb));
  else return 0;
}

template<class CHAINE, typename T> char* TabChaine<CHAINE,T> :: indice2mot(int indice) const
{
  int lgr= indice2lgrmot(indice);
  int index= indice2index(indice);
  //  printf("indice=%d,lgr=%d,lgr=%d,index=%d\n",indice,lgr,
  //  strlen(alphabet->code2mot((indice-index),lgr)),index);
  //  char *mot;
  //  mot= alphabet->code2mot((indice-index),lgr);
  //  return mot;
  return (alphabet->code2mot((indice-index),lgr));
}

template<class CHAINE, typename T> void TabChaine<CHAINE,T> :: incremente(int indice, int n)
{
  if (indice < nbrevaleurs) VAL[indice] += n;
  if (indice2lgrmot(indice)==1) VAL[0]  += n;
}

template<class CHAINE, typename T> void TabChaine<CHAINE,T> :: initialisation()
{
  for (int i=0; i<nbrevaleurs;i++) {
    VAL[i]=0;
  }
}

template<class CHAINE, typename T> void TabChaine<CHAINE,T> :: initialisation(double GCrate)
{
  T precval;
  //  char tampon;  // derniere lettre du mot
  int der=0;
  VAL[0]=1;
  for (int i=1; i<nbrevaleurs ; i++) {
    precval= VAL[indiceprefixe(i)];
    der = indiceprefixe(indice2mot(i),1,indice2lgrmot(i)-1,indice2lgrmot(i))-1;
    //    tampon = alphabet->lettre[indiceprefixe(indice2mot(i),1,
    //    indice2lgrmot(i)-1,indice2lgrmot(i)) - 1];
    VAL[i] = precval * ( ((der==1) || (der==2)) ? GCrate/2 : (1-GCrate)/2 ) ;
    //    VAL[i] = precval * ( ((strcmp(&tampon,"G")) &&
    //    (strcmp(&tampon,"C"))) ? GCrate : 1-GCrate ) ;
    //    printf("i=%d,mot=%s,tampon=%c\n",i,indice2mot(i),
    //    alphabet->lettre[der]);
  }
}

// permet d'eviter des probas de 0 (et du coup des log de 0...)
template<class CHAINE, typename T> void TabChaine<CHAINE,T> :: pseudocount (int pseudocomptes)
{
  for (int i=0;i<nbrevaleurs;i++)
    incremente(i,pseudocomptes);
}

template<class CHAINE, typename T> void TabChaine<CHAINE,T> :: seq2compte(char* seq, int parcodon)
{
  int i,j;
  if ( (seq!=NULL)&&(strlen(seq)>lgrmax)) {
    for(i=0; i<lgrmax-1 ;i++) { // begining of the seq (before the order)
      for(j=0;j<=i;j++) incremente(mot2indice(seq,j+1,i-j));
    }
    
    for(i=lgrmax-1 ; i<strlen(seq) ; i++) {
      for(j=0;j<lgrmax;j++) {
	incremente(mot2indice(seq,(lgrmax-j),i-(lgrmax-j-1)));
      }
    }
  }
}

template<class CHAINE, typename T> void TabChaine<CHAINE,T> :: seq2compte(char* seq, int debut, int fin, int parcodon)
{
  int i,j;
  if ((seq!=NULL)&&(strlen(seq)>lgrmax)) {
    for (i=debut; ((i<lgrmax-1)&&(i<=fin))  ; i++) {  
      // the amount context is not sufficient for the order
      // this case should not happen
      for(j=0;j<=i;j++) incremente(mot2indice(seq,j+1,i-j));
    }
    for(i= ((debut >= lgrmax-1)? debut : lgrmax-1); i<= fin ; i++) {
      for(j=0;j<lgrmax;j++) {
	incremente(mot2indice(seq,(lgrmax-j),i-(lgrmax-j-1)));
	//	fprintf(stdout,"incremente2 %d\n",mot2indice(seq,(lgrmax-j),i-(lgrmax-j-1)));
      }
    }
  }
}

// lit une sequence format fasta, retourne 0 si pb, 1 si OK, et si un
// caractere n'est pas prevu (hors alphabet), retourne 2 et Sequence=NULL
template<class CHAINE, typename T> int TabChaine<CHAINE,T> :: fichier2seq (FILE *fp, char* &Sequence)
{
  //  fprintf(stderr,"fichier2seq: entree,seq=%s\n",Sequence);
  const unsigned int  INCR_SIZE = 10000;
  const unsigned int  INIT_SIZE = 10000;
  char  *seqname, Line [1000];
  int  Len,Size;
  char  Ch;
  
  while  ((Ch = fgetc (fp)) != EOF && Ch != '>');
  // None found
  if  (Ch == EOF)
    return 0;
  
  // Skip the header
  fgets (Line, 1000 , fp);
  Len = strlen(Line);
  assert (Line [Len - 1] == '\n');
  seqname = strtok (Line, " \t\n");
  
  // Read the sequence
  Sequence = (char *) malloc (INIT_SIZE * sizeof(char));
  Size = INIT_SIZE;
  Len = 0;
  
  while  ((Ch = fgetc (fp)) != EOF && Ch != '>') {
    if  (isspace (Ch))
      continue;
    
    if  (Len >= Size)
      {
	Size += INCR_SIZE;
	Sequence = (char *) realloc (Sequence, sizeof(char)*Size);
      }
    Ch = toupper (Ch);
    
    if ( strcmp((const char *)alphabet->lettre,"ACGT")==0) {
      switch  (Ch)
	{
	case  'A' : case  'C' : case  'G' : case  'T' :
	  break;
	case  'S' : Ch = 'C';
	  break;
	case  'W' : Ch = 'A';
	  break;
	case  'R' : Ch = 'A';
	  break;
	case  'Y' : Ch = 'C';
	  break;
	case  'M' : Ch = 'A';
	  break;
	case  'K' : Ch = 'G';
	  break;
	case  'B' : Ch = 'C';
	  break;
	case  'D' : Ch = 'A';
	  break;
	case  'H' : Ch = 'A';
	  break;
	case  'V' : Ch = 'A';
	  break;
	case  'N' : Ch = 'A';
	  break;
	default :
	  fprintf (stderr, "Unexpected character `%c\' fonction fichier2seq ((multi)fasta nucleique requis)-> sequence %s rejetee\n", Ch,seqname);
	  Sequence=NULL;
	  return 2;
	}
    }
    else{
      if ( strcmp((const char *)alphabet->lettre,"ACDEFGHIKLMNPQRSTVWY*")==0) {
	if ( alphabet->operator[](Ch) == alphabet->taille) {
	  fprintf (stderr, "Unexpected character `%c\' fonction fichier2seq ((multi)fasta proteique requis)-> sequence %s rejetee\n", Ch,seqname);
	  Sequence=NULL;
	  return 2;
	}
      }
      else{
	if (alphabet->operator[](Ch) == alphabet->taille) {
	  fprintf (stderr, "Unexpected character `%c\' fonction fichier2seq-> sequence %s rejetee\n", Ch,seqname);
	  Sequence=NULL;
	  return 2;
	}
      }
    }
    Sequence [Len++] = Ch;
  }
  
  // TMP!! Cas particulier: sequence proteique a terminer par un stop
  // (sinon P(TGA)=0)!!
  // A traiter peut etre plus proprement...
  if ( strcmp((const char *)alphabet->lettre,"ACDEFGHIKLMNPQRSTVWY*")==0) {
    if  (Len >= Size)
      {
	Size += 1;
	Sequence = (char *) realloc (Sequence, sizeof(char)*Size);
      }
    Sequence [Len++] = '*';
  }
  
  // prepare for next call
  if (Ch != EOF) 
    ungetc(Ch,fp);
  
  Sequence[Len] = '\0';
  
  //   fprintf(stderr,"fichier2seq: sortie,seq=%s\n",Sequence);
  return  1;
}

template<class CHAINE, typename T> void TabChaine<CHAINE,T> :: fichier2compte (FILE *fp, int parcodon)
{
  char* Sequence=NULL;
  while (fichier2seq(fp,Sequence)) {
    seq2compte(Sequence,parcodon);
    free(Sequence);
    Sequence=NULL;
  };
}

template<class CHAINE, typename T> void TabChaine<CHAINE,T> :: fichier2compte (FILE *fp, int debut, int fin, int parcodon)
{
  char* Sequence=NULL;
  while (fichier2seq(fp,Sequence)) {
    seq2compte(Sequence,debut,fin,parcodon);
    free(Sequence);
    Sequence=NULL;
  };
}

template<class CHAINE, typename T> void TabChaine<CHAINE,T> :: sauve2fichier (FILE *fp)
{
  int    i;
  unsigned short  buffer;
  
  (void) fwrite(&lgrmax, sizeof(int), 1, fp);
  (void) fwrite(&alphabet->taille, sizeof(int), 1, fp);
  (void) fwrite(&nbrevaleurs, sizeof(int), 1, fp);
  
  for (i = 0;  i < nbrevaleurs;  i ++) {
    //buffer = (unsigned short) floor((double) (VAL[i]/BINARY_PRECISION+0.5));
    //(void) fwrite(&buffer,  sizeof(unsigned short), 1, fp);
    (void) fwrite(&VAL[i],  sizeof(T), 1, fp);
  }
  return;
}


template<class CHAINE, typename T> int TabChaine<CHAINE,T> :: chargefichier (FILE *fp)
{
  int  A, N;
  unsigned int M;
  char endian;
  int Ok;
  
  endian = 0;
  M = A = N = 0;
  
  Ok = fread (&M, sizeof(int), 1, fp);
  if (!Ok) return 1;
  
  if (M == LEndianReverse(lgrmax))
    {
      endian = 1;
      M = lgrmax;
    }
  
  Ok = fread (&A, sizeof(int), 1, fp);
  if (!Ok) return 1;
  
  if (endian) A = LEndianReverse(A);
  
  Ok = fread (&N, sizeof(int), 1, fp);
  if (!Ok) return 1;
  
  if (endian) N = LEndianReverse(N); 
  
  if ( ((signed)M != lgrmax) || (A !=(signed)alphabet->taille) || ((signed)N != nbrevaleurs) ) fprintf(stderr,"markov.cc : Incompatibility between model expected and read in function chargefichier: M=%d, lgrmax=%d, A=%d,alphabet->taille =%d, N=%d, nbrevaleurs=%d\n",M,lgrmax,A,alphabet->taille,N,nbrevaleurs);
  assert ((signed)M == lgrmax && A == (signed)alphabet->taille && (signed)N == nbrevaleurs);
  
  // compatibilite des matrices entre machines garantie uniquement
  // pour le type usi
  if (typeid(T) == typeid(unsigned short int)) {
    Ok = fread (VAL, sizeof(unsigned short), nbrevaleurs, fp);
    if (Ok != nbrevaleurs) return 1;
    
    if (endian)
      for (int i = 0;  i < nbrevaleurs;  i ++)
	VAL[i] = SEndianReverse(VAL[i]);
  }
  
  else{
    Ok = fread (VAL, sizeof(T), nbrevaleurs, fp);
    if (Ok != nbrevaleurs) return 1;
  }
  return 0;
}

template<class CHAINE, typename T> int TabChaine<CHAINE,T> :: indiceprefixe (char* mot) const
{
  int lgr;
  lgr=strlen(mot);
  return (indexlgrmots[lgr]+alphabet->mot2code(mot,lgr,0,(lgr-2)));
}

// donne l'indice du premier "declinant" ou mot de longueur lgr
// contenant le prefixe seq[deb-fin]. Exemple:
// avec (seq ATGGC, lgr 3, deb 1, fin 2) trouve le premier mot de 3 lettres 
// de prefixe TG (="declinant" de TGN), soit l'indice de TGA!
template<class CHAINE, typename T> int TabChaine<CHAINE,T> :: indiceprefixe (char* seq, int lgr, int deb, int fin) const
{
  return (indexlgrmots[lgr]+alphabet->mot2code(seq, lgr, deb, fin));
}

// Plus simple, avec l'indice, renvoie l'indice du mot de taille l-1
template<class CHAINE, typename T> int TabChaine<CHAINE,T> :: indiceprefixe (int indice) const
{
  return (int) ((indice-1)/alphabet->taille);
}

template<class CHAINE, typename T> int TabChaine<CHAINE,T> :: usi() const
{
  if (typeid(T) == typeid(unsigned short int)) return 1;
  return 0;
}

// necessite en argument une classe de meme type chaine et de meme taille
// pour passer des frequences aux probabilites
template<class CHAINE, typename T> void TabChaine<CHAINE,T> :: compte2probas (const TabChaine<CHAINE, int>* comptage, int occurence_min)
{
  int i,j,cumul,indicedeclinant;
  char* mot= new char[lgrmax];
  if (nbrevaleurs == comptage->nbrevaleurs) {
    if ((typeid(T) == typeid(unsigned short int)) || (typeid(T) == typeid(int)) ) VAL[0] = 1;
    else VAL[0] = 1.0;
    for(i=1;i<nbrevaleurs;i++) {
      cumul=0;
      mot= indice2mot(i);
      indicedeclinant=indiceprefixe(indice2mot(i),indice2lgrmot(i),0,indice2lgrmot(i)-2);
      for(j=indicedeclinant; j< indicedeclinant+alphabet->taille; j++) {
	cumul+= comptage->VAL[j];
	// Ex. avec adn: Si mot(i)=CG, cumule N(CA),N(CC),N(CG) et N(CT)
      }
      //      if (indice2lgrmot(i) == lgrmax) {
      //	fprintf(stdout,"Proba de %s = ( N(%s)=%d / Somme de n(%s) a n(%s)= %d )= %f\n",indice2mot(i),indice2mot(i),comptage->VAL[i],indice2mot(indicedeclinant),indice2mot(indicedeclinant+alphabet->taille-1),cumul,(double)comptage->VAL[i] / (double)cumul);
      //	fflush(stdout);
      //    }
      if (cumul < occurence_min) { // pas assez d'info pour calculer une proba
	if (indice2lgrmot(i)==1) { // et en plus on n'est qu'a l'ordre 0
	  fprintf(stderr,"WARNING!! Not enough information to compute P(%s) AT ORDER ZERO!\n",mot);
	}
	else {
	  VAL[i]= VAL[mot2indice(mot,indice2lgrmot(i)-1,1)];
	  //fprintf(stderr,"Warning : not enough information to compute P(%s), taking P(%s)\n",mot,indice2mot(mot2indice(mot,indice2lgrmot(i)-1,1)));
	}
      }
      else {
	if ( (typeid(T) == typeid(unsigned short int)) || (typeid(T) == typeid(int)) )
	  VAL[i]= ((cumul==0) ? 0 : real2usi( (double)comptage->VAL[i] / (double)cumul) );
	else
	  VAL[i]=((cumul==0) ? 0 : (double)comptage->VAL[i] / (double)cumul ) ;
      }
    }
  }
  else{
    fprintf(stderr,"\n\n\nERROR !!\n\n\nError in markov.cc (method compte2probas)\n\n\nERROR !!\n\n\n");
  }
  delete [] mot;
}

// Convertisseurs de types double/USI (double/unsigned short int)
template<class CHAINE, typename T>  unsigned short int TabChaine<CHAINE,T> :: real2usi (double x) {
  return ( (unsigned short int) floor ( (double)(x*65535.0+0.5) ) );
}

template<class CHAINE, typename T> double TabChaine<CHAINE,T> :: usi2real (unsigned short int n) {
  return ( (double)n / 65535.0 );
}

template<class CHAINE, typename T> int TabChaine<CHAINE,T> :: nombredegc (int indice) const
{
  int compt=0;
  char* mot= new char[lgrmax];
  mot= indice2mot(indice);
  for(int i=0;i<lgrmax;i++) {
    if ((mot[i]=='G')||(mot[i]=='C')) compt++;
  }
  return compt;
}

///////////////////////////////////////////////////////
//   D'UNE SEQ AUX PROBABILITES ET LOG de PROBAS  /////
///////////////////////////////////////////////////////

// proba d'apparition de la lettre en position i selon l'ordre maximum
template<class CHAINE, typename T> T TabChaine<CHAINE,T> :: proba (char *seq, int i) const
{
  //  printf("chaine %s de lgr %d; Proba(%c) en position %d: %lf\n",seq,lgrseq,seq[i],i, ((i > lgrmax-1)? (VAL[mot2indice(seq,lgrmax,i-lgrmax+1)]):(VAL[mot2indice(seq,i+1,i)])));
  if (i > lgrmax-1) return (VAL[mot2indice(seq,lgrmax,i-lgrmax+1)]);
  else          return (VAL[mot2indice(seq,i+1,0)]);
}

// proba d'apparition de la lettre en position i selon un ordre particulier
template<class CHAINE, typename T> T TabChaine<CHAINE,T> :: proba (char *seq, int i, int ordre) const
{
  if (i > ordre) return (VAL[mot2indice(seq,ordre+1,i-ordre)]);
  else          return (VAL[mot2indice(seq,i+1,0)]);
}

// proba d'apparition de la sequence seq[deb-fin]
template<class CHAINE, typename T> T TabChaine<CHAINE,T> :: probaseq (char *seq, int deb, int fin) const
{
  T prob=1.0;
  for (int i=deb; i<=fin; i++) {
    prob*=proba(seq,i);
  }
  return prob;
}
// proba d'apparition de la sequence seq // ATTENTION!! utilise strlen, seq doit de finir par 0!
template<class CHAINE, typename T> T TabChaine<CHAINE,T> :: probaseq (char *seq) const
{
  T prob=1.0;
  for (int i=0; i< strlen(seq); i++) {
    prob*=proba(seq,i);
  }
  return prob;
}

// proba d'apparition de la sequence seq[deb-fin] avec un ordre particulier
template<class CHAINE, typename T> T TabChaine<CHAINE,T> :: probaseq (char *seq, int deb, int fin, int ordre) const
{
  T prob=1.0;
  for (int i=deb; i<=fin; i++) {
    prob*=proba(seq,i,ordre);
  }
  return prob;
}

///////// LOG de probas ////////
template<class CHAINE, typename T> T TabChaine<CHAINE,T> :: logproba (char *seq, int i) const
{
  if (i > lgrmax-1) return (log( VAL[mot2indice(seq,lgrmax,i-lgrmax+1)]));
  else          return (log(VAL[mot2indice(seq,i+1,0)]));
}
// logproba d'apparition de la lettre en position i selon un ordre particulier
template<class CHAINE, typename T> T TabChaine<CHAINE,T> :: logproba (char *seq, int i, int ordre) const
{
  if (i > ordre) return (log(VAL[mot2indice(seq,ordre+1,i-ordre)]));
  else          return (log(VAL[mot2indice(seq,i+1,0)]));
}

// logproba d'apparition de la sequence seq[deb-fin]
template<class CHAINE, typename T> T TabChaine<CHAINE,T> :: logprobaseq (char *seq, int deb, int fin) const
{
  T score=0;
  for (int i=deb; i<=fin; i++) {
    score+=logproba(seq,i);
  }
  return score;
}
// logproba d'apparition de la sequence seq[deb-fin] avec un ordre particulier
template<class CHAINE, typename T> T TabChaine<CHAINE,T> :: logprobaseq (char *seq, int deb, int fin, int ordre) const
{
  T score=0;
  for (int i=deb; i<=fin; i++) {
    score+=logproba(seq,i,ordre);
  }
  return score;
}

// logproba d'apparition de la sequence seq // ATTENTION!! utilise strlen, seq doit de finir par 0!
template<class CHAINE, typename T> T TabChaine<CHAINE,T> :: logprobaseq (char *seq) const
{
  T score=0;
  for (int i=0; i< strlen(seq); i++) {
    score+=logproba(seq,i);
  }
  return score;
}

/////////////////////////////////////
//   RAPPORT DE VRAISEMBLANCE  /////
/////////////////////////////////////

template<class CHAINE, typename T> double TabChaine<CHAINE,T> :: 
rapdevrai (char *seq, int n, const TabChaine<CHAINE, double> &autremodele, int ordre) const
{
  return ( proba(seq,n,ordre)/autremodele.proba(seq,n,ordre) );
}

template<class CHAINE, typename T> double TabChaine<CHAINE,T> :: 
lograpdevrai (char *seq, int n, const TabChaine<CHAINE, double> &autremodele, int ordre) const
{
  return ( logproba(seq,n,ordre) - autremodele.logproba(seq,n,ordre)  );
}

template<class CHAINE, typename T> double TabChaine<CHAINE,T> :: 
lograpdevrai (char *seq, int deb, int fin, const TabChaine<CHAINE, double> &autremodele, int ordre) const
{
  int i;
  double val=0.0;
  for (i=deb; i<=fin; i++) {
    val+= lograpdevrai(seq, i, autremodele, ordre);
  }
  return val;
}

//--------------------------
// Fonction qui pour une TabChaine, un indice de codon(de 0 a 63) et un codegenetique,
// cumule les valeurs des cases des codons synonymes et renvoie le total.
template<class CHAINE, typename T> T TabChaine<CHAINE,T> :: cumuleVAL (int indice) const
{
  char* codegenetique=CODEGENETIQUE;
  T cumul=0;
  for (int i=0 ; i<64 ; i++) {
    if ( codegenetique[i] == codegenetique[indice] )
      cumul += VAL[i+21];
    //    printf("i=%d, aa= %c vs %c, cumul=%lf\n",i,codegenetique[i],codegenetique[indice],cumul);
  }
  return cumul;
}

//--------------------------
// Fonction qui pour une TabChaine, un Acide Amine (code une lettre) et un codegenetique,
// cumule les valeurs des cases des codons synonymes et renvoie le total.
template<class CHAINE, typename T> T TabChaine<CHAINE,T> :: cumuleVAL (char aa) const
{
  const char* codegenetique=CODEGENETIQUE;
  T cumul=0;
  for (int i=0 ; i<64 ; i++) {
    if ( codegenetique[i] == aa )
      cumul += VAL[i+21];
    //    printf("i=%d, aa= %c vs %c, cumul=%lf\n",i,codegenetique[i],codegenetique[indice],cumul);
  }
  return cumul;
}

/////////////////////////////////////////////////////////
//////                CLASSE TabChaine             //////
/////////////////////////////////////////////////////////
//////                      FIN                    //////
/////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////
//////                CLASSE UsageCode             //////
/////////////////////////////////////////////////////////
//////                    DEBUT                    //////
/////////////////////////////////////////////////////////


// constructeur par defaut (pas d'argmt nec.)
UsageCode :: UsageCode() : TabChaine <ChaineADN,int> (2, new ChaineADN) , nbreaa(20), nbrecodons(64), offset(21)
{
  //  codegenetique= "                     KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
  codegenetique= CODEGENETIQUE;
  //  codegenetique= "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
  usagecode = new double[nbrecodons+1];
  for (int i=0;i<nbrecodons;i++) {
    usagecode[i]=0.0;
  }
  initialisation();
}
UsageCode :: ~UsageCode()
{
  delete usagecode;
}
void UsageCode :: affichage() const
{
  int i;
  printf("Affichage de la classe UsageCode:\n");
  if(alphabet != NULL) alphabet->affichage();
  printf("lgrmax=%d, taille alphabet=%d, nbrevaleurs=%d\n",lgrmax,alphabet->taille,nbrevaleurs);
  printf("nbre d'aa:%d, nbre de codons (sans stops):%d\n",nbreaa,nbrecodons);
  for(i=0; i< nbrecodons; i++) {
    printf ("codon:%s occurence:%d aa:%c occurence:%d usage:%f\n", indice2mot(i+offset), VAL[i+offset], codegenetique[i], cumuleaa(i), usagecode[i]);
  }
}

int UsageCode :: cumuleaa(int aa) const
{
  int occurence=0;
  for (int i=0;i<nbrecodons;i++) {
    if (codegenetique[i] == codegenetique[aa]) // cherche tous les codons de l'aa
      occurence+= VAL[i+offset]; // occurence de l'aa += occurence du codon
  }
  return occurence;
}

void UsageCode :: compte2usage ()
{
  for (int i=0;i<nbrecodons;i++) {
    usagecode[i]= ( (cumuleaa(i)==0) ? 0.0 : (double)VAL[i+offset] / (double)cumuleaa(i)) ;
  }
}

/////////////////////////////////////////////////////////
//////     CLASSE ProtMat : MATRICE PROTEIQUE      //////
/////////////////////////////////////////////////////////
//////                    DEBUT                    //////
/////////////////////////////////////////////////////////

// constructeur1 (avec taille explicitee)
ProtMat :: ProtMat(char* aa, int _nbreaa) : TabChaine <Chaine,int> (1, new Chaine(aa)) , nbreaa(_nbreaa), offset(_nbreaa)
{
}
// constructeur2 (plus court, mais utilise strlen)
ProtMat :: ProtMat(char* aa) : TabChaine <Chaine,int> (1, new Chaine(aa)) , nbreaa(strlen(aa)), offset(nbreaa)
{
}

//---------------------------------------------
// FONCTION qui lit un fichier MATRICE PROTEIQUE (type Blosum ou Pam)
// et cree la classe ProtMat correspondante, et renvoie 1 si PB.
// (1er armgt= pointeur vers un fichier ouvert de type Blosum ou Pam
//  2eme argument= pointeur VIDE de type ProtMat)
// (passage par reference du pointeur)
int fichier2protmat(FILE *fp, ProtMat* &MAT)
{
  char Line [MAX_LINE];
  int i, j, n;
  char tampon;
  int nmaxAA=50;
  char* STR= new char[nmaxAA+1]; // MAX d'AA pour la matrice
  tampon = '#';

  // on saute les commentaires
  fscanf (fp, "%c", &tampon); 
  while (tampon=='#') {
    if (fgets (Line, MAX_LINE, fp) == NULL) return(1);
    fscanf (fp, "%c", &tampon); 
  }
  ungetc(tampon,fp);

  // on est a la ligne contenant les Acides Amines (une lettre chacun)
  n=0;
  STR[n]='\0';
  while ( (tampon != '\n') && ( ! feof(fp) ) ) { // comptage et memorisation des AA
    tampon=fgetc(fp);
    if (isspace(tampon)) continue;
    if (n>=nmaxAA) {fprintf(stderr,"error in PROTMAT file, blosum/pam format required (too many AA in first line)\n");return(1);}
    STR[n+1] = '\0';
    STR[n] = tampon;
    n++; 
    //    printf("n= %d acides amines=%s\n",n,STR);
  }
  //  printf("lecture de %d acides amines: %s\n",n,STR);
  
  // CREATION de la matrice
  MAT= new ProtMat(STR);
  MAT->initialisation();
  //  MAT->affichage();

  // remplissage des valeurs
  i=-1; // valeurs scannees
  j=-1;  // valeurs a enregistrer
  while ( ! (feof(fp) ) ) {
    fscanf(fp,"%s",STR);
    i++ ;
    if ( (i % (1+MAT->alphabet->taille)) == 0 ) continue;  // 1ere colonne= lettres des AA
    //    printf("MAT[%d][%d]=%s\n",ligne,col,STR);
    j++;
    if (j>=MAT->nbrevaleurs) {fprintf(stderr,"error in PROTMAT file, blosum/pam format required\n");return(1);}
    MAT->VAL[j+MAT->offset+1] = atoi(STR);
  }
  // MAT->affichage(2);
  return(0);
  
}
