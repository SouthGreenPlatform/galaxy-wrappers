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
// $Id: markov.h,v 1.12 2008-07-02 11:19:44 tschiex Exp $
// ------------------------------------------------------------------
// File:     markov.h
// Contents: class Chaine, template TabChaine
// ------------------------------------------------------------------

#ifndef MARKOV
#define MARKOV

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "../../Const.h"

#ifndef CLASSchaine
#define CLASSchaine
class Chaine
{
 public:
  const unsigned int taille;
  unsigned char* const lettre;
  Chaine();
  Chaine(const char* const data);
  ~Chaine();
  void affichage() const;
  char operator [] (int i) const;
  int operator [] (char c) const;
  char* code2mot(int code, int lgr) const;
  int mot2code(const char* const mot) const;
  int mot2code(const char* const mot, int lgr) const;
  int mot2code(const char* const mot,int lgr, int deb) const;
  int mot2code(const char* const mot,int lgr, int deb, int fin) const;
};

class ChaineADN : virtual public Chaine
{
 public:
  ChaineADN();
  ~ChaineADN();
};

class ChainePROT : virtual public Chaine
{
 public:
  ChainePROT();
  ~ChainePROT();
};

class ChainePROT21 : virtual public Chaine
{
 public:
  ChainePROT21();
  ~ChainePROT21();
};

class ChainePROT24 : virtual public Chaine
{
 public:
  ChainePROT24();
  ~ChainePROT24();
};
#endif

#ifndef CLASStabchaine
#define CLASStabchaine
template<class CHAINE, typename T> class TabChaine
{
 public:
  // createurs/destructeur
  TabChaine();
  ~TabChaine();
  TabChaine(int ordre, CHAINE* _chaine); // ordre etant en fait l'ordre max.
  // donnees
  const int lgrmax;       // longeur du mot le plus long (=ordremax+1)
  const CHAINE* alphabet; // on peut faire des Modeles de Markov prot, nuc...
  int* indexlgrmots;      // tab. dont les indices sont les lgrs des mots et
                          // les valeurs les indices pour le tableau VAL du
                          // 1er mot de la longueur correspondante; ex avec
                          // adn(4 lettres): indexlgrmots[0,1,2,3]=0,1,5,21
                          // de VAL[5] a VAL[21] correspondent les mots de
                          // longeur 2 (21-5 ou 4x4 mots) les valeurs
                          // correspondant aux mots de longueur 3 commencent
                          // a VAL[21]
  int nbrevaleurs;  // taille de VAL
  T* VAL;           // tab. contenant soit les occurences des mots
                    // (template int), soit les probas d'apparition de la
                    // derniere letter du mot connaissant le debut
                    // (template double) methodes
  int indice2index(int indice) const;
  int indice2lgrmot(int indice) const;
  T indice2VAL (int indice) const;
  int mot2indice(char* mot) const;                // use strlen
  int mot2indice(char* mot,int lgr) const;
  int mot2indice(char* mot,int lgr, int deb) const;
  char* indice2mot(int indice) const;
  void affichage(int l=-1);  // par defaut, n'affiche pas les valeurs,
                             // l==0 -> ttes les valeurs,
                             // l!=0 -> valeurs pour les mots de lgr l
  void affichagevaleurs(int l=0);
  void incremente (int indice,int n=1); // attention: pour <int>,
                                        // modifie VAL[indice]
  void initialisation ();    // attention: pour <int>, remet le compteur a 0.
  void initialisation (double GCrate);  // init. sachant un GC%
  void pseudocount (int pseudocomptes=1); // incremente toutes les
                                          // valeurs de pseudocomptes
  int fichier2seq(FILE *fp,char* &Sequence);
  // appel du style: 
  // char* Sequence;
  // "while (fichier2seq(fp,Sequence)){
  // seq2compte(Sequence);free(Sequence) }"
  void seq2compte(char* seq, int parcodon=0);
  void seq2compte(char* seq, int debut, int fin, int parcodon=0);
  //  void seq2compte(char* seq, int debut, int fin, int parcodon=0);
  //  void seq2compte(char* seq, int debut, int fin, int phase,
  //                  int parcodon=0); A faire...
  void fichier2compte (FILE *fp, int parcodon=0);
  void fichier2compte (FILE *fp, int debut, int fin, int parcodon=0);
  // sauvegarde; ATTENTION! uniquement pour <double>
  void sauve2fichier (FILE *fp);
  // lecture d'un fichier contenant un modele sauvegarde
  int chargefichier (FILE *fp);
  int indiceprefixe (char* mot) const;
  int indiceprefixe (char* seq, int lgr, int deb, int fin) const;
  int indiceprefixe (int indice) const;
  void compte2probas (const TabChaine<CHAINE, int>* comptage,
		      int occurence_min=100);
  int nombredegc (int indice) const;
  
  // Convertisseurs de types double (e.g. double) et  usi (unsigned short int):
  unsigned short int real2usi (double x); 
  double usi2real (unsigned short int n);
  int usi() const;
  
  T proba (char *seq, int i, int ordre)               const;
  T proba (char *seq, int i)                          const;
  T probaseq (char *seq, int deb, int fin, int ordre) const;
  T probaseq (char *seq, int deb, int fin)            const;
  T probaseq (char *seq)                              const;//a eviter(strlen)!
  T logproba (char *seq, int i, int ordre)            const;
  T logproba (char *seq, int i)                       const;
  T logprobaseq (char *seq, int deb, int fin, int ordre) const;
  T logprobaseq (char *seq, int deb, int fin)         const;
  // aurait pu compacter les 2 methodes precedentes:
  // T logprobaseq (char *seq, int deb, int fin, int ordre=lgrmax-1)  const;
  T logprobaseq (char *seq)                           const;//a eviter(strlen)!
  double rapdevrai    (char *seq, int n, const TabChaine<CHAINE,
		       double> &autremodele, int ordre) const;
  double lograpdevrai (char *seq, int n, const TabChaine<CHAINE,
		       double> &autremodele, int ordre) const;
  double lograpdevrai (char *seq, int deb, int fin, const TabChaine<CHAINE,
		       double> &autremodele, int ordre) const;
  T cumuleVAL(int indice) const;
  T cumuleVAL(char aa) const;
};
#endif

#ifndef CLASSusagecode
#define CLASSusagecode
class UsageCode : public TabChaine <ChaineADN,int>
  // herite de TabChaine
  // utilisation: 
  // UsageCode uc; uc.initialisation(); fichier= fopen ("cds.tfa", "rt"); 
  // uc.fichier2compte(fichier,1); // lecture par codon de la CDS
  //  uc.compte2usage(); // remplit le tableau d'usage du code
{
 public:
  const int nbreaa;      // ty. 20
  const int nbrecodons;  // ty. 64 avec les stops
  const int offset;      // ty.21, cad le nbre de cases du
                         // tableau val inutiles (lgr des mots<3)
  const char *codegenetique; // char[65], une lettre (de l'a.a.) par codon.
                       // pourra faire l'objet d'une classe, avec le code
                       // standard par defaut
  // codegenetique= "                     KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";  

  // codegenetique= "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
  double *usagecode; //usage des codons=proba d'utiliser un codon sachant l'aa
  UsageCode();
  ~UsageCode();
  void affichage() const;
  int cumuleaa(int aa) const; 
  // a partir d'un (indice de) codon, compte combien de fois l'aa
  // correspondant a ete lu en cherchant tous les codons correspondant a
  // l'acide amine donne et en cumulant leurs occurences
  void compte2usage (); 
  // pour passser du comptage (fait eg. par fichier2compte) a usagecode:
};
#endif

#ifndef CLASSprotmat
#define CLASSprotmat
class ProtMat : public TabChaine <Chaine,int>
  // herite de TabChaine
{
 public:
  const int nbreaa; // ty. 20 mais 24 pour BLOSUM62
  const int offset; // 21, cad le nbre de cases du tableau val inutiles
                    // (lgr des mots<3)
  //  ProtMat();
  //  ~ProtMat();
  // CONSTRUCTEURS:
  ProtMat(char* aa);                // use strlen
  ProtMat(char* aa, int _nbreaa);
};
#endif

#endif
