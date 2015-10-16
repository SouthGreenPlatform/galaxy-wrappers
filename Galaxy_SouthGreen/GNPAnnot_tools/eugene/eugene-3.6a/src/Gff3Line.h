/** \file   Gff3line
 *  \author Sebastien Letort
 *  \date   10 novembre 2006
 *  Ce fichier contient tout ce qu'il faut pour gerer les lignes du gff3
 *  On trouve les enumerations des codes SOFA et SO,
 *  la classe Sofa, qui sert a g�rer les codes
 *  la classe Gff3Line et sa classe interne Attributes
 */

#ifndef __GFF3Line_H__
#define __GFF3Line_H__

//les bibliotheques de la STL
#include <vector>
#include <string>
#include <fstream>  //ofstream
#include <iomanip>  //settfill et setw

//nos bibliotheques
#include "Prediction_cte.h"  //acces aux constantes
#include "System.h"  //to_string()

//sofa->j'ai pas trouv� de raison de mettre des entiers plutot que des chaines
//sauf l'enum
/** \enum
 *  la liste des codes SOFA et des codes SO (moins ceux deja present dans sofa)
 *  \link http://song.cvs.sourceforge.net/song/ontology/sofa.obo?view=log
 *  \link http://song.cvs.sourceforge.net/song/ontology/so.obo?view=log
 */
enum
{
  SOFA_GENE= 704,           ///<  un gene
  SOFA_MRNA = 234,          ///<  ARNm episse "does not contain introns"
  SOFA_EXON = 147,          ///<  un exon de l'ARNm episse
  SOFA_INTRON = 188,        ///<  un intron de l'ARNm
  SOFA_INTERGEN = 605,      ///<  region entre 2 genes
  SOFA_5_UTR = 204,         ///<  5'UTR (de l'ARNm episse)
  SOFA_3_UTR = 205,         ///<  UTR3' (de l'ARNm episse)
  SOFA_CDS = 316,            ///<  un CDS avec un start ET un stop
  SOFA_NCRNA      = 655,            ///< non protein coding rna

  SO_5_EXON = 200,       ///<  un exon contenant CDS/UTR 5'
  SO_3_EXON = 202,       ///<  un exon contenant  CDS/UTR 3'
  SO_5_CODING_EXON = 196,       ///<  un exon contenant CDS 5'
  SO_3_CODING_EXON = 197,       ///<  un exon contenant  CDS 3'
  SO_5_NONCODING_EXON = 445,    ///<  un exon de 5'UTR (du pre-ARNm)
  SO_3_NONCODING_EXON = 444,    ///<  un exon d'UTR3' (du pre-ARNm)
  SO_SINGLE_EXON = 5845,        ///<  un exon 'solitaire'
  SO_INTERIOR_EXON = 4,         ///<  un exon interieur codant (partie d'un CDS)
  SO_UTR_INTRON = 446          ///<  un intron d'UTR

};

/**  \class  Sofa
 *  \brief  Cette classe doit permettre de traiter les codes Sofa
 *  Pour le moment, sont regroupees ici les fonctions qui permettent de gerer
 *  les codes sofa. Du coup elles sont statiques !
 *  \todo    Faire ça proprement !!!
 */
class Sofa
{
  public :
  Sofa(){}
  ~Sofa(){}
  ///renvoie le nom correspondant au code sofa
  static std::string getName(unsigned int code_so);
  ///renvoie la chaine "SO:XXXXXXX" avec xxx le code sofa
  static std::string codeToString(const int code_so);
};

/** \class  GFF3Line
 *  \author  Sebastien Letort
 *  \date    13 novembre 2006
 *  \brief  cette classe sert a ecrire une ligne pour le GFF3
 *  \todo    A terme cette classe doit permettre de lire une ligne gff3
 *  la gestion de la classe interne Attribute est encore en discution
 */
class Gff3Line
{
  private :
//definition de la classe attribute nested class ou inner class ?
// cette classe n'a pas de raison d'exister en dehors d'une ligne de gff3
    /** \class  Attribute
     *  \brief  Pour le moment cette classe est une simple string
     *  \todo   le faire bien, doit changer quand la lecture sera codee
     */
    class Attribute
    {
      std::string liste_;
      
      public :
        Attribute(std::string attribut)
        {  
          liste_ = attribut;  
        }
        ~Attribute()
        {}

        ///ajoute un attribut a la liste
        void addAttribute(std::string attribut)
        {  
          liste_ += ";" + attribut;  
        }
        
        ///renvoie la liste d'attributs sous forme d'une chaine
        std::string str()
        {
          return liste_;  
        }
    };

//les attributs:
    std::string seq_id_;
    std::string source_;
    int type_sofa_;
    int start_, end_;
    float score_;
    char strand_;
//    unsigned short phase_;
    int phase_;  //'int' pour gerer le INDEFINI
    Attribute* attr_;

    /** \enum
     *  cette enumeration montre l'ordre d'apparition des elements dans la ligne
     *  utilise par la fonction print uniquement
     */
    enum elements
    {
      SEQ_ID = 0,
      SOURCE = 1,
      TYPE = 2,
      START = 3,
      END = 4,
      SCORE = 5,
      STRAND = 6,
      PHASE = 7,
      ATTRIBUTES = 8,
      NB_ELEMENTS_GFF3 = 9
    };

  public :
//les constructeurs
    ///constructeur pour une sequence donnee
    Gff3Line(std::string seq_id="");
    //ctr par recopie, recopie tous sauf les attributs
    Gff3Line(const Gff3Line& l)
    {
      seq_id_ = l.seq_id_;  source_ = l.source_;
      type_sofa_ = l.type_sofa_;
      start_ = l.start_;  end_ = l.end_;
      score_ = l.score_;  strand_ = l.strand_;
      phase_ = l.phase_;
    }
//le destructeur
    ~Gff3Line()
    {
      if(attr_)
        delete attr_;
    }

//les methodes
    ///\todo  faire une fn de verif permettrait de ne pas initialiser des param
    ///\todo  de facon inadequate (type doit etre != de '.')
    void check() const
    {}

    void print(std::ofstream& out, bool sofa_name=true) const;

    int getEnd() const
    { return end_;  }
    int getType() const
    { return type_sofa_;  }

    //la premiere ligne obligatoire
    static std::string header()
    {  return "##gff-version 3\n";  }
    //la premiere ligne obligatoire, suivi de la 2e
    static std::string header(
            char* seq, unsigned int start, unsigned int end)
    {
      std::ostringstream buff;
      buff  << "##gff-version 3\n"
            << "##sequence-region " << seq
            << " " << to_string(start)
            << " " << to_string(end)
            << std::endl;
      return buff.str();
    }
    static std::string endComplexFeature()
    {  return "###\n";  }

//les methodes inline pour l'initialisation
    ///affecte seq_id_
    void setSeqId(std::string str)
    {  seq_id_ = str;  }

    ///affecte le type_sofa_
    void setType(int type)
    {  type_sofa_ = type;  }

    /// affecte le champs start
    void setStart(unsigned int start)
    {  start_ = start;  }

    /// affecte le champs end
    void setEnd(unsigned int end)
    {  end_ = end;  }

    /// affecte le champs score
    void setScore(float score)
    {  score_ = score;  }

    /// affecte le champs strand
    void setStrand(char strand)
    {  strand_ = strand;  }

    /// affecte le champs phase
    void setPhase(char phase)
    {  phase_ = phase;  }

    /// affecte le champs attributes, provisoire
    /// \todo  faire une fonction plus contrainte avec un objet Attribute ?
    void setAttribute(std::string str)
    {  attr_ = new Attribute(str);  }  //verifier la chaine ?

    /// ajoute un attribut au champ attributes
    /// \todo  faire une fonction plus contrainte ?
    void addAttribute(std::string str)
    {  attr_->addAttribute(str);  }
};

#endif //__GFF3Line_H__
