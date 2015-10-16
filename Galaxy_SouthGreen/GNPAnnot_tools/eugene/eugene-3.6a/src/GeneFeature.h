
#ifndef GENEFEATURE_H
#define GENEFEATURE_H

#include <string>
#include <iostream>

//Feature class:
#include "Locus.h"
#include "Attributes.h"
using namespace std; 

/**
 * class GeneFeature
 */


class GeneFeature
{
  private:
    string    id_;
    string    seqid_;
    string    source_;
    string    type_;
    Locus *   locus_;
    int       length_;
    double    score_;
    char      phase_;
    Attributes * attributes_;
    bool      valid_ ;
    
    void      ParseLine( char * line, int );
    
  public:
  
    // Constructors/Destructors
    GeneFeature ( );
    GeneFeature ( char * , int); 
    GeneFeature ( const GeneFeature & gene) ;
    //GeneFeature (string seqId, string source, string type, int start, int end, double score, char strand, char phase, char * attributes);
    virtual ~GeneFeature ( );
    
    void    setValid ( bool valid);
    void    setSeqId ( string id );
    
    bool    getValid ( )  const ;
    string  getSeqId ( )  const ;
    string  getString ( ) const ;
    char    getPhase ( )  const ;
    string  getId ( )     const ;
    string  getParent ( ) const ;
    string  getType ( )   const ;
    Locus * getLocus ( )  const ;
    double  getScore ( )  const ;
    int     getLength ( ) const ;
    bool    hasTarget ( )  const ;
    Attributes * getAttributes () const ;
    
    friend ostream& operator<<(ostream& o, const GeneFeature & geneFeature )
    {
      return o <<  geneFeature.getString();
    }

};

#endif // GENEFEATURE_H
