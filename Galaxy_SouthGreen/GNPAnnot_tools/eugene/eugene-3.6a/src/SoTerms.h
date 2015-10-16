
#ifndef SOTERMS_H
#define SOTERMS_H

#include <string>
#include <map>
#include "Const.h"
#include "System.h"
#include <iostream>
using namespace std; 

/**
 * class SoTerms
    
    example of enregistrement
    [Term]
    id: SO:0000000
    name: Sequence_Ontology
    subset: SOFA
    is_a: SO:... ! 

 */


class SoTerms
{
  private:
    map <string,string>  idToName_;
    map <string,string>  nameToId_;
    map <string, string> isA_; 

  public:
  
    // Constructors/Destructors
    SoTerms();
    ~SoTerms ( );
    void loadSoFile();
    void loadFile( char * filename );
    bool existsId   (string name);
    bool existsName (string name);
    int size();
    void Print ();
    string getIdFromName(string name);
    string getIsAFromId (string id);
    bool   isANcRNA     (string id);

};
#endif // SOTERMS_H
