
#ifndef ATTRIBUTES_H
#define ATTRIBUTES_H

#include <string>
#include <vector>
#include <iostream>

//Feature class:
#include "Target.h"
#include "Const.h"
    
/**
  * class Attributes
  */


class Attributes
{


private:
 
  //From gff3 definition
  string    id_;
  string    name_;
  string    alias_;
  string    parent_;
  Target *  target_;
  // on an alignment list of the Gap.
//  vector < Gap * > gaps_;
  
  string derivesFrom_;
  string note_;
  string dbxref_;
  string ontologyTerm_;
  string database_;  
  //For eugene definition
  int         length_;
  int         nbExon_;

public:
  

  // Constructors/Destructors
  Attributes  ( );
  Attributes  ( std::string line ); 
  Attributes  ( const Attributes & attribut ) ;
  ~Attributes ( );
  
  void InitializeDefault ( void);
  
  void        setId ( std::string id );
  void	      setOntologyTerm ( std::string term );
  std::string getId ( ) const ;
  std::string getString ( ) const ;
  std::string getParent ( ) const ;
  std::string getOntologyTerm ( ) const ;
  Target *    getTarget ( ) const ;
  bool        hasTarget ( ) const ;
  // Puplic attribute accessor methods
  //  
  friend std::ostream& operator<<(std::ostream& o, const Attributes & attributes )
  {
	return o << attributes.getString();
  }
	

  //Should be static members , 
  static const int NB_ATTRIBUTES_NAMES_=18;
  static std::vector<std::string> attributeNames_; // "ID=","Note=","Target="
  void InitializeAttributeNames()
  {
	attributeNames_.push_back("ID");
	attributeNames_.push_back("Name");
	attributeNames_.push_back("Alias");
	attributeNames_.push_back("Parent");
	attributeNames_.push_back("Target");
	attributeNames_.push_back("Gap");
	attributeNames_.push_back("Derives_from");
	attributeNames_.push_back("Note");
	attributeNames_.push_back("Dbxref");
	attributeNames_.push_back("Ontology_term");
	attributeNames_.push_back("database");
	attributeNames_.push_back("is_full_length");
	attributeNames_.push_back("target_length");
	attributeNames_.push_back("target_sequence");
	attributeNames_.push_back("nb_exon");
	attributeNames_.push_back("length");
	attributeNames_.push_back("frame_hit");
	attributeNames_.push_back("score_hit");
  }

};


  // Private attributes
	//From GFF3 specification 
		//ID=match00001;Note=0907A18;Target=sp_P54263_SYN_THETH 6 34;Ontology_term=SO:000012;
	//From sensor : 
		//database=cDNAMt;
		//In class Target : is_full_length=1;target_length=231;
	//Eugene output: 
		//nb_exon=7;length=3598
#endif // ATTRIBUTES_H
