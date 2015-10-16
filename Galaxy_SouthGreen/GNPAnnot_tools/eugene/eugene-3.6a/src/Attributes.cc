/***************************************************************************
 *   Copyright (C) 2007 by Noirot   *
 *   cnoirot@milieu   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "Attributes.h"

//Attribut de classe;
    std::vector<std::string> Attributes::attributeNames_;
    
// -----------------------
//  Default constructeur
// -----------------------

Attributes::Attributes()
{
  InitializeDefault();
}


Attributes::Attributes ( std::string line ) 
{
  if (attributeNames_.empty()) {
    InitializeAttributeNames();
  }	
  
  InitializeDefault();

  char * oneAttributeString=new char[line.length()+1];;
  char att[MAX_GFF_LINE];
  strcpy(att,line.c_str());

  oneAttributeString = strtok (att,"=;");
    char targetString[MAX_GFF_LINE];
    strcpy (targetString,"");
    int isFullLength=-1;
    int targetLength=-1;
    int targetFrameHit=0;
    int targetScoreHit=-1;
    string targetSequence="";
    // Parse attributes 
    //  
    while ( oneAttributeString != NULL )
    {
      //Note=0907A18;
      char value[MAX_GFF_LINE]="";
      int i=0;
      while (i<NB_ATTRIBUTES_NAMES_ && strcmp (value,"") ==0 && oneAttributeString != NULL)  
      {
	// cherche dans le tableau attributeNames_ l'indice correspondand a l'attribut en cours.
	if ( strcmp (attributeNames_[i].c_str(),oneAttributeString) ==0  ) 
	{
	  oneAttributeString = strtok (NULL,"=;");
	  strcpy(value,oneAttributeString);
	  // else at end of loop, add 1, with no reason.
	  i--;
        }
	i++;
      } //end while
	
      if  ( i >= NB_ATTRIBUTES_NAMES_){
	oneAttributeString = strtok (NULL,"=;");
        continue;
      }

      //en fonction de lindice i on sait sur quel attribut on se situe. 
      if ( attributeNames_[i] == "ID" && ( strcmp(value,"")!=0 ) ) {
	id_=to_string(value);
      }
      if ( attributeNames_[i] == "Name" && ( strcmp(value,"")!=0 )){
	name_=to_string(value);
      }
      if ( attributeNames_[i] == "Alias" && ( strcmp(value,"")!=0 )){
	alias_=to_string(value);
      }
      if ( attributeNames_[i] == "Parent" && ( strcmp(value,"")!=0 )){
	parent_=to_string(value);
	
      }
      if ( attributeNames_[i] == "Target" && ( strcmp(value,"")!=0 )) {
	strcpy (targetString,value);
      }
      if ( attributeNames_[i] == "Gap" && ( strcmp(value,"")!=0 )) {
	//gap_=value;
      }
      if ( attributeNames_[i] == "Derives_from" && ( strcmp(value,"")!=0 )) {
	derivesFrom_=to_string(value);
      }
      if ( attributeNames_[i] == "Note" && ( strcmp(value,"")!=0 )){
	note_=to_string(value);
      }
      if ( attributeNames_[i] == "Dbxref" && ( strcmp(value,"")!=0 )){
	dbxref_=to_string(value);
      }
      if ( attributeNames_[i] == "Ontology_term" && ( strcmp(value,"")!=0 )) {
	ontologyTerm_=to_string(value);
      }
      if ( attributeNames_[i] == "database" && ( strcmp(value,"")!=0 )) {
	database_=to_string(value);
      }
      if ( attributeNames_[i] == "is_full_length" && ( strcmp(value,"")!=0 )) {
	isFullLength=atoi(value);
      }
      if ( attributeNames_[i] == "target_length" && ( strcmp(value,"")!=0 )) {
	targetLength=atoi(value);
      }
      if ( attributeNames_[i] == "target_sequence" && ( strcmp(value,"")!=0 )) {
	targetSequence=to_string(value);
      }
      if ( attributeNames_[i] == "nb_exon" && ( strcmp(value,"")!=0 )) {
	nbExon_= atoi(value);
      }
      if ( attributeNames_[i] == "length" && ( strcmp(value,"")!=0 )) {
	length_= atoi(value);
      }
      if ( attributeNames_[i] == "frame_hit" && ( strcmp(value,"")!=0 )) {
	targetFrameHit= atoi(value);
      }
      if ( attributeNames_[i] == "score_hit" && ( strcmp(value,"")!=0 )) {
	targetScoreHit= atoi(value);
      }
      oneAttributeString = strtok (NULL,"=;");
    } //end for all substring
    
    // check target attributes:
    if ( strcmp(targetString,"")!=0 ) { 
      char targetName[MAX_LINE];
      char targetStart[MAX_LINE],targetEnd[MAX_LINE],targetStrand=' ';
      //name
      char* targetAtt = strtok (targetString," ");
      strcpy(targetName,targetAtt);
      //start
      targetAtt = strtok (NULL," ");
      strcpy(targetStart,targetAtt);
      //end
      targetAtt = strtok (NULL," ");
      strcpy(targetEnd,targetAtt);
      //end
      targetAtt = strtok (NULL," ");
      if (targetAtt != NULL)
      {
	targetStrand=targetAtt[0];
      }
      target_ = new Target (targetName, targetSequence, targetLength, isFullLength, atoi(targetStart), atoi(targetEnd), targetStrand,targetFrameHit,targetScoreHit);
    }
    delete [] oneAttributeString;
}

Attributes::Attributes ( const Attributes &  attribut ) 
{
  id_=attribut.id_;
  name_=attribut.name_;
  alias_=attribut.alias_;
  note_=attribut.note_;
  parent_=attribut.parent_;
  if (attribut.target_ != NULL )
  {
    target_=new Target( attribut.target_);
  }
  database_=attribut.database_;
  derivesFrom_=attribut.derivesFrom_;
  ontologyTerm_=attribut.ontologyTerm_;
  length_=attribut.length_;
  nbExon_=attribut.nbExon_;
}


// -----------------------
//  Destructeur
// -----------------------
Attributes::~Attributes ( ) 
{ 
  if  ( target_ != NULL )
  {
    delete target_;
  }
}
//Initialize
void Attributes::InitializeDefault ( void)
{
  id_="";
  name_="";
  alias_="";
  note_="";
  parent_="";
  target_=NULL;
  database_="";
  derivesFrom_="";
  ontologyTerm_="";
  length_=-1;
  nbExon_=-1;
  //gaps_ = new std::vector();
  
}

// -----------------------
//Accessor
// -----------------------

// Id
//
void Attributes::setId ( std::string id ) 
{
	id_ = id;
}

void Attributes::setOntologyTerm ( std::string term )
{
	ontologyTerm_ = term;
}

std::string Attributes::getId ( ) const 
{
	return id_;
}

std::string Attributes::getParent ( ) const 
{
  return parent_;
}

std::string Attributes::getOntologyTerm () const 
{
  return ontologyTerm_;
}

Target * Attributes::getTarget () const 
{
  return target_;
}

bool Attributes::hasTarget ( ) const 
{
  bool res=false;
  if (target_ != NULL )
  {
    if ( target_->hasLocus() )
        res=true;
  }
  return res;
}

//getString
std::string Attributes::getString ( ) const 
{	
	std::string my_attributes ="";
	if ( id_ != "" ) {
		my_attributes +=attributeNames_[0]+"="+id_+";";
	}
	if ( name_ != "" ) {
	  my_attributes +=attributeNames_[1]+"="+name_+";";
	}
	if ( alias_ != "" ) {
	  my_attributes +=attributeNames_[2]+"="+alias_+";";
	}
 	if ( parent_ != ""  ) {
	  my_attributes +=attributeNames_[3]+"="+parent_+";";
 	}
	if ( target_ != NULL ) {
		my_attributes +=target_->getString();
	}
/*	if ( !gaps_.empty() ) {
		//my_attributes +=gaps_.getString();
	}*/
	if ( derivesFrom_ != ""  ) {
	  my_attributes +=attributeNames_[6]+"="+derivesFrom_+";";
	}
	if ( note_ != "" ) {
	  my_attributes +=attributeNames_[7]+"="+note_+";";
	}
	if ( dbxref_ != "" ) {
	  my_attributes +=attributeNames_[8]+"="+dbxref_+";";
	}
	if ( ontologyTerm_ != "" ) {
	  my_attributes +=attributeNames_[9]+"="+ontologyTerm_+";";
	}
	if ( database_ != "" ) {
	  my_attributes +=attributeNames_[10]+"="+database_+";";
	}
	if ( nbExon_ != -1 ) {
	  my_attributes +=attributeNames_[14]+"="+to_string(nbExon_)+";";
	}
	if ( length_ != -1 ) {
	  my_attributes +=attributeNames_[15]+"="+to_string(length_)+";";
	}
	return my_attributes;
	
} 

