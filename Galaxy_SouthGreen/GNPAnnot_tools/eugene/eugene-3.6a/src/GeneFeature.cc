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
#include "GeneFeature.h"
#include "GeneFeatureSet.h"
// -----------------------
//  Default constructeur
// -----------------------
GeneFeature::GeneFeature()
{
  id_=".";
  seqid_=".";
  source_=".";
  type_=".";
  locus_=NULL;
  score_=0.0;
  phase_='.';
  attributes_=NULL;
  valid_= true ;
  length_=0;
}


//Constructor from a gff3 line
// prerequis : no \n at the end of the line.
GeneFeature::GeneFeature ( char * line, int lineno) 
{
  id_=".";
  seqid_=".";
  source_=".";
  type_=".";
  locus_=NULL;
  score_=0.0;
  phase_='.';
  attributes_=NULL;
  valid_= true ;
  length_=0;
  
  string line_str =to_string(line);
  //chomp : 
  size_t pos=line_str.find("\n") ;
  if (pos != string::npos)
  {
    line[pos]='\0';
  }
  
  if ( line[0] != '#' )
  {  
    ParseLine(line, lineno); 
  }
}

//Constructor from a gff3 line

GeneFeature::GeneFeature ( const GeneFeature & gene)
{
  id_=gene.id_;
  seqid_=gene.seqid_;
  source_=gene.source_;
  type_=gene.type_;
  if (gene.locus_ != NULL)
  {
    locus_=new Locus( *gene.locus_ );
  }
  score_=gene.score_;
  phase_=gene.phase_;
  if (gene.attributes_ != NULL)
  {
    attributes_=new Attributes ( *gene.attributes_);
  }
  valid_= gene.valid_ ;
  length_=gene.length_;
//   fprintf (stderr,"%p end constructeur par recopy \n",this);
//   fflush (stderr);

}
// ---------------------------------------
//  ParseLine : Parse a gff3 line with a token, check SO/SOFA code and Parent Attribut.
// ---------------------------------------

void GeneFeature::ParseLine ( char * line, int linenum ) 
{
  
  int size = strlen(line);
  char * token;
  char delims[] = "\t";
  token = strtok (line,delims);
  int i =0;
  int start, end;
  char strand;
  attributes_=NULL;
  locus_=NULL;
  length_=0;
  while ( token != NULL && i<9)
  {
    //cout << "Token "<<i<< " valeur: "<<token<<endl;
    switch (i)
    {
      case 0 : 
      { //sequence id
	seqid_= to_string (token);
	break;
      }
      case 1 : 
      { //source
	source_= to_string (token);
	break;
      }
      case 2 : 
      { //type : SOFA/SO
	type_= to_string (token);
	if (  ( ! GeneFeatureSet::soTerms_->existsId(type_) ) && ( !GeneFeatureSet::soTerms_->existsName(type_) ) )  
	{
	  cout << "WARNING : " << type_ << " is not SO/SOFA referenced terms in line " << linenum << endl;
	  valid_=false;
	}
	break;
      }
      case 3 : 
      {
	start= atoi (token);
	break;
      }
      case 4 : 
      {
	end= atoi (token);
	break;
      }
      case 5 : 
      {
	score_= atof (token);
	break;
      }
      case 6 : 
      {
	char * tmp= new char [MAX_LINE];
	strcpy (tmp, token);
	strand= tmp[0];
	delete [] tmp;
	break;
      }
      case 7 : 
      {
	char * tmp= new char [MAX_LINE];
	strcpy (tmp, token);
	phase_= tmp[0];
	delete [] tmp;
	break;
      }
      case 8 : 
      {
	//cout << "Attributes : " << token <<endl;
	char * attributeString = new char[MAX_GFF_LINE];
	strcpy (attributeString , token );
	attributes_ = new Attributes (attributeString);
	id_=attributes_->getId();
	if ( id_ == "" )
	{
	  cout <<"WARNING : Feature has no ID in line " << linenum <<endl;
	  valid_=false;
	}
	delete [] attributeString;
	break;
      }
    }

    token = strtok (NULL,delims);
    i++;
  }
  locus_  = new Locus (start, end, strand);
  length_ = end-start+1;
  if (i < 9) {
    cout <<"WARNING : incomplete GFF feature in line " << linenum <<endl;
    valid_ = false;
  }
}



// -----------------------
//  Destructeur
// -----------------------
GeneFeature::~GeneFeature ( ) 
{ 
  if (locus_ != NULL)
  {
    delete locus_;
  }
  if (attributes_ != NULL)
  {
    delete attributes_;
  }
}

// -----------------------
// Accessor
// -----------------------

// Sequence Id
//
void GeneFeature::setSeqId ( string seqid ) 
{
  seqid_ = seqid;
}
string GeneFeature::getSeqId ( ) const 
{
  return seqid_;
}

// Id of feature , copy from attributes
//

string GeneFeature::getId ( ) const 
{
  return id_;
}

// Parent attributes get from attributes class if exists
//
string GeneFeature::getParent() const 
{
  string res="";
  if (attributes_ != NULL )
  {
    res= attributes_->getParent();
  }
  return res;
}

// Type
//
string GeneFeature::getType() const 
{
  return type_;
}

double GeneFeature::getScore() const 
{
  return score_;
}
// Valid attribut.
//
void GeneFeature::setValid (bool valid)
{
  valid_= valid;
}
bool GeneFeature::getValid ( ) const 
{
  return valid_;
}
// ---------------------------------------
//  getString : Build a string gff3 format
// ---------------------------------------
string GeneFeature::getString ( ) const 
{
  string geneFeature = seqid_+"\t"+source_+"\t"+type_+"\t";
  if ( locus_ != NULL )
  {
    geneFeature += to_string(locus_->getStart()) + "\t" + to_string(locus_->getEnd()) + "\t" + to_string(score_) + "\t" + to_string(locus_->getStrand()) + "\t";
  }
  else
  {
    geneFeature += ".\t.\t" + to_string(score_) + "\t.\t";
  }
  geneFeature += to_string(phase_)+"\t";
  if ( attributes_ != NULL )
  {
    geneFeature +=attributes_->getString();
  }
  return geneFeature;
} 

//Locus
Locus * GeneFeature::getLocus ( )  const 
{
  return locus_;
}

char GeneFeature::getPhase ( ) const 
{
  return phase_;
}

int GeneFeature::getLength ( ) const 
{
  return (length_);
}

Attributes * GeneFeature::getAttributes ( ) const 
{
  return (attributes_);
}

bool    GeneFeature::hasTarget ( )  const 
{
    bool res=false;
    if (attributes_ != NULL)
    {
      if ( attributes_->hasTarget() )
      {
	res=true;
      }
    }
    return res;
}
