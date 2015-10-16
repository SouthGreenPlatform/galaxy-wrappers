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
#include "SoTerms.h"
#include "Param.h"

extern Parameters PAR;

// -----------------------
//  Default constructeur
// -----------------------
SoTerms::SoTerms(  )
{
}

SoTerms::~SoTerms ( )
{
  idToName_.clear();
  nameToId_.clear();
}


void SoTerms::loadSoFile( )
{
  // Get the SoTerm file path
  char * filenameSoTerms = PAR.getC("Gff3.SoTerms", 0, 0);
  char * soTerms         = new char[FILENAME_MAX+1];
  strcpy(soTerms, PAR.getC("eugene_dir"));
  strcat(soTerms, filenameSoTerms );

  loadFile(soTerms);

  delete []soTerms;
}


void SoTerms::loadFile( char * filename )
{
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    fprintf(stderr, "Class SoTerms : cannot open start file %s\n", filename);
    exit(2);
  }
  char tag[60] , value[100];
  
  int i=1;
  int j=0;
  char  line[MAX_LINE];
  while(fp  &&  fgets (line, MAX_LINE, fp) != NULL) 
  {
    j++;
    if (line[0] == 'i' && line[1] == 'd') 
    {
      i = sscanf(line, "id: %s", &value);
      if (i > 0) 
      {
	char soId[60];
	char soName[60];
	strcpy (soId, value );
	fgets (line, MAX_LINE, fp);
	i = sscanf(line, "name: %s", &value);
	strcpy (soName, value );
	idToName_[to_string(soId)]=to_string(soName);
	nameToId_[to_string(soName)]=to_string(soId);
	const char* resFgets = " ";
	while (fp  && strncmp (line,"[Term]",6)  != 0 && fgets (line, MAX_LINE, fp)!=NULL)
	{
	 
	  //cout << "sscanf : synonym : "<<line<< "comparaison" <<strncmp (line,"[Term]",6) <<endl;
	  if (strncmp (line,"synonym",7)  == 0  )
	  {
	    
	    //synonym: "amplicon" RELATED []
	    //synonym: "PCR product" EXACT []
	    string tmpStr = to_string(line);
	    //if (tmpStr.find("EXACT") != string::npos )
	    //{
	    int begin=tmpStr.find_first_of('"')+1;
	    
	    string synomynString = tmpStr.substr(begin,tmpStr.find_last_of('"')-begin); //entre " "
	     // tmpStr = synomynString.substr(0,synomynString.find('"')); //delete string after synomyn value
	    int pos=synomynString.find(" ");
	      
	      while ( pos != string::npos  )
	      { 
		synomynString.at(pos) = '_';
		pos=synomynString.find(" ");
	      }
	     // cout << "synonym : "<<synomynString << endl;
	      nameToId_[synomynString]=to_string(soId);
	  
	  //  }
	  }
	  else if (strncmp (line,"is_a", 4)  == 0  ) // Save the isA value
	  {
             string tmpStr = to_string(line);
             int    begin  = tmpStr.find_first_of(' ') + 1;
             string isA    = tmpStr.substr(begin, tmpStr.find_last_of('!') - begin - 1);
             // Save the 'is-a' value of the current so term
             isA_[soId] = isA;
          }
	}
      }
    }
  }
  fclose(fp);
}

void SoTerms::Print (void)
{
  map<string,string>::iterator it;
  for ( it=idToName_.begin() ; it != idToName_.end(); it++ )
  {
     cout << "idToName_ Key : " << (*it).first << " Value : " << (*it).second << endl;
  }
  for ( it=nameToId_.begin() ; it != nameToId_.end(); it++ )
  {
    cout << "nameToId_ Key : " << (*it).first << " Value : " << (*it).second << endl;
  }
  for ( it=isA_.begin() ; it != isA_.end(); it++ )
  {
    cout << "isA_ Key : " << (*it).first << " Value : " << (*it).second << endl;
  }
}

bool SoTerms::existsId(string id)
{
  map<string,string>::iterator it=idToName_.find(id);
  bool res = true;
  if (it == idToName_.end()) {
    res=false;
  }
  return res;
}

bool SoTerms::existsName(string name)
{
  int pos=name.find(" ");
  while ( pos != string::npos  )
  { 
    name.at(pos) = '_';
    pos=name.find(" ");
  }
  map<string,string>::iterator it=nameToId_.find(name);
  bool res = true;
  
  if (it == nameToId_.end()) {
    res=false;
  }
  return res;
}

int SoTerms::size (void)
{
  return idToName_.size();
}

string SoTerms::getIdFromName (string name)
{
  int pos=name.find(" ");
  while ( pos != string::npos  )
  { 
    name.at(pos) = '_';
    pos=name.find(" ");
  }
  return nameToId_[name];
}

/* Return the 'is_a' value of the soterm id */
string SoTerms::getIsAFromId(string id)
{
  map<string,string>::iterator iter = isA_.find(id);
  if ( iter != isA_.end() ) 
  {
      return iter->second;
  }

  return "";
}

/* Return true if the id represents a term which is a ncRNA 
(directly the ncRNA term or derived from ncRNA) */
bool SoTerms::isANcRNA(string id)
{
	if (id == "SO:0000655") // Code So for NcRNA
		return true;
	else 
	{
		string parentId = getIsAFromId(id); // get the parent term id
		if (parentId == "")
			return false;
		else 
			return SoTerms::isANcRNA(parentId);
	}
}
