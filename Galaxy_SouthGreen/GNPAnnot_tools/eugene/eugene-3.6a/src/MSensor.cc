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
// $Id: MSensor.cc,v 1.34 2007-08-14 08:05:31 cnoirot Exp $
// ------------------------------------------------------------------
// File:     MSensor.cc
// Contents: classes UseSensor, MasterSensor
// ------------------------------------------------------------------

#include "MSensor.h"


// reserve memory for static variables
bool                        MasterSensor::IsInitialized = false;
std::string                 MasterSensor::PluginsDir;
std::vector <std::string>   MasterSensor::LoadedSensorsList;
std::vector <std::string>   MasterSensor::MSSensorsList;
std::vector <SensorLoader*> MasterSensor::dllList;

inline bool Prior(const UseSensor *A, const UseSensor *B)
{ return(A->Priority < B->Priority); }

/*************************************************************
 **                        UseSensor                        **
 *************************************************************/
// -------------------------
//  Default constructor.
// -------------------------
UseSensor :: UseSensor ()
{
  Priority = 0;
  Name[0] = '\000';
}

// -------------------------
//  Default constructor.
// -------------------------
UseSensor :: UseSensor (int priority, char name[FILENAME_MAX+1])
{
  Priority = priority;
  strcpy(Name,name);
}




/*************************************************************
 **                       MasterSensor                      **
 *************************************************************/

extern Parameters PAR;


// ------------------------
//  Default destructor.
// ------------------------
MasterSensor :: ~MasterSensor (void)
{
  // Need to be commented because a segmentation fault sometimes occured 
  // for example in Units tests
  // delete created instances with new
//   unsigned int i;
//   for(i=0; i<theSensors.size(); i++) {
//     delete theSensors[i];
//     theSensors[i] = NULL;
//   }
  theSensors.clear();
  
//   for(i=0; i<dllList.size(); i++) {
//     delete dllList[i];
//     dllList[i] = NULL;
//   }
  dllList.clear();
}

// ------------------------
//  Init Master.
// ------------------------
void MasterSensor :: InitMaster (DNASeq *X)
{
  char *key_name;
  int  val_prior;
  int  j;
  unsigned int i;
  char* c = new char[FILENAME_MAX+1];
  std::vector <UseSensor*> msList;
  std::string use_name;

  if (!IsInitialized) {
    PluginsDir = (std::string)PAR.getC("eugene_dir")+"/"+PLUGINS_DIR+"/";

    // On récupère les couples nom de sensor/priorité du .par
    PAR.ResetIter();
    while( PAR.getUseSensor(&key_name, &val_prior) ) 
      msList.push_back( new UseSensor(val_prior, key_name) );

    // On les tri
    sort(msList.begin(), msList.end(), Prior);
    
    // Update the list of sensors defined at the top level
    for (i=0; i<msList.size(); i++)  MSSensorsList.push_back( (std::string) msList[i]->Name );
    
    // delete instances of UseSensor
    for (i=0; i<msList.size(); i++)  delete msList[i];
    
    IsInitialized = true;
  } 

  // create instance(s) of used sensors
  for (i=0; i<MSSensorsList.size(); i++) {
    use_name = MSSensorsList[i] + ".use";
    strcpy(c, use_name.c_str()); // to avoid passing a const as argument to PAR.getI
    if (PAR.probeKey(c))
      for (j=0; j<PAR.getI(c); j++){
	theSensors.push_back( MakeSensor( MSSensorsList[i].c_str(),j, X) );
	theSensors.back()->Init(X);
      }
  }
  delete [] c;
}

// ------------------------
// Create an instance of a sensor and load the .so before if necessary
// ------------------------
Sensor* MasterSensor :: MakeSensor (std::string name, int n, DNASeq *X)
{
  Sensor* s = NULL;
  int dll_index;

  dll_index = LoadSensor (name);
  s = dllList[dll_index]->MakeSensor( n, X);

  return s;
}


// ------------------------
// load the sensor s  if no yet done (j is just for print)
// return index of sensor in dllList
// ------------------------
int MasterSensor :: LoadSensor (std::string name)
{
  bool is_loaded = false;
  int dll_index = 0;
  std::string complete_name;

  // is the sensor yet loaded
  for (unsigned int i=0; i< LoadedSensorsList.size(); i++)
    if ( name == LoadedSensorsList[i] ) {
      is_loaded = true;
      dll_index = i;
      i = LoadedSensorsList.size();
    }

  if (!is_loaded) {
    LoadedSensorsList.push_back( name );
    dll_index = LoadedSensorsList.size() - 1; 
    complete_name = PluginsDir + name + ".so";
    char * nom = new char [complete_name.length()+1];
    strcpy (nom, complete_name.c_str() );
    dllList.push_back( new SensorLoader ( nom ) ); 
    delete [] nom;
    if(!dllList[dll_index]->LastError()) {
      fprintf(stderr,"Loading %.21s", name.c_str());
      for(int k=name.size(); k < 22; k++) fprintf(stderr,".");
      fprintf(stderr,"done\n");
    } else {
      fprintf(stderr,"Error: ignored plugin (invalid or not found) : %s\n",complete_name.c_str());
      exit(2);
    }
  }

  return dll_index;
}



// ------------------------
//  Init. the sensors.
// ------------------------
void MasterSensor :: InitSensors (DNASeq *X)
{
  for(int i=0; i<(int)theSensors.size(); i++)
    theSensors[i]->Init(X);
}

// --------------------------
//  Get informations at pos.
// --------------------------
void MasterSensor :: GetInfoAt (DNASeq *X, int pos, DATA *d)
{
  int i;
  
  for(i=0; i < DATA::LastSigType;  i++)
    d->sig[i].Clear();
  
  for(i=0; i< DATA::LastContentsType; i++) d->contents[i] = 0.0;
  
  d->EstMatch = 0; // WARNING temporaire : EST -> on est dans intron
  
  for(i=0; i<(int)theSensors.size(); i++) 
    theSensors[i]->GiveInfo(X, pos, d);
  
  for (i=0; i< DATA::LastSigType; i++)
    d->sig[i].SetToDefault();
}

// --------------------------
//  Print informations at pos.
// --------------------------
void MasterSensor :: PrintDataAt (DNASeq *X, int pos, DATA *d, FILE *OUT)
{
  int i,j;
  fprintf(OUT, "%6d %c", 1+pos, (*X)[pos]);

  for(i=0; i<DATA::LastContentsType; i++)
    fprintf(OUT, " %5.3g",d->contents[i]);

  for(i=0; i< Signal::LastEdge;  i++) {
    fprintf(OUT, " ||");
    for (j=0; j< DATA::LastSigType; j++)
      fprintf(OUT, " %5.3g",d->sig[j].weight[i]);
  }
  fprintf(OUT, "\n");
}

// --------------------------------------------
//  Get special info at special pos.
//  Retourne 1 si les sensors sont porteurs
//  d'infos de type "type" 0 sinon.
// --------------------------------------------
int MasterSensor :: GetInfoSpAt (unsigned char type,
				 DNASeq *X, int pos, DATA *d)
{
  int i;
  int info = 0;  // Aucune info
  
  for(i=0; i< DATA::LastSigType;  i++)
    d->sig[i].Clear();
  
  for(i=0; i<DATA::LastContentsType; i++) d->contents[i] = 0.0;
  
  d->EstMatch = 0; // WARNING temporaire : EST -> on est dans intron

  for(i=0; i<(int)theSensors.size(); i++) 
    if (theSensors[i]->type & type) {
      theSensors[i]->GiveInfo(X, pos, d);
      info = 1;
    }

  for (i=0; i< DATA::LastSigType; i++)
    d->sig[i].SetToDefault();

  return info;
}

// --------------------------
//  Post analyse the sensors.
// --------------------------
void MasterSensor :: PostAnalyse (Prediction *pred, FILE *MISC_INFO)
{
  for(int i=0; i<(int)theSensors.size(); i++)
    theSensors[i]->PostAnalyse(pred, MISC_INFO);
}

