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
// $Id: Sensor.NcRNA.cc,v 1.3 2010-02-01 16:18:41 sallet Exp $
// ------------------------------------------------------------------
// File:     Sensor.NcRNA.cc
// Contents: Sensor
// ------------------------------------------------------------------

#include "Sensor.NcRNA.h"
#include "../../MSensor.h"


extern MasterSensor *MS;
extern Parameters   PAR;

/*************************************************************
 **                        SensorNcRNA                        **
 *************************************************************/
// ----------------------
//  Default constructor.
// ----------------------
SensorNcRNA :: SensorNcRNA (int n, DNASeq *X) : Sensor(n)
{ 
}

// ----------------------
//  Default destructor.
// ----------------------
SensorNcRNA :: ~SensorNcRNA ()
{
}

// --------------------------------
//  Init ncRNA
// --------------------------------
void SensorNcRNA :: Init (DNASeq *X)
{
  // Read parameters
  strcpy(tStartNpc, PAR.getC("NcRNA.TStartNpc*", GetNumber()));
  strcpy(tStopNpc,  PAR.getC("NcRNA.TStopNpc*",  GetNumber()));
  strcpy(npcRna,    PAR.getC("NcRNA.NpcRna*",    GetNumber()));

  char tempname[FILENAME_MAX+1];
  fileExt      = PAR.getC("NcRNA.FileExtension",    GetNumber());
  inputFormat_ = to_string(PAR.getC("NcRNA.format", GetNumber(),1));
  
  fprintf(stderr, "Reading %s file....", fileExt);
  fflush(stderr);
  strcpy(tempname, PAR.getC("fstname"));
  strcat(tempname, ".");
  strcat(tempname, fileExt);

  if ( inputFormat_ == "GFF3" )
  {
    strcat(tempname,".gff3");
    GeneFeatureSet * geneFeatureSet = new GeneFeatureSet (tempname);
    ReadNcRNAGff3(*geneFeatureSet,  X->SeqLen);
    delete geneFeatureSet;
  }
  else 
  {
	fprintf(stderr, "Error in NcRNA format: GFF3 file is required.\n");
	exit(2);
  }

  fprintf(stderr, "done\n");
  fflush(stderr);	

}

// -----------------------
//  GiveInfo.
// -----------------------
void SensorNcRNA:: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  int i, j;
  float k;

  // Get signals information
  if (!vSig.empty()) 
  {
    iSig = 0;
    while (iSig < (int)vSig.size() )
    {
      if ((vSig[iSig]->pos == pos)) {
        i = vSig[iSig]->type;         // signal type
        j = vSig[iSig]->edge;         // ~ strand
        k = atof(vSig[iSig]->score);  // Score ]-oo +oo)[
        d->sig[i].weight[j] += k;
      }
      iSig++;
    }
  }

  // Get content information
  if(!vCon.empty()) 
  {
    iCon = 0;
    while(iCon < (int)vCon.size()  &&  (pos > vCon[iCon]->end || pos < vCon[iCon]->start) ) iCon++;
    while ( iCon < (int)vCon.size() )
    {
	if ( (pos >= vCon[iCon]->start)  &&  (pos <= vCon[iCon]->end) )
	  d->contents[vCon[iCon]->type] += *(vCon[iCon]->score);
        iCon++;
    }
  }
}


// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorNcRNA :: Plot(DNASeq *TheSeq)
{}

// ------------------
//  Post analyse
// ------------------
void SensorNcRNA :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}

void SensorNcRNA :: ReadNcRNAGff3(GeneFeatureSet & geneFeatureSet, int len)
{
  int   startC, endC, edge;   // C -> content
  int   startS, endS;         // S -> signal       need for reverse
  char  strand;
  string idSo;
  char  feature[50];
  vector<GeneFeature *>::iterator it = geneFeatureSet.getIterator();
  int nbGeneFeature                  = geneFeatureSet.getNbFeature();
  int j = 0;

  for ( int i=0 ; i < nbGeneFeature ; i++, it++ )
  {
    j++;
    startC = (*it)->getLocus()->getStart();
    endC   = (*it)->getLocus()->getEnd();
    strcpy (feature, (*it)->getType().c_str());
    strand = (*it)->getLocus()->getStrand();

    idSo   = (*it)->getType(); //recup code SO (third column of gff3)
    if ( idSo.find("SO:") == string::npos ) // If its not the SO code
    {
      string tmp = GeneFeatureSet::soTerms_->getIdFromName(idSo);
      idSo       = tmp;
    }
    /* Strand ? */
    if (strand == '+')
    {
      edge    = 0;
      startS = startC;
      endS   = endC;
    }
    else if (strand == '-')
    {
      edge   = 1;
      startS = endC   + 1;
      endS   = startC - 1;
    }
    else 
    {
      fprintf(stderr, "WARNING: feature %s line %d strand unknown"
	  " => ignored.\n", feature, j);
      continue;
    }
    // check it is a ncRNA
    if ( !GeneFeatureSet::soTerms_->isANcRNA(idSo) )
    {
	fprintf(stderr, "WARNING: feature %s line %d is not a non coding protein RNA"
	  " => ignored.\n", feature, j);
      continue;
    }
    startC--;
    endC--;

    vCon.push_back( new Contents(startC, endC, DATA::RNAF+edge, new float(atof(npcRna))));
    vSig.push_back( new Signals (startS-1, DATA::tStartNpc, edge, tStartNpc) );
    vSig.push_back( new Signals (endS,     DATA::tStopNpc,  edge, tStopNpc)  );
  }
}
