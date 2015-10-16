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
// $Id: Sensor.GSplicer.cc,v 1.16 2009-01-12 14:33:28 sallet Exp $
// ------------------------------------------------------------------
// File:     Sensor.GSplicer.cc
// Contents: Sensor GSplicer  
// ------------------------------------------------------------------

#include "Sensor.GSplicer.h"

extern Parameters PAR;

/*************************************************************
 **                       SensorGSplicer
 *************************************************************/

// ----------------------
// Default constructor.
// ----------------------
SensorGSplicer :: SensorGSplicer (int n, DNASeq *X) : Sensor(n)
{
  char tempname[FILENAME_MAX+1];

  type = Type_Acc|Type_Don;


  fflush(stderr);
  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".Gsplicer");

  
  inputFormat_ = to_string(PAR.getC("GSplicer.format", GetNumber(),1));

  if ( inputFormat_ == "GFF3" )
  {
    fprintf(stderr, "Reading splice site file (GeneSplicer)........forward, reverse ");
      strcat(tempname,".gff3");
      ReadGSplicerGff3(tempname, X->SeqLen);
  }
  else
  {
    fprintf(stderr, "Reading splice site file (GeneSplicer)........forward, reverse ");
    ReadGSplicer(tempname, X->SeqLen);
  }
  fflush(stderr);
 
  fprintf(stderr,"done\n");
  CheckSplices(X,vPosAccF, vPosDonF, vPosAccR, vPosDonR);
}

// ----------------------
//  Default destructor.
// ----------------------
SensorGSplicer :: ~SensorGSplicer ()
{
  // Clear the data structures 
  vPosAccF.clear();  vPosAccR.clear();
  vPosDonF.clear();  vPosDonR.clear();
  vValAccF.clear();  vValAccR.clear();
  vValDonF.clear();  vValDonR.clear();
}

// ----------------------
//  Init.
// ----------------------
void SensorGSplicer :: Init (DNASeq *X)
{
  coefAcc = PAR.getD("GSplicer.coefAcc*",GetNumber());
  penAcc  = PAR.getD("GSplicer.penAcc*", GetNumber());
  coefDon = PAR.getD("GSplicer.coefDon*",GetNumber());
  penDon  = PAR.getD("GSplicer.penDon*", GetNumber());
  
  iAccF = iDonF = 0;
  iAccR = iDonR = 0;
  PositionGiveInfo = -1;

  if (PAR.getI("Output.graph")) Plot(X);
}

// -------------------------------------
//  Read GeneSplicer GFF3 file.
// -------------------------------------
void SensorGSplicer :: ReadGSplicerGff3(char name[FILENAME_MAX+1], int SeqLen)
{
  GeneFeatureSet * geneFeatureSet = new GeneFeatureSet (name); 

  vector< GeneFeature *>::iterator it = geneFeatureSet->getIterator();
  int nbElement=geneFeatureSet->getNbFeature();
  //geneFeatureSet->printFeature();
  int i=0;
  while ( i<nbElement )
  {
    //(*it)->second();
    GeneFeature * tmpFeature = *it;
    string idSo=tmpFeature->getType();
    if ( idSo.find("SO:") == string::npos )
    {
      string tmp=GeneFeatureSet::soTerms_->getIdFromName(idSo);
      idSo=tmp;
    }
     // Forward
    
    if ( tmpFeature->getLocus()->getStrand() == '+' ) 
    {
    
      if (idSo == "SO:0000163") //donor
      {
	vPosDonF.push_back( tmpFeature->getLocus()->getStart()-1 );
	vValDonF.push_back( tmpFeature->getScore() );
	
      }
      else 
      {
	if (idSo == "SO:0000164")  //acceptor
	{
	  vPosAccF.push_back(tmpFeature->getLocus()->getEnd() );
  	  vValAccF.push_back( tmpFeature->getScore() );
	}
      }
    }
    // Reverse
    if ( tmpFeature->getLocus()->getStrand() == '-' ) {
      if(idSo == "SO:0000163") {
	vPosDonR.push_back( tmpFeature->getLocus()->getEnd() );
	vValDonR.push_back( tmpFeature->getScore() );
      }
      else {
	if (idSo == "SO:0000164")  //acceptor
	{
	  vPosAccR.push_back( tmpFeature->getLocus()->getStart()-1 );
	  vValAccR.push_back( tmpFeature->getScore() );
	}
      }
    }
    it++;
    i++;
  }
  delete geneFeatureSet;
}

// ------------------------
//  GiveInfo.
// ------------------------
void SensorGSplicer :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  bool update = false;

  if ( (PositionGiveInfo == -1) || (pos != PositionGiveInfo+1) ) update = true; // update indexes on vectors
  PositionGiveInfo = pos;

  // Accepteur Forward
  if(!vPosAccF.empty()) {
    if (update) 
      iAccF = lower_bound(vPosAccF.begin(), vPosAccF.end(), pos)-vPosAccF.begin();
    
    if((iAccF<(int)vPosAccF.size()) && (vPosAccF[iAccF] == pos)) {
      d->sig[DATA::Acc].weight[Signal::Forward] += (vValAccF[iAccF]*coefAcc) - penAcc;
      iAccF++;
    }
  }

  // Accepteur Reverse
  if(!vPosAccR.empty()) {
    if (update) 
      iAccR = lower_bound(vPosAccR.begin(), vPosAccR.end(), pos)-vPosAccR.begin();
    
    if((iAccR<(int)vPosAccR.size()) && (vPosAccR[iAccR] == pos)) {
      d->sig[DATA::Acc].weight[Signal::Reverse] += (vValAccR[iAccR]*coefAcc) - penAcc;
      iAccR++;
    }
  }
 
  // Donneur Forward
  if(!vPosDonF.empty()) {
    if (update)
      iDonF = lower_bound(vPosDonF.begin(), vPosDonF.end(), pos)-vPosDonF.begin();

    if ((iDonF<(int)vPosDonF.size()) && (vPosDonF[iDonF] == pos)) {
      d->sig[DATA::Don].weight[Signal::Forward] += (vValDonF[iDonF]*coefDon) - penDon;
      iDonF++;
    }
  }
    
  // Donneur Reverse
  if(!vPosDonR.empty()) {
    if (update) 
      iDonR = lower_bound(vPosDonR.begin(), vPosDonR.end(), pos)-vPosDonR.begin();

    if ((iDonR<(int)vPosDonR.size()) && (vPosDonR[iDonR] == pos)) {
      d->sig[DATA::Don].weight[Signal::Reverse] += (vValDonR[iDonR]*coefDon) - penDon;
      iDonR++;
    }
  }

}

// ----------------------------
//  Plot Sensor information.
// ----------------------------
void SensorGSplicer :: Plot(DNASeq *X)
{
  for (int i =0; i < (int)vPosAccF.size(); i++)
    PlotAcc(vPosAccF[i],1,Norm(log(vValAccF[i]),20.0));
  
  for (int i =0; i < (int)vPosDonF.size(); i++)
    PlotDon(vPosDonF[i],1,Norm(log(vValDonF[i]),20.0));
  
  for (int i =0; i < (int)vPosAccR.size(); i++)
    PlotAcc(vPosAccR[i],-1,Norm(log(vValAccR[i]),20.0));

  for (int i =0; i < (int)vPosDonR.size(); i++)
    PlotDon(vPosDonR[i],-1,Norm(log(vValDonR[i]),20.0));
}

// ---------------------
//  Plot normalisation.
// ---------------------
double SensorGSplicer :: Norm(double x, double n) {
  return (((n)+(Max(-(n),x)))/(n));
}

// ------------------
//  Post analyse.
// ------------------
void SensorGSplicer :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}

// -------------------------------------
//  Read GeneSplicer file.
// -------------------------------------
void SensorGSplicer :: ReadGSplicer(char name[FILENAME_MAX+1], int SeqLen)
{
  FILE *fp;
  char type[10];
  int i,j,nt1,nt2;
  double score;
 
  if (!(fp = fopen(name, "r"))) {
    fprintf(stderr, "cannot open splice sites file %s\n",  name);
    exit(2);
  }
  
  i=0;
  while(!feof(fp)){
    i++;
    j = fscanf(fp,"%d %d %lf %*s %s\n", &nt1, &nt2, &score, type);
    
    if (j < 4) {
      fprintf(stderr, "\nError in splice sites file! %s, line %d\n", name, i);
      exit(2);
    }
    // Forward
    if ( (nt2-nt1) == 1 ) {
      if(type[0]=='d') {
	vPosDonF.push_back( nt1-1 );
	vValDonF.push_back( score );
      }
      else {
	vPosAccF.push_back( nt2 );
	vValAccF.push_back( score );
      }
    }
    // Reverse
    if ( (nt2-nt1) == -1 ) {
      if(type[0]=='d') {
	vPosDonR.push_back( nt1 );
	vValDonR.push_back( score );
      }
      else {
	vPosAccR.push_back( nt2-1 );
	vValAccR.push_back( score );
      }
    }
  }
  fclose(fp);
}
