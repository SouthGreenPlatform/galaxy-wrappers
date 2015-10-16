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
// $Id: Sensor.SPred.cc,v 1.33 2009-01-12 15:01:36 sallet Exp $
// ------------------------------------------------------------------
// File:     Sensor.SPred.cc
// Contents: Sensor SPred
// ------------------------------------------------------------------

#include "Sensor.SPred.h"
 
extern Parameters PAR;

#define NORM(x,n) (((n)+(Max(-(n),x)))/(n))

/*************************************************************
 **                     SensorSplicePredictor               **
 *************************************************************/

// ----------------------
//  Default constructor.
// ----------------------
SensorSPred :: SensorSPred (int n, DNASeq *X) : Sensor(n)
{
  char tempname[FILENAME_MAX+1];

  type = Type_Acc|Type_Don;
  
  fprintf(stderr, "Reading splice site file (Splice Predictor)...");  
  fflush(stderr);
  strcpy(tempname,PAR.getC("fstname"));
  strcat(tempname,".spliceP");
  
  inputFormat_ = to_string(PAR.getC("SPred.format", GetNumber(),1));

  if ( inputFormat_ == "GFF3" )
  {
    strcat(tempname,".gff3");
    ReadSPredGff3(tempname, X->SeqLen);
    fprintf(stderr,"forward, reverse done\n");
    fflush(stderr);
    
  }
  else
  {
    
    ReadSPredF(tempname, X->SeqLen);
    fprintf(stderr,"forward,");
    fflush(stderr);
  
    strcpy(tempname,PAR.getC("fstname"));
    strcat(tempname,".splicePR");
    ReadSPredR(tempname, X->SeqLen);
    fprintf(stderr," reverse done\n");
  }
  CheckSplices(X,vPosAccF, vPosDonF, vPosAccR, vPosDonR);
  // vectors for reverse are put in the increasing order
  reverse(vPosAccR.begin(), vPosAccR.end()); 
  reverse(vValAccR.begin(), vValAccR.end());
  reverse(vPosDonR.begin(), vPosDonR.end());
  reverse(vValDonR.begin(), vValDonR.end());
}

// ----------------------
//  Default destructor.
// ----------------------
SensorSPred :: ~SensorSPred ()
{
  vPosAccF.clear();  vPosAccR.clear();
  vPosDonF.clear();  vPosDonR.clear();
  vValAccF.clear();  vValAccR.clear();
  vValDonF.clear();  vValDonR.clear();
}

// ----------------------
//  Init SPred.
// ----------------------
void SensorSPred :: Init (DNASeq *X)
{
  accP = PAR.getD("SPred.accP*",GetNumber());
  accB = PAR.getD("SPred.accB*",GetNumber());
  donP = PAR.getD("SPred.donP*",GetNumber());
  donB = PAR.getD("SPred.donB*",GetNumber());

  iAccF = iDonF = iAccR = iDonR = 0;
  PositionGiveInfo = -1;

  if (PAR.getI("Output.graph")) Plot(X);
}

// -------------------------------------
//  Read Splice Predictor forward file.
// -------------------------------------
void SensorSPred :: ReadSPredF(char name[FILENAME_MAX+1], int SeqLen)
{
  FILE *fp;
  char buf[FILENAME_MAX];
  double par[3],strength;
  char *type;
  int i=1,j, k;
  int prevkA = -1, prevkD = -1;
  
  fp = FileOpen(NULL,name,"r");
  
  k = -1;
  while (1) {

    if (!fgets(buf,FILENAME_MAX-1,fp))
      break;
    
    // vieux format ou -p5 ?
    if (buf[0] == 'A' || buf[0] == 'D') { //old
      j = sscanf(buf+9,"%d",&k);
      j += sscanf(buf+38,"%lf",&par[0]);
      j += sscanf(buf+45,"%lf",&par[1]);
      j += sscanf(buf+52,"%lf",&par[2]);
      type = buf;
    }
    else
      {
	j = sscanf(buf,"%*d %*s %*s %d %*f %*f %*f %lf %lf %lf",&k, par,par+1,par+2);
	type = buf+2;
      }

    strength = par[0];

    // erreur: on a pas tout lu ou ca ne croit pas
      if ((j < 4) || ((type[0] == 'D') && (k < prevkD)) || ((k < prevkA)))      {
	// fprintf(stderr,"%d %c %d %d %d\n",j,type[0],k,prevkD,prevkA);
 	fprintf(stderr, "\nError in splice sites file %s, line %d\n", name, i);
 	exit(2);
      }
      
      if (type[0] == 'D' && strength != 0.0 && k>prevkD) {
	prevkD = k;
	vPosDonF.push_back( k-1 );
	vValDonF.push_back( strength );
      }  else
	if (strength != 0.0 && k > prevkA)  {
	  prevkA = k;
	  vPosAccF.push_back( k );
	  vValAccF.push_back( strength );
	}
      i++;
  }
  
  if (k == -1) fprintf(stderr,"WARNING: Empty splice predictor file !\n");
  fclose(fp);
}

// -------------------------------------
//  Read Splice Predictor reverse file.
// -------------------------------------
void SensorSPred :: ReadSPredR(char name[FILENAME_MAX+1], int SeqLen)
{
  FILE *fp;
  char buf[FILENAME_MAX];
  double strength;
  double par[3];
  char* type;
  int i = 1, j, k;
  int prevkA = INT_MAX, prevkD = INT_MAX;

  fp = FileOpen(NULL,name,"r");

  k = -1;
  while (1) {
    if (!fgets(buf,FILENAME_MAX-1,fp))
      break;

    // vieux format ou -p5 ?
    if (buf[0] == 'A' || buf[0] == 'D') { //old
      j = sscanf(buf+9,"%d",&k);
      j += sscanf(buf+38,"%lf",&par[0]);
      j += sscanf(buf+45,"%lf",&par[1]);
      j += sscanf(buf+52,"%lf",&par[2]);
      type = buf;
    }
    else
      {
	j = sscanf(buf,"%*d %*s %*s %d %*f %*f %*f %lf %lf %lf",&k, par,par+1,par+2);
	type = buf+2;
      }

    strength = par[0];

    // on ne lit pas tout on ca l'index position ne decroit pas
    if ((j < 4) || (type[0] == 'D' && k > prevkD) || (k > prevkA)) {
      fprintf(stderr, "\nError in splice sites file %s, line %d\n", name, i);
      exit(2);
    }
          
    if (strength != 0.0) 
      if (type[0] == 'D') {
	prevkD = k;
	vPosDonR.push_back( k );
	vValDonR.push_back( strength );
      }
      else {
	prevkA = k;
	vPosAccR.push_back( k-1 );
	vValAccR.push_back( strength );
      }
    i++;
  }
  
  if (k == -1) fprintf(stderr,"WARNING: Empty splice predictor file !\n");
  fclose(fp);
}


// ------------------------
//  GiveInfo signal SPred.
// ------------------------
void SensorSPred :: GiveInfo (DNASeq *X,int pos, DATA *d)
{
  bool update = false;
  double f;

  if ( (PositionGiveInfo == -1) || (pos != PositionGiveInfo+1) ) update = true; 
  PositionGiveInfo = pos;

  // Accepteur Forward
  if(!vPosAccF.empty()) {
    if (update) 
      iAccF = lower_bound(vPosAccF.begin(), vPosAccF.end(), pos)-vPosAccF.begin();
    
    if((iAccF<(int)vPosAccF.size()) && (vPosAccF[iAccF] == pos)) {
      f = pow(vValAccF[iAccF], accB) * accP;
      d->sig[DATA::Acc].weight[Signal::Forward] += log(f);
      d->sig[DATA::Acc].weight[Signal::ForwardNo] += log(1.0-f);
      iAccF++;
    }
  }
  
  // Accepteur Reverse
  if (!vPosAccR.empty()) {
    if (update) 
      iAccR = lower_bound(vPosAccR.begin(), vPosAccR.end(), pos)-vPosAccR.begin();

    if((iAccR<(int)vPosAccR.size()) && (vPosAccR[iAccR] == pos)) {
      f = pow(vValAccR[iAccR], accB) * accP;
      d->sig[DATA::Acc].weight[Signal::Reverse] += log(f);
      d->sig[DATA::Acc].weight[Signal::ReverseNo] += log(1.0-f);
      iAccR++;
    }
  }

  // Donneur Forward
  if(!vPosDonF.empty()) {
    if (update) 
      iDonF = lower_bound(vPosDonF.begin(), vPosDonF.end(), pos)-vPosDonF.begin();

    if ((iDonF<(int)vPosDonF.size()) && (vPosDonF[iDonF] == pos)) {
      f = pow(vValDonF[iDonF], donB) * donP;
      d->sig[DATA::Don].weight[Signal::Forward] += log(f);
      d->sig[DATA::Don].weight[Signal::ForwardNo] += log(1.0-f);
      iDonF++;
    }
  }
  
  // Donneur Reverse
  if(!vPosDonR.empty()) {
    if (update) 
      iDonR = lower_bound(vPosDonR.begin(), vPosDonR.end(), pos)-vPosDonR.begin();

    if((iDonR<(int)vPosDonR.size()) && (vPosDonR[iDonR] == pos)) {
      f = pow(vValDonR[iDonR], donB) * donP;
      d->sig[DATA::Don].weight[Signal::Reverse] += log(f);
      d->sig[DATA::Don].weight[Signal::ReverseNo] += log(1.0-f);
      iDonR++;
    }
  }

}


// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorSPred :: Plot(DNASeq *X)
{
  double f;

  for (int i =0; i < (int)vPosAccF.size(); i++) {
    f = pow(vValAccF[i], accB) * accP;
    PlotAcc(vPosAccF[i], 1, NORM(log(f),20.0));
  }
  
  for (int i =0; i < (int)vPosDonF.size(); i++) {
    f = pow(vValDonF[i], donB) * donP;
    PlotDon(vPosDonF[i], 1, NORM(log(f),20.0));
  }
  
  for (int i =0; i < (int)vPosAccR.size(); i++) {
    f = pow(vValAccR[i], accB) * accP;
    PlotAcc(vPosAccR[i], -1, NORM(log(f),20.0));
  }
 
  for (int i =0; i < (int)vPosDonR.size(); i++) {
    f = pow(vValDonR[i], donB) * donP;
    PlotDon(vPosDonR[i], -1, NORM(log(f),20.0));
  }
}

// ------------------
//  Post analyse
// ------------------
void SensorSPred :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}

// -------------------------------------
//  Read ReadNG2Gff3 GFF3 file.
// -------------------------------------
void SensorSPred :: ReadSPredGff3(char name[FILENAME_MAX+1], int SeqLen)
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
	  vPosAccF.push_back( tmpFeature->getLocus()->getEnd() );
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
	  vPosAccR.push_back( tmpFeature->getLocus()->getStart()-1);
	  vValAccR.push_back( tmpFeature->getScore() );
	}
      }
    }
    it++;
    i++;
  }
  delete geneFeatureSet;
}


