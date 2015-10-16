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
// $Id: Sensor.AnnotaStruct.cc,v 1.24 2010-02-01 16:18:41 sallet Exp $
// ------------------------------------------------------------------
// File:     Sensor.AnnotaStruct.cc
// Contents: Sensor AnnotaStruct
// ------------------------------------------------------------------

// *********************************************************************
// IMPORTANT WARNING : this plugin relies on pointeurs to char and float 
// for high level informations to be adequately and simply updated during
// Init() for parameter optimization. Only low level information should
// be specifically allocated. 
// *********************************************************************

#include "Sensor.AnnotaStruct.h"
#include <functional>
#include <stdio.h>
#include <stdlib.h>

// ************************************
// Comparators for Signal by Position
// ************************************
bool BySigPos(const Signals *A, const Signals *B)
{ return (A->pos < B->pos); };


struct BySigPosF : public std::binary_function<Signals *,Signals *,bool> {
    bool operator()(Signals * x, Signals * y) { return BySigPos(x, y); }
};

bool ByConStart(const Contents *A, const Contents *B)
{ return (A->start < B->start); };

bool ByConEnd(const Contents *A, const Contents *B)
{ return (A->end < B->end); };

//************************************************************
// Copy C string with associated malloc
//************************************************************
inline char * CharCopy(const char* source)
{
   return strcpy((char *)malloc(strlen(source)+1),source);
}

/*************************************************************
 **                     SensorAnnotaStruct                  **
 *************************************************************/

extern Parameters PAR;

// ----------------------
// Default constructor.
// ----------------------
SensorAnnotaStruct :: SensorAnnotaStruct (int n, DNASeq *X) : Sensor(n)
{
  char tempname[FILENAME_MAX+1];
  // all types are possible.
  type = Type_Any;

// Inilialize the inline parameters
  char exonRead[20], intronRead[20], cdsRead[20], npcRnaRead[20];
  char startRead[20] ,stopRead[20] ,accRead[20] ,donRead[20] ,tStartRead[20] ,tStopRead[20],  tStartNpcRead[20], tStopNpcRead[20] ;
  strcpy(exonRead,PAR.getC("AnnotaStruct.Exon*",        GetNumber()));
  strcpy(intronRead,PAR.getC("AnnotaStruct.Intron*",    GetNumber()));
  strcpy(cdsRead,PAR.getC("AnnotaStruct.CDS*",         GetNumber()));
  strcpy(npcRnaRead, PAR.getC("AnnotaStruct.npcRNA*",  GetNumber()));
  strcpy(startRead,PAR.getC("AnnotaStruct.Start*",     GetNumber()));
  strcpy(stopRead,PAR.getC("AnnotaStruct.Stop*",       GetNumber()));
  strcpy(accRead,PAR.getC("AnnotaStruct.Acc*",         GetNumber()));
  strcpy(donRead,PAR.getC("AnnotaStruct.Don*",         GetNumber()));
  strcpy(tStartRead,PAR.getC("AnnotaStruct.TrStart*",  GetNumber()));
  strcpy(tStopRead,PAR.getC("AnnotaStruct.TrStop*",    GetNumber()));
  strcpy(tStartNpcRead,PAR.getC("AnnotaStruct.TrStartNpc*", GetNumber()));
  strcpy(tStopNpcRead,PAR.getC("AnnotaStruct.TrStopNpc*",   GetNumber()));

  if (exonRead[0] != 'i') exonInline = 0;
  else	exonInline = 1;

  if (intronRead[0] != 'i') intronInline = 0;
  else intronInline = 1;

  if (cdsRead[0] != 'i') cdsInline = 0;
  else cdsInline = 1;

  if (npcRnaRead[0] != 'i') npcRnaInline = 0;
  else  npcRnaInline = 1;

  if (startRead[0] != 'i') startInline = 0;
  else startInline= 1;
  
  if (stopRead[0] != 'i') stopInline = 0;	
  else stopInline = 1;
  
  if (accRead[0] != 'i') accInline = 0;
  else accInline = 1;
  
  if (donRead[0] != 'i') donInline=0;
  else donInline= 1;
  
  if (tStartRead[0] != 'i') tStartInline = 0;
  else tStartInline = 1;

  if (tStopRead[0] != 'i') tStopInline = 0;
  else tStopInline = 1;

  if (tStartNpcRead[0] != 'i') tStartNpcInline = 0;
  else tStartNpcInline = 1;

  if (tStopNpcRead[0] != 'i') tStopNpcInline = 0;
  else tStopNpcInline = 1;


  fileExt       = PAR.getC("AnnotaStruct.FileExtension", GetNumber());
  transFeatName = to_string(PAR.getC("AnnotaStruct.TranscriptFeature", GetNumber()));
  inputFormat_  = to_string(PAR.getC("AnnotaStruct.format", GetNumber(),1));
  
  fprintf(stderr, "Reading %s file....", fileExt);
  fflush(stderr);
  strcpy(tempname, PAR.getC("fstname"));
  strcat(tempname, ".");
  strcat(tempname, fileExt);
  if ( inputFormat_ == "GFF3" )
  {
    strcat(tempname,".gff3");
    GeneFeatureSet * geneFeatureSet = new GeneFeatureSet (tempname);
    ReadAnnotaStructGff3(*geneFeatureSet,  X->SeqLen);
    delete geneFeatureSet;
  }
  else
  {
    ReadAnnotaStruct(tempname, X->SeqLen);
  }
  fprintf(stderr, "done\n");
  fflush(stderr);
  
  
  // Sort vSig by pos
  sort(vSig.begin(), vSig.end(), BySigPos);
  
  // Sort vCon by end and by start
  sort(vCon.begin(), vCon.end(), ByConEnd);
  stable_sort(vCon.begin(), vCon.end(), ByConStart);
  
  // Delete redundancy
  for(int i=0; i<(int)vSig.size(); i++) {
    if(i == 0) continue;
    if(vSig[i]->pos ==vSig[i-1]->pos  && vSig[i]->type ==vSig[i-1]->type &&
       vSig[i]->edge==vSig[i-1]->edge && vSig[i]->score==vSig[i-1]->score) {
      vSig.erase(i+vSig.begin());
      i--;
    }
  }

  for(int i=0; i<(int)vCon.size(); i++) {
    if(i == 0) continue;
    if(vCon[i]->start==vCon[i-1]->start && vCon[i]->end  ==vCon[i-1]->end &&
       vCon[i]->type ==vCon[i-1]->type  && vCon[i]->score==vCon[i-1]->score) {
      vCon.erase(i+vCon.begin());
      i--;
    }
  }

  // Patch because of non start/stop sites due to non full gene predicted
  // by fgeneshpasa.  If start stop is non canonical -> remove it
  for(int i=(int)vSig.size()-1; i>=0; i--) {
    if (vSig[i]->type == DATA::Start)
      if(vSig[i]->edge) {
	if (!X->IsEStart(vSig[i]->pos-1, -1)) { 
		vSig.erase(i+vSig.begin());
		fprintf(stderr,"Annotastruct: bad ATG on reverse strand\n");
	}
      }
      else {
	if (!X->IsEStart(vSig[i]->pos, 1))    { 
		vSig.erase(i+vSig.begin());
		fprintf(stderr,"Annotastruct: bad ATG on forward strand\n");
	}
      }
    if (vSig[i]->type ==  DATA::Stop)
      if(vSig[i]->edge) {
	if (!X->IsStop(vSig[i]->pos+2, -1)) { vSig.erase(i+vSig.begin()); }
      }
      else {
	if (!X->IsStop(vSig[i]->pos-3, 1))  { vSig.erase(i+vSig.begin()); }
      }
  }

  // FOR DEBUG
  std::vector <int> vPosAccF, vPosAccR;
  std::vector <int> vPosDonF, vPosDonR;
  std::vector <int> vPosStaF, vPosStaR;
  std::vector <int> vPosStoF, vPosStoR;
  for(int i=0; i<(int)vSig.size(); i++) {
    //vSig[i]->PrintS();
    switch (vSig[i]->type) {
    case DATA::Start : 
      if(vSig[i]->edge) vPosStaR.push_back(vSig[i]->pos);
      else vPosStaF.push_back(vSig[i]->pos);
      break;
    case DATA::Stop : 
      if(vSig[i]->edge) vPosStoR.push_back(vSig[i]->pos);
      else vPosStoF.push_back(vSig[i]->pos);
      break;
    case DATA::Acc :
      if(vSig[i]->edge) vPosAccR.push_back(vSig[i]->pos);
      else vPosAccF.push_back(vSig[i]->pos);
      break;
    case DATA::Don :
      if(vSig[i]->edge) vPosDonR.push_back(vSig[i]->pos);
    else vPosDonF.push_back(vSig[i]->pos);
      break;
    }
  }
  CheckStart  (X,vPosStaF, vPosStaR);
  CheckStop   (X,vPosStoF, vPosStoR);
  CheckSplices(X,vPosAccF, vPosDonF, vPosAccR, vPosDonR);
  vPosStaF.clear(); vPosStaR.clear();
  vPosStoF.clear(); vPosStoR.clear();
  vPosAccF.clear(); vPosDonF.clear(); vPosAccR.clear(); vPosDonR.clear();
  
  if (PAR.getI("Output.graph")) Plot(X);
}

// ----------------------
//  Default destructor.
// ----------------------
SensorAnnotaStruct :: ~SensorAnnotaStruct ()
{
  // Clear the data structures
  vSig.clear();
  vCon.clear();
}

// ----------------------
//  Init.
// ----------------------
void SensorAnnotaStruct :: Init (DNASeq *X)
{
  char exonRead[20], intronRead[20], cdsRead[20], npcRnaRead[20];
  char startRead[20], stopRead[20] ,accRead[20] ,donRead[20] ,tStartRead[20] ,tStopRead[20], tStartNpcRead[20], tStopNpcRead[20];

  strcpy(exonRead,   PAR.getC("AnnotaStruct.Exon*",   GetNumber()));
  strcpy(intronRead, PAR.getC("AnnotaStruct.Intron*", GetNumber()));
  strcpy(cdsRead,    PAR.getC("AnnotaStruct.CDS*",    GetNumber()));
  strcpy(npcRnaRead, PAR.getC("AnnotaStruct.npcRNA*", GetNumber()));

  if (exonRead[0] != 'i') exonPAR   = atof(exonRead);
  else exonPAR = 0;

  if (intronRead[0] != 'i') intronPAR = atof(intronRead);
  else intronPAR = 0;

  if (cdsRead[0] != 'i') cdsPAR = atof(cdsRead);
  else cdsPAR = 0;

  if (npcRnaRead[0] != 'i') npcRnaPAR = atof(npcRnaRead);
  else npcRnaPAR = 0;

  strcpy(startRead,PAR.getC("AnnotaStruct.Start*",     GetNumber()));
  strcpy(stopRead,PAR.getC("AnnotaStruct.Stop*",       GetNumber()));
  strcpy(accRead,PAR.getC("AnnotaStruct.Acc*",         GetNumber()));
  strcpy(donRead,PAR.getC("AnnotaStruct.Don*",         GetNumber()));
  strcpy(tStartRead,PAR.getC("AnnotaStruct.TrStart*",  GetNumber()));
  strcpy(tStopRead,PAR.getC("AnnotaStruct.TrStop*",    GetNumber()));
  strcpy(tStartNpcRead,PAR.getC("AnnotaStruct.TrStartNpc*", GetNumber()));
  strcpy(tStopNpcRead,PAR.getC("AnnotaStruct.TrStopNpc*",  GetNumber()));

  strcpy(startPAR,     PAR.getC("AnnotaStruct.StartType",   GetNumber()));
  strcpy(stopPAR,      PAR.getC("AnnotaStruct.StopType",    GetNumber()));
  strcpy(accPAR,       PAR.getC("AnnotaStruct.AccType",     GetNumber()));
  strcpy(donPAR,       PAR.getC("AnnotaStruct.DonType",     GetNumber()));
  strcpy(tStartPAR,    PAR.getC("AnnotaStruct.TrStartType", GetNumber()));
  strcpy(tStopPAR,     PAR.getC("AnnotaStruct.TrStopType",  GetNumber()));
  strcpy(tStartNpcPAR, PAR.getC("AnnotaStruct.TrStartNpcType", GetNumber()));
  strcpy(tStopNpcPAR,  PAR.getC("AnnotaStruct.TrStopNpcType",  GetNumber()));	

  if (startRead[0]     != 'i') strcat(startPAR,     startRead);
  if (stopRead[0]      != 'i') strcat(stopPAR,      stopRead);
  if (accRead[0]       != 'i') strcat(accPAR,       accRead);
  if (donRead[0]       != 'i') strcat(donPAR,       donRead);
  if (tStartRead[0]    != 'i') strcat(tStartPAR,    tStartRead);
  if (tStopRead[0]     != 'i') strcat(tStopPAR,     tStopRead);
  if (tStartNpcRead[0] != 'i') strcat(tStartNpcPAR, tStartNpcRead);
  if (tStopNpcRead[0]  != 'i') strcat(tStopNpcPAR,  tStopNpcRead);

  PosSigGiveInfo = -1;
  PosConGiveInfo = -1;
  iSig = 0;
  iCon = 0;
  
}

// --------------------------
//  Read start forward file.
// --------------------------
void SensorAnnotaStruct :: ReadAnnotaStruct(char name[FILENAME_MAX+1], int len)
{
  FILE  *fp;
  char  line[MAX_LINE];
  int   i;
  int   startC, endC, edge;   // C -> content
  int   startS, endS;         // S -> signal       need for reverse
  char  strand;
  char  phase[2];
  float scF;
  char  feature[20];
  char  scoreC[20];
  int   frame = -1;
 
  fp = FileOpen(NULL, name, "r", PAR.getI("EuGene.sloppy"));
  
  int j=0;
  while(fp  &&  fgets (line, MAX_LINE, fp) != NULL) {
    j++;
    if (line[0] != '#') {
      // GFF line : seqn source feature start end score strand phase
      i = sscanf(line, "%*s %*s %s %d %d %s %c %s",
		 feature, &startC, &endC, scoreC, &strand, phase);
      if (i < 6) {
	if (i==-1) {
	  if(j==1)
	    fprintf(stderr,"WARNING: empty AnnotaStruct file !...");
	}
	else {
	  fprintf(stderr, "Error in AnnotaStruct file %s, line %d.\n",name,j);
	  exit(2);
	}
      }
      else {
	
	/* Score ? */
	if (strcmp(scoreC, ".")) {
	  if (scoreC[0] == 'p' || scoreC[0] == 's') scF = atof(scoreC+1);
	  else                                      scF = atof(scoreC);
	}
	else scF = 0.0;

	/* Phase ? */
	if (strcmp(phase, ".")) {
	  if (strand == '+') { frame = (startC - 1)   % 3; }
	  else               { frame = (len   - endC) % 3; }
	  frame = (frame + atoi(phase)) % 3;
	}

	/* Strand ? */
	if      (strand == '+') {
	  edge = 0;
	  startS = startC;
	  endS   = endC;
	}
	else if (strand == '-') {
	  edge = 1;
	  startS = endC+1;
	  endS   = startC-1;
	}
	else {
	  fprintf(stderr, "WARNING: feature %s line %d strand unknown"
		  " => ignored.\n", feature, j);
	  continue;
	}
	startC--;
	endC--;
	
	/* Parse by feature */
	// Low level (contents OR signals)
	if     (!strcmp(feature, "trStart"))
	  vSig.push_back(new Signals(startS-1, DATA::tStart, edge, CharCopy(scoreC)));
	else if(!strcmp(feature, "trStop"))
	  vSig.push_back(new Signals(startS,   DATA::tStop,  edge, CharCopy(scoreC)));
	else if(!strcmp(feature, "start"))
	  vSig.push_back(new Signals(startS-1, DATA::Start,  edge, CharCopy(scoreC)));
	else if(!strcmp(feature, "stop"))
	  vSig.push_back(new Signals(startS,   DATA::Stop,   edge, CharCopy(scoreC)));
	else if(!strcmp(feature, "acc"))
	  vSig.push_back(new Signals(startS,   DATA::Acc,    edge, CharCopy(scoreC)));
	else if(!strcmp(feature, "don"))
	  vSig.push_back(new Signals(startS-1, DATA::Don,    edge, CharCopy(scoreC)));
	else if(!strcmp(feature, "ins"))
	  vSig.push_back(new Signals(startS,   DATA::Ins,    edge, CharCopy(scoreC)));
	else if(!strcmp(feature, "del"))
	  vSig.push_back(new Signals(startS,   DATA::Del,    edge, CharCopy(scoreC)));
	else if(!strcmp(feature, "trStartNpc"))
	  vSig.push_back(new Signals(startS-1, DATA::tStartNpc, edge, CharCopy(scoreC)));
        else if(!strcmp(feature, "trStopNpc"))
	  vSig.push_back(new Signals(startS, DATA::tStopNpc, edge, CharCopy(scoreC)));
	else if(!strcmp(feature, "exon"))
	  PushInCon(startC, endC, new float(scF), strand, phase, frame);
	else if(!strcmp(feature, "intron"))
	  vCon.push_back(new Contents(startC, endC, DATA::IntronF+edge, new float(scF)));
	else if(!strcmp(feature, "utr5"))
	  vCon.push_back(new Contents(startC, endC, DATA::UTR5F+edge,  new float(scF)));
	else if(!strcmp(feature, "interg"))
	  vCon.push_back(new Contents(startC, endC, DATA::InterG,  new float(scF)));	  
	else if(!strcmp(feature, "utr3"))
	  vCon.push_back(new Contents(startC, endC, DATA::UTR3F+edge,  new float(scF)));
	else if(!strcmp(feature, "utr")) {
	  vCon.push_back(new Contents(startC, endC, DATA::UTR5F+edge,  new float(scF)));
	  vCon.push_back(new Contents(startC, endC, DATA::UTR3F+edge,  new float(scF)));
	}
	else if (!strcmp(feature, "ncrna")) 
          vCon.push_back(new Contents(startC, endC, DATA::RNAF+edge,  new float(scF)));
	else if(!strcmp(feature, "intronutr"))
	  vCon.push_back(new Contents(startC,endC,DATA::IntronUTRF+edge,new float(scF)));
	
	// High level (contents OR/AND signals)
	else if(!strcmp(feature, "E.Init")) 
	{
	  vSig.push_back(new Signals(startS-1, DATA::Start, edge, startPAR));
	  vSig.push_back(new Signals(endS,     DATA::Don,   edge, donPAR));
	  PushInCon(startC, endC, &cdsPAR, strand, phase, frame);
	}
	else if(!strcmp(feature, "E.Intr")) {
	  vSig.push_back(new Signals(startS-1, DATA::Acc,   edge, accPAR));
	  vSig.push_back(new Signals(endS,     DATA::Don,   edge, donPAR));
	  PushInCon(startC, endC, &cdsPAR, strand, phase, frame);
	}
	else if(!strcmp(feature, "E.Term")) {
	  vSig.push_back(new Signals(startS-1, DATA::Acc,   edge, accPAR));
	  vSig.push_back(new Signals(endS,     DATA::Stop,  edge, stopPAR));
       	  PushInCon(startC, endC, &cdsPAR, strand, phase, frame);
	}
	else if(!strcmp(feature, "E.Sngl")) {
	  vSig.push_back(new Signals(startS-1, DATA::Start, edge, startPAR));
	  vSig.push_back(new Signals(endS,     DATA::Stop,  edge, stopPAR));
	  PushInCon(startC, endC, &cdsPAR, strand, phase, frame);
	}
	else if(!strcmp(feature, "UTR5")) {
	  vSig.push_back(new Signals (startS-1,DATA::tStart,edge, tStartPAR));
	  vCon.push_back(new Contents(startC,endC,DATA::UTR5F+edge, &cdsPAR));
	}
	else if(!strcmp(feature, "UTR3")) {
	  vSig.push_back(new Signals (endS,    DATA::tStop, edge, tStopPAR));
	  vCon.push_back(new Contents(startC,endC,DATA::UTR3F+edge, &cdsPAR));
	}
	else if(!strcmp(feature, "UTR")) {
	  vSig.push_back(new Signals(startS-1, DATA::tStart, edge,tStartPAR));
	  vSig.push_back(new Signals(endS,     DATA::tStop,  edge,tStopPAR));
	  vCon.push_back(new Contents(startC,endC,DATA::UTR5F+edge, &cdsPAR));
	  vCon.push_back(new Contents(startC,endC,DATA::UTR3F+edge, &cdsPAR));
	}
	else if(!strcmp(feature, "Intron")) {
	  vCon.push_back(new Contents(startC,endC,DATA::IntronF+edge,&intronPAR));
	}
	else if (strcmp(feature, "ncRNA")) {
	  vCon.push_back(new Contents(startC,endC,DATA::RNAF+edge,&npcRnaPAR));
	  vSig.push_back(new Signals(startS-1, DATA::tStartNpc, edge, tStartNpcPAR));
          vSig.push_back(new Signals(endS,     DATA::tStopNpc,  edge, tStopNpcPAR));
	}
	
	else if(!strcmp(feature, "E.Any")) {
	  PushInCon(startC, endC, &exonPAR, strand, phase, frame);
	  vCon.push_back(new Contents(startC,endC,DATA::UTR5F+edge, &exonPAR));
	  vCon.push_back(new Contents(startC,endC,DATA::UTR3F+edge, &exonPAR));
	}
	else if(!strcmp(feature, "Intron.Any")) {
	  vSig.push_back(new Signals (startS-1, DATA::Don, edge, donPAR));
	  vSig.push_back(new Signals (endS,     DATA::Acc, edge, accPAR));
	  vCon.push_back(new Contents(startC,endC,DATA::IntronF+edge,&intronPAR));
	  vCon.push_back(new Contents(startC,endC,DATA::IntronUTRF+edge,&intronPAR));
	}
	else if(!strcmp(feature, "E.First")) {
	  PushInCon(startC, endC, &exonPAR, strand, phase, frame);
	  vSig.push_back(new Signals (startS-1,DATA::tStart, edge, tStartPAR));
	  vCon.push_back(new Contents(startC,endC,DATA::UTR5F+edge, &exonPAR));
	}
	else if(!strcmp(feature, "E.Last")) {
	  PushInCon(startC, endC, &exonPAR, strand, phase, frame);
	  vSig.push_back(new Signals (endS,  DATA::tStop, edge, tStopPAR));
	  vCon.push_back(new Contents(startC,endC,DATA::UTR3F+edge, &exonPAR));
	}
	else if(!strcmp(feature, "E.Extreme")) {
	  PushInCon(startC, endC, &exonPAR, strand, phase, frame);
	  vSig.push_back(new Signals(startS-1, DATA::tStart, edge,tStartPAR));
	  vSig.push_back(new Signals(endS,     DATA::tStop,  edge,tStopPAR));
	  vCon.push_back(new Contents(startC,endC,DATA::UTR5F+edge, &exonPAR));
	  vCon.push_back(new Contents(startC,endC,DATA::UTR3F+edge, &exonPAR));
	}
	else
	  fprintf(stderr, "WARNING: feature %s line %d unknown => ignored.\n",
		  feature, j);
      }
    }
  }
  if (fp) fclose(fp);
}

//----------------------------------------------------------------------------------------------------------
// FillOntologyTerms : assigns SOTerms to GeneFeatures based on parent information to a transcript feature
//----------------------------------------------------------------------------------------------------------
void SensorAnnotaStruct ::FillOntologyTerm(GeneFeatureSet & geneFeatureSet)
{
   vector<GeneFeature *>::iterator it = geneFeatureSet.getIterator();
   int nbGeneFeature=geneFeatureSet.getNbFeature();
   int i=0;
   
   for ( i=0 ; i < nbGeneFeature ; i++, it++ )
   {
      // cout << "<" << (*it)->getType() << "> =? <" << transFeatName << ">\n";
      // Do we have a transcript feature ?
      if ((*it)->getType() == transFeatName)
      {
	 // cout << "YES\n";
         // Then scan the CDS and mark them with the correct SOTerm
	 // get the children of the "transcript" feature
         vector <GeneFeature *>& childrenArray = geneFeatureSet.getChildren((*it)->getId());
	 int c = childrenArray.size();
	 // cout << "size children: " << c << "\n";

	if (c == 0) continue; // No children... Ignore this.

         if (c == 1) // Sngl exon gene
	 {
	    // CDS with no SOTerm
	    if ((childrenArray[0]->getType() == "CDS") && 
	        (childrenArray[0]->getAttributes()->getOntologyTerm() == ""))
	    {
		childrenArray[0]->getAttributes()->setOntologyTerm("SO:0005845");
		// cout << "Setting SO:0005845 - Sngl to " << *(childrenArray[0]) << "\n";
	    }
            continue;
	 }

	 // Get the strandness
	 char strand = (*it)->getLocus()->getStrand();
	 int first,last,incr,lastCDSIdx = -1;
	 bool firstMet = false;

	 if (strand == '+')
	 {
	    first = 0;
	    last = c-1;
	    incr = 1;
	 } 
 	 else
         {
	    first = c-1;
	    last = 0;
	    incr = -1;
	 }

	 for (int j = first; j <= last; j += incr)
	 {
	    // cout << "Children " << j << " Type " << childrenArray[j]->getType() << " Ontol " << childrenArray[j]->getAttributes()->getOntologyTerm() << "\n";


	    // CDS with no SOTerm
	    if ((childrenArray[j]->getType() == "CDS") && 
	        (childrenArray[j]->getAttributes()->getOntologyTerm() == ""))
	    {
		if (!firstMet) // First CDS
		{
		    firstMet = true;
		    childrenArray[j]->getAttributes()->setOntologyTerm("SO:0000196");
		    // cout << "Setting SO:0000196 - First to " << *(childrenArray[j]) << "\n";
		    // TODO Check ATG !
		}
		else
		{
		    lastCDSIdx = j; // Internal Exon or may be LastExon
		    childrenArray[j]->getAttributes()->setOntologyTerm("SO:0000004");
		    // cout << "Setting SO:0000004 - Internal to " << *(childrenArray[j]) << "\n";
		    //  TODO Check Splices
		}
	    }
	 }
	 if (lastCDSIdx >= 0)
	 {
	    childrenArray[lastCDSIdx]->getAttributes()->setOntologyTerm("SO:0000197");
	    // cout << "Setting SO:0000197 - Last to " << *(childrenArray[lastCDSIdx]) << "\n";
	    // TODO Check Stop 
	 }
      }
   }
}


//gff3
void SensorAnnotaStruct ::ReadAnnotaStructGff3(GeneFeatureSet & geneFeatureSet , int len)
{
  int   startC, endC, edge;   // C -> content
  int   startS, endS;         // S -> signal       need for reverse
  char  strand;
  char  phase[2];
  float scF;
  char  feature[50];
  string  onthology_term;
  string idSo;  // id of the so term (ex: SO:0000100)

  char  scoreC[20];
  int   frame = -1;
  
  FillOntologyTerm(geneFeatureSet);
  int j=0;
  vector<GeneFeature *>::iterator it = geneFeatureSet.getIterator();
  
  int nbGeneFeature=geneFeatureSet.getNbFeature();
  int i=0;
  for ( i=0 ; i < nbGeneFeature ; i++, it++ )
  {
    j++;
    strcpy (scoreC, "s");
    strcpy (feature, (*it)->getType().c_str());
    startC         = (*it)->getLocus()->getStart();
    endC           = (*it)->getLocus()->getEnd();
    /* Score ? */
    scF            = (*it)->getScore();
    //score by default for ins/del;
    strcat ( scoreC, to_string(scF).c_str() ); 
    strand         = (*it)->getLocus()->getStrand();
    strcpy (phase, to_string((*it)->getPhase()).c_str() );
    onthology_term = (*it)->getAttributes()->getOntologyTerm();
    idSo           = (*it)->getType(); //recup code SO (third column of gff3)
    if ( idSo.find("SO:") == string::npos ) // If its not the SO code
    {
      string tmp = GeneFeatureSet::soTerms_->getIdFromName(idSo);
      idSo       = tmp;
    }

    /* Phase ? */
    if (strcmp(phase, "."))
    {
      if (strand == '+') { frame = (startC - 1) % 3; }
      else               { frame = (len - endC) % 3; }
      frame = (frame + atoi(phase)) % 3;
    }

    /* Strand ? */
    if (strand == '+') 
    {
      edge = 0;
      startS = startC;
      endS   = endC;
    }
    else if (strand == '-') 
    {
      edge = 1;
      startS = endC+1;
      endS   = startC-1;
    }
    else 
    {
      fprintf(stderr, "WARNING: feature %s line %d strand unknown"
	  " => ignored.\n", feature, j);
      continue;
    }
    startC--;
    endC--;
    
    // Low level signals : use inline values
    if      ( idSo =="SO:0000315" ) //trStart : transcription_start_site 
	vSig.push_back(new Signals(startS-1, DATA::tStart, edge, GetScoreC(DATA::tStart,scF,true)));

    else if ( idSo =="SO:0000616" ) //trStop : transcription_end_site
       	vSig.push_back(new Signals(startS,   DATA::tStop,  edge, GetScoreC(DATA::tStop,scF, true)));

    else if ( idSo =="SO:0000318" ) //start : start_codon
        vSig.push_back(new Signals(startS-1, DATA::Start,  edge, GetScoreC(DATA::Start,scF, true)));

    else if ( idSo =="SO:0000319" ) //stop : stop_codon
        vSig.push_back(new Signals(startS,   DATA::Stop,   edge, GetScoreC(DATA::Stop,scF, true)));
    
    else if ( idSo == "SO:0000164") //acc : three_prime_splice_site
        vSig.push_back(new Signals(startS,   DATA::Acc,    edge, GetScoreC(DATA::Acc,scF, true)));
    
    else if ( idSo == "SO:0000163") //don : five_prime_splice_site
        vSig.push_back(new Signals(startS-1, DATA::Don,    edge, GetScoreC(DATA::Don,scF, true)));

    else if ( idSo == "SO:0000366") //ins : insertion_site
      vSig.push_back(new Signals(startS,   DATA::Ins,    edge, CharCopy(scoreC)));

    else if ( idSo == "SO:0000687") //del : deletion_junction
      vSig.push_back(new Signals(startS,   DATA::Del,    edge, CharCopy(scoreC)));

    // High level (contents OR/AND signals)
    else if ( idSo == "SO:0000316" && onthology_term == "SO:0000196") 
      //CDS && five_prime_coding_exon_region == E.Init
    {
      vSig.push_back(new Signals(startS-1, DATA::Start, edge, GetScoreC(DATA::Start,scF, startInline)));
      vSig.push_back(new Signals(endS,     DATA::Don,   edge, GetScoreC(DATA::Don,  scF, donInline)));
      PushInCon(startC, endC, (cdsInline ? new float(scF) : &cdsPAR), strand, phase, frame);
    }

    else if ( idSo == "SO:0000316" && onthology_term == "SO:0000004") 
      //CDS && interior_coding_exon == E.Intr
    {
      vSig.push_back(new Signals(startS-1, DATA::Acc, edge,   GetScoreC(DATA::Acc,scF, accInline)));
      vSig.push_back(new Signals(endS,     DATA::Don,   edge, GetScoreC(DATA::Don,scF, donInline)));
      PushInCon(startC, endC, (cdsInline ? new float(scF) : &cdsPAR), strand, phase, frame);
    }

    else if ( idSo == "SO:0000316" && onthology_term == "SO:0000197") //CDS && three_prime_coding_exon_region == E.Term
    {
      vSig.push_back(new Signals(startS-1, DATA::Acc,   edge, GetScoreC(DATA::Acc,scF, accInline)));
      vSig.push_back(new Signals(endS,     DATA::Stop,  edge, GetScoreC(DATA::Stop,scF, stopInline)));
      PushInCon(startC, endC, (cdsInline ? new float(scF) : &cdsPAR), strand, phase, frame);
    }
    else if ( idSo == "SO:0000316" && onthology_term == "SO:0005845") //CDS && single_exon == "E.Sngl"
    {
      vSig.push_back(new Signals(startS-1, DATA::Start, edge, GetScoreC(DATA::Start,scF, startInline)));
      vSig.push_back(new Signals(endS,     DATA::Stop,  edge, GetScoreC(DATA::Stop,scF, stopInline)));
      PushInCon(startC, endC, (cdsInline ? new float(scF) : &cdsPAR), strand, phase, frame);
    }
    else if ( idSo == "SO:0000316") //CDS
      PushInCon(startC, endC, (cdsInline ? new float(scF) : &cdsPAR), strand, phase, frame);

    else if ( idSo == "SO:0000204" ) //five_prime_UTR
    {
      vSig.push_back(new Signals (startS-1,DATA::tStart,edge, GetScoreC(DATA::tStart,scF, tStartInline)));
      vCon.push_back(new Contents(startC,endC,DATA::UTR5F+edge, (cdsInline ? new float(scF) : &cdsPAR)));
    }

    else if ( idSo == "SO:0000205" ) //three_prime_UTR
    {
      vSig.push_back(new Signals (endS,    DATA::tStop, edge, GetScoreC(DATA::tStop,scF, tStopInline)));
      vCon.push_back(new Contents(startC,endC,DATA::UTR3F+edge, (cdsInline ? new float(scF) : &cdsPAR)));
    }
    else if ( idSo == "SO:0000203" ) //UTR
    {
      vSig.push_back(new Signals(startS-1, DATA::tStart, edge, GetScoreC(DATA::tStart,scF, tStartInline)));
      vSig.push_back(new Signals(endS,     DATA::tStop,  edge, GetScoreC(DATA::tStop,scF, tStopInline)));

      vCon.push_back(new Contents(startC,endC,DATA::UTR5F+edge, (cdsInline ? new float(scF) : &cdsPAR)));
      vCon.push_back(new Contents(startC,endC,DATA::UTR3F+edge, (cdsInline ? new float(scF) : &cdsPAR)));
    }

    else if ( idSo == "SO:0000188" && onthology_term == "SO:0000191") //intron not UTR !
      vCon.push_back(new Contents(startC,endC,DATA::IntronF+edge, (intronInline ? new float(scF) : &intronPAR)));

    else if ( idSo == "SO:0000188" ) // "Intron.Any"
    {
      vSig.push_back(new Signals (startS-1, DATA::Don, edge, GetScoreC(DATA::Don,scF, donInline)));
      vSig.push_back(new Signals (endS,     DATA::Acc, edge, GetScoreC(DATA::Acc,scF, accInline)));
      vCon.push_back(new Contents(startC,endC,DATA::IntronF+edge, (intronInline ? new float(scF) : &intronPAR)));
      vCon.push_back(new Contents(startC,endC,DATA::IntronUTRF+edge, (intronInline ? new float(scF) : &intronPAR)));
    }

    else if ( idSo == "SO:0000147" && onthology_term == "SO:0000200") //"E.First"
    {
      PushInCon(startC, endC, (exonInline ? new float(scF) : &exonPAR), strand, phase, frame);
      vCon.push_back(new Contents(startC,endC,DATA::UTR5F+edge, (exonInline ? new float(scF) : &exonPAR)));
      vSig.push_back(new Signals (startS-1,DATA::tStart, edge, GetScoreC(DATA::tStart,scF, tStartInline)));
    }

    else if ( idSo == "SO:0000147" && onthology_term == "SO:0000202") // "E.Last"
    {
      PushInCon(startC, endC, (exonInline ? new float(scF) : &exonPAR), strand, phase, frame);
      vCon.push_back(new Contents(startC,endC,DATA::UTR3F+edge, (exonInline ? new float(scF) : &exonPAR)));
      vSig.push_back(new Signals (endS,  DATA::tStop, edge, GetScoreC(DATA::tStop,scF, tStopInline)));
    }
    else if (idSo == "SO:0000147" && onthology_term.find("SO:0000202")!= string::npos && onthology_term.find("SO:0000200")!= string::npos) 
    {
      PushInCon(startC, endC, (exonInline ? new float(scF) : &exonPAR), strand, phase, frame);
      vCon.push_back(new Contents(startC,endC,DATA::UTR5F+edge, (exonInline ? new float(scF) : &exonPAR)));
      vCon.push_back(new Contents(startC,endC,DATA::UTR3F+edge, (exonInline ? new float(scF) : &exonPAR)));
      vSig.push_back(new Signals(startS-1, DATA::tStart, edge, GetScoreC(DATA::tStart,scF, tStartInline)));
      vSig.push_back(new Signals(endS,     DATA::tStop,  edge, GetScoreC(DATA::tStop,scF, tStopInline)));
    }
    else if ( idSo == "SO:0000147" ) //E.Any 
    {
      PushInCon(startC, endC, (exonInline ? new float(scF) : &exonPAR), strand, phase, frame);
      vCon.push_back(new Contents(startC,endC,DATA::UTR5F+edge, (exonInline ? new float(scF) : &exonPAR)));
      vCon.push_back(new Contents(startC,endC,DATA::UTR3F+edge, (exonInline ? new float(scF) : &exonPAR)));
    }
    else if (GeneFeatureSet::soTerms_->isANcRNA(idSo)) // ncRNA or a kind of ncRNA (tRNA, rRNA, ...)
    {
      vCon.push_back( new Contents(startC, endC, DATA::RNAF+edge, (npcRnaInline ? new float(scF) : &npcRnaPAR)));
      vSig.push_back( new Signals (startS-1, DATA::tStartNpc, edge, GetScoreC(DATA::tStartNpc, scF, tStartNpcInline)));
      vSig.push_back( new Signals (endS,     DATA::tStopNpc,  edge, GetScoreC(DATA::tStopNpc, scF, tStopNpcInline))  );
    }
    else if ((*it)->getType() != transFeatName)
      fprintf(stderr, "WARNING: feature %s line %d unknown => ignored.\n",
	      feature, j);
    
    //fprintf(stderr, "END : feature %s line %d idSO : %s, Ontomlogy_term: %s.\n",
	   // feature,j, idSo.c_str(), onthology_term.c_str() );
    
  }
   
}

char *SensorAnnotaStruct :: GetScoreC (int type , float scF, bool inlineScore)
{
    char * scoreC = (inlineScore ? new char[20] : NULL);

	if (type == DATA::tStart)
		if (!inlineScore) return tStartPAR;
		else
		{
	    		strcpy (scoreC ,tStartPAR );
	 		strcat (scoreC , to_string(scF).c_str()); // ajout du score lu
			return scoreC;
		}


	else if (type == DATA::tStop)
		if (!inlineScore) return tStopPAR;
		else
		{
			strcpy (scoreC ,tStopPAR ); 
			strcat (scoreC , to_string(scF).c_str()); // ajout du score lu
			return scoreC;
		}

	else if (type == DATA::Start)
		if (!inlineScore) return startPAR;
		else
		{
			strcpy (scoreC ,startPAR );
			strcat (scoreC , to_string(scF).c_str());
			return scoreC;
		}

	else if (type == DATA::Stop)
		if (!inlineScore) return stopPAR;
		else
		{
			strcpy (scoreC ,stopPAR );
			strcat (scoreC , to_string(scF).c_str());
			return scoreC;
		}

	else if (type == DATA::Acc)
		if (!inlineScore) return accPAR;
		else
		{
			strcpy (scoreC ,accPAR );
			strcat (scoreC , to_string(scF).c_str());
			return scoreC;
		}

	else if (type == DATA::Don)
		if (!inlineScore) return donPAR;
		else
		{
			strcpy (scoreC ,donPAR );
			strcat (scoreC , to_string(scF).c_str()); 
			return scoreC;
		}
	else if (type == DATA::tStartNpc)
		if (!inlineScore) return tStartNpcPAR;
		else
		{
	    		strcpy (scoreC ,tStartNpcPAR );
	 		strcat (scoreC , to_string(scF).c_str()); // ajout du score lu
			return scoreC;
		}
	else if (type == DATA::tStopNpc)
		if (!inlineScore) return tStopNpcPAR;
		else
		{
	    		strcpy (scoreC ,tStopNpcPAR );
	 		strcat (scoreC , to_string(scF).c_str()); // ajout du score lu
			return scoreC;
		}
}
// ----------------
//  push_back con.
// ----------------
void SensorAnnotaStruct :: PushInCon(int d, int e, float *sc,
				     char st, char p[2], int f)
{
  int k;
  if (st == '-') k = 3;
  else           k = 0;
  if (! strcmp(p, ".") ) {
    vCon.push_back(new Contents(d, e, DATA::ExonF1+k, sc));
    vCon.push_back(new Contents(d, e, DATA::ExonF2+k, sc));
    vCon.push_back(new Contents(d, e, DATA::ExonF3+k, sc));
  }
  else
    vCon.push_back(new Contents(d, e, f+k, sc));
}

// ------------------------
//  GiveInfo.
// ------------------------
void SensorAnnotaStruct :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  Signals TmpSig;
  bool  update = false;
  int   i, j, iConTMP;
  float k;
  
  /* Signals */
  // update indexes on vector ?
  if((PosSigGiveInfo == -1) || (pos != PosSigGiveInfo+1)) update = true;
  PosSigGiveInfo = pos;
  if (!vSig.empty()) {
    if(update) {
      TmpSig.pos = pos;
      iSig = lower_bound(vSig.begin(), vSig.end(),
			 &TmpSig, BySigPosF()) - vSig.begin();
    }
    while((iSig<(int)vSig.size()) && (vSig[iSig]->pos == pos)) {
      i = vSig[iSig]->type;
      j = vSig[iSig]->edge;
      if (vSig[iSig]->score[0] == 'p'  ||  vSig[iSig]->score[0] == 's')
	k = atof(vSig[iSig]->score + 1);
      else
	k = atof(vSig[iSig]->score);
      if (vSig[iSig]->score[0] == 'p') {     // Probability [0 1]
	d->sig[i].weight[j]   += log(k);
	d->sig[i].weight[j+2] += log(1-k);
      }
      else                                   // Score ]-oo +oo[
	d->sig[i].weight[j] += k;
      iSig++;
    }
  }

  /* Contents */
  // update indexes on vector ?
  if((PosConGiveInfo == -1) || (pos != PosConGiveInfo+1)) update = true;
  PosConGiveInfo = pos;
  
  if(!vCon.empty()) {
    if(update) {
      iCon = 0;
      while(iCon < (int)vCon.size()  &&  pos > vCon[iCon]->end) iCon++;
    }
    iConTMP = iCon;
    while((iConTMP < (int)vCon.size())   &&
	  (pos >= vCon[iConTMP]->start)  &&  (pos <= vCon[iConTMP]->end)) {
      d->contents[vCon[iConTMP]->type] += *(vCon[iConTMP]->score);
      iConTMP++;
    }
    while(iCon < (int)vCon.size()  &&  pos > vCon[iCon]->end) iCon++;
  }
}

// ----------------------------
//  Plot Sensor information.
// ----------------------------
void SensorAnnotaStruct :: Plot(DNASeq *X)
{
}

// ------------------
//  Post analyse.
// ------------------
void SensorAnnotaStruct :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}
