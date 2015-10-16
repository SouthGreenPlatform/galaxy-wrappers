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
// $Id: Sensor.Tester.cc,v 1.25 2009-08-28 12:48:46 sallet Exp $
// ------------------------------------------------------------------
// File:     Sensor.Tester.cc
// Contents: Sensor Tester
// ------------------------------------------------------------------

// MacOS-X kludge. cmath undefines these macros. Turn them into inlines 
#include <math.h>
inline int (isinf)(double r) { return isinf(r); }
inline int (isnan)(double r) { return isnan(r); }

#include <iomanip>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "Sensor.Tester.h"

#include "../../MSensor.h"
#include "../../SensorIF.h"

extern Parameters PAR;
extern MasterSensor* MS;


/*************************************************************
 **                      SensorTester                       **
 *************************************************************/
// ----------------------
//  Default constructor.
// ----------------------
SensorTester :: SensorTester (int n, DNASeq *X) : Sensor(n)
{
  Todo = PAR.getC("Tester.Make");

  type = Type_None;
  inputFormat_ = to_string(PAR.getC("Tester.format", GetNumber(),1));

  if (Todo == "SPSN") {
    DATA Data;
    double dF, dR, dFdon, dRdon;
    int NbNonCanonicalDetected, NbNonCanonicalAnnotated; 
    bool is_posF, is_posR, is_posFdon, is_posRdon; 
    bool is_annotF, is_annotR, is_annotFdon, is_annotRdon;

    IsSPSN = true;

    if (!IsInitialized) {
      SensorName = "Sensor." + (std::string) PAR.getC("Tester.Sensor"); 
      SensorInstance = PAR.getI("Tester.Sensor.Instance"); 
      MinNumbers = PAR.getI("Tester.SPSN.MinNumbers"); 
      std::string evaluated_name = (std::string) PAR.getC("Tester.SPSN.Eval");
      if (evaluated_name == "START")  EvaluatedType = EVALUATED_START; 
      else if (evaluated_name == "STOP") EvaluatedType = EVALUATED_STOP;
      else if (evaluated_name == "SPLICE") EvaluatedType = EVALUATED_SPLICE;
      else {std::cerr<<"ERROR: The value "<<evaluated_name<<" set to the Tester.SPSN.Eval parameter is not accepted.\n"; exit(2);}

      // initialize once some class (static) attributs
      NbToDoSequence = (int) PAR.getD("NbSequence");
      NbDoneSequence = 0;
      for (int k=0; k< NbToDoSequence; k++) {
      	std::vector< std::vector<double> > v1; Scores.push_back( v1 );
      	std::vector< std::vector<bool> > v2;   IsAnnotateds.push_back( v2 );
      	std::vector< std::vector<bool> > v3;   IsAPositions.push_back( v3 );
      }

      IsInitialized = true;
    }

    sensor = MS->MakeSensor( SensorName, SensorInstance, X); sensor->Init(X);
    //gene = ReadGFFAnnotation(X);
    ReadAnnotation(X);
    SensorType = sensor->type;

    NbNonCanonicalDetected = 0;
    NbNonCanonicalAnnotated = 0;
    for (int i=0; i<X->SeqLen; i++) {
      for(int j=0; j<DATA::LastSigType; j++) Data.sig[j].Clear();
      sensor->GiveInfo(X,i,&Data);

      if (EvaluatedType == EVALUATED_START) {
	dF = Data.sig[DATA::Start].weight[Signal::Forward];
	dR = Data.sig[DATA::Start].weight[Signal::Reverse];
	is_posF = ( X->IsEStart(i,1) != 0 );
	is_posR = ( X->IsEStart(i-1,-1) != 0 );
	is_annotF = gene->IsState(DATA::Start,i,'+');
	is_annotR = gene->IsState(DATA::Start,i,'-');
	if (is_posF || is_posR) {
	  UpdateThreshold(dF); UpdateThreshold(dR); 
	  std::vector<double> v4; v4.push_back(dF); v4.push_back(dR);
	  Scores[NbDoneSequence].push_back( v4 );
	  std::vector<bool> v5; v5.push_back( is_annotF ); v5.push_back( is_annotR );
	  IsAnnotateds[NbDoneSequence].push_back( v5 );
	  std::vector<bool> v6; v6.push_back( is_posF ); v6.push_back( is_posR );
	  IsAPositions[NbDoneSequence].push_back( v6 );
	} else {
	  if (dF!=0) NbNonCanonicalDetected += 1; if (dR!=0) NbNonCanonicalDetected += 1;
	  if (is_annotF) NbNonCanonicalAnnotated += 1; if (is_annotR) NbNonCanonicalAnnotated += 1;
	}

      } else if (EvaluatedType == EVALUATED_STOP) {
	dF = Data.sig[DATA::Stop].weight[Signal::Forward];
	dR = Data.sig[DATA::Stop].weight[Signal::Reverse];
	is_posF = ( X->IsStop(i-3,1) != 0 );
	is_posR = ( X->IsStop(i+2,-1) != 0 );
	is_annotF = gene->IsState(DATA::Stop,i,'+');
	is_annotR = gene->IsState(DATA::Stop,i,'-');
	if (is_posF || is_posR) {
	  UpdateThreshold(dF); UpdateThreshold(dR); 
	  std::vector<double> v4; v4.push_back(dF); v4.push_back(dR);
	  Scores[NbDoneSequence].push_back( v4 );
	  std::vector<bool> v5; v5.push_back( is_annotF ); v5.push_back( is_annotR );
	  IsAnnotateds[NbDoneSequence].push_back( v5 );
	  std::vector<bool> v6; v6.push_back( is_posF ); v6.push_back( is_posR ); 
	  IsAPositions[NbDoneSequence].push_back( v6 );
	} else {
	  if (dF!=0) NbNonCanonicalDetected += 1; if (dR!=0) NbNonCanonicalDetected += 1;
	  if (is_annotF) NbNonCanonicalAnnotated += 1; if (is_annotR) NbNonCanonicalAnnotated += 1;
	}

      } else if (EvaluatedType == EVALUATED_SPLICE) {
	dF = Data.sig[DATA::Acc].weight[Signal::Forward];
	dR = Data.sig[DATA::Acc].weight[Signal::Reverse];
	dFdon = Data.sig[DATA::Don].weight[Signal::Forward];
	dRdon = Data.sig[DATA::Don].weight[Signal::Reverse];
	is_posF = ( ((*X)[i-2]=='a' && (*X)[i-1]=='g') );
	is_posR = ( ((*X)(i)=='g' && (*X)(i+1)=='a') );
	is_posFdon = ( ((*X)[i]=='g' && ( (*X)[i+1]=='t' || (*X)[i+1]=='c' )) );
	is_posRdon = ( ((*X)(i-1)=='g' && ( (*X)(i-2)=='t' || (*X)(i-2)=='c' )) );
	is_annotF = gene->IsState(DATA::Acc,i,'+');
	is_annotR = gene->IsState(DATA::Acc,i,'-');
	is_annotFdon = gene->IsState(DATA::Don,i,'+');
	is_annotRdon = gene->IsState(DATA::Don,i,'-');
	if (is_annotF && !is_posF) NbNonCanonicalAnnotated += 1;
	else if (is_annotR && !is_posR) NbNonCanonicalAnnotated += 1;
	else if (is_annotFdon && !is_posFdon) NbNonCanonicalAnnotated += 1;
	else if (is_annotRdon && !is_posRdon) NbNonCanonicalAnnotated += 1;
	else if (dF!=0 && !is_posF) NbNonCanonicalDetected += 1;
	else if (dR!=0 && !is_posR) NbNonCanonicalDetected += 1;
	else if (dFdon!=0 && !is_posFdon) NbNonCanonicalDetected += 1;
	else if (dRdon!=0 && !is_posRdon) NbNonCanonicalDetected += 1;
	else if (is_posF || is_posR || is_posFdon || is_posRdon) {
	  UpdateThreshold(dF); UpdateThreshold(dR); 
	  UpdateThreshold(dFdon); UpdateThreshold(dRdon); 
	  std::vector<double> v4; v4.push_back(dF); v4.push_back(dR); 
	                          v4.push_back(dFdon); v4.push_back(dRdon);
	  Scores[NbDoneSequence].push_back( v4 );
	  std::vector<bool> v5; v5.push_back( is_annotF ); v5.push_back( is_annotR ); 
                                v5.push_back( is_annotFdon ); v5.push_back( is_annotRdon );
	  IsAnnotateds[NbDoneSequence].push_back( v5 );
	  std::vector<bool> v6; v6.push_back( is_posF ); v6.push_back( is_posR );
	                        v6.push_back( is_posFdon ); v6.push_back( is_posRdon ); 
	  IsAPositions[NbDoneSequence].push_back( v6 );
	}
      } else 
	{std::cerr<<"ERROR: bad type of sensor to analyze with Tester in SensorTester::SensorTester.\n"; exit(2);}
    }
    
    if (NbNonCanonicalDetected != 0)
      std::cout<<"BE CAREFUL: in the sequence "<<X->Name<<", "<<NbNonCanonicalDetected
	       <<" SITES DETECTED BY THE SENSOR ARE NON CANONICAL\n"
	       <<"AND NOT TAKEN INTO ACCOUNT TO ESTIMATE SENSIBILITY AND SPECIFICITY.\n";
    if (NbNonCanonicalAnnotated != 0)
      std::cout<<"BE CAREFUL: in the sequence "<<X->Name<<", "<<NbNonCanonicalAnnotated
	       <<" SITES ARE ANNOTATED.\n"
	       <<"AND NOT TAKEN INTO ACCOUNT TO ESTIMATE SENSIBILITY AND SPECIFICITY.\n";

    NbDoneSequence+=1; 
    // the last sequence, compute specificity and sensibility
    if (NbDoneSequence == NbToDoSequence)  AnalyzeSPSN();

  } else if (Todo == "TEST") {
    IsSPSN = false;

    if (!IsInitialized) {
      std::string OutputFile;

      SensorName = (std::string) PAR.getC("Tester.Sensor"); 
      OutputFile = (std::string)PAR.getC("Output.Prefix")+"test."+SensorName+".gff";
      SensorInstance = PAR.getI("Tester.Sensor.Instance"); 
  
      fp = new FILE;
      // On verif outputFile n'existe pas
      if (!(fp = fopen(OutputFile.c_str(), "r"))) {
	if (!(fp = fopen(OutputFile.c_str(), "w"))) 
	  {fprintf(stderr, "cannot open %s output file\n", OutputFile.c_str());exit(2);}
	else
	  fprintf(fp,"SeqName\t Source\tFeature\t  Start\t    End\t  Score"
		  "\t Strand\t  Frame\t    T/F\t  State\n");
      }
      else {
	fclose(fp);
	fprintf(stderr,"WARNING: test output file \"%s\" exist\n", OutputFile.c_str());exit(2);
      }
      IsInitialized = true;
    }
    ReadAnnotation(X);
    sensor = MS->MakeSensor( "Sensor."+SensorName, SensorInstance, X );
    sensor->Init(X);

    // Sequence name without extension
    strcpy(seqName, BaseName(PAR.getC("fstname")));
    if (char * suffix = rindex(seqName,'.')) *suffix = 0;

  } else
    {std::cerr<<"ERROR: bad value given at the Tester.Make parameter.\n"; exit(2);}
}


// --------------------
//  Default destructor : free memory
// --------------------
SensorTester :: ~SensorTester (void)
{
  if (!IsSPSN) {
    fclose(fp);
    delete sensor;
    delete gene;
  }
}


// -------------------
//  Init from new seq
// -------------------
void SensorTester :: Init (DNASeq *X)
{
}


// ----------
//  GiveInfo
// ----------
void SensorTester :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
  if (!IsSPSN) {    
    DATA Data;
    const char *predSigType = new char[10];
    const char *tf;
    const char *truthState;
    for(int j=0; j<DATA::LastSigType; j++) Data.sig[j].Clear();
    sensor->GiveInfo(X,pos,&Data);

    for(int j=DATA::Start; j<=DATA::Don; j++) {
      // Forward
      if (Data.sig[j].weight[Signal::Forward] != 0.0) {
	tf = SigType_TF(j, pos, &predSigType);
	
	truthState = State(pos);
	if (!strcmp(truthState, "ExonR") || !strcmp(truthState, "IntronR"))
	  tf = "False";
	
	fprintf(fp,"%s\t%7.7s\t%7s\t%7d\t      ."
		"\t%7.2f\t      +\t      .\t%7s\t%7s\n",
		seqName, SensorName.c_str(), predSigType, pos,
		Data.sig[j].weight[Signal::Forward],
		tf, truthState);
      }
      // Reverse
      if(Data.sig[j].weight[Signal::Reverse] != 0.0) {
	tf = SigType_TF(j, pos, &predSigType);
	
	truthState = State(pos);
	if (!strcmp(truthState, "ExonF") || !strcmp(truthState, "IntronF"))
	  tf = "False";
	
	fprintf(fp,"%s\t%7.7s\t%7s\t%7d\t      ."
		"\t%7.2f\t      -\t      .\t%7s\t%7s\n",
		seqName, SensorName.c_str(), predSigType, pos,
		Data.sig[j].weight[Signal::Reverse],
		tf, truthState);
      }
    }
    delete [] predSigType;  //DOES NOT WORK, WHY ??????
  } 
}


// -------------------------
//  Plot Sensor information
// -------------------------
void SensorTester :: Plot(DNASeq *X)
{
}


// --------------
//  Post analyse
// --------------
void SensorTester :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
}


void SensorTester ::ReadAnnotation ( DNASeq *X)
{
    std::cerr <<"Reading coordinates file......................"; fflush(stderr);
    char tempname[FILENAME_MAX+1];

    strcpy(tempname,PAR.getC("fstname"));
    strcat(tempname,".gff");

    if ( inputFormat_ == "GFF3" )
    {
	strcat(tempname,".gff3");
	
	GeneFeatureSet * geneFeatureSet = new GeneFeatureSet (tempname);
	ReadGFF3Annotation(*geneFeatureSet, X);
	delete geneFeatureSet;
    }
    else
    {
	ReadGFFAnnotation(tempname , X);
    }
}


// -----------------
//  Read coord file
//  Lecture du fichier .gff (format gff) pour charger un objet
//  prediction "simplifi�" (on distingue 6 �tats : IG, UTR, ExonF,
//  ExonR, IntronF et IntronR) -> but quel �tat r�el � telle pos ?
//  Mots clefs pour le .gff : Features UTR5, UTR3, E.Init, E.Intr,
//  E.Term, E.Sngl
// -----------------
// BE CAREFUL:
//         GFF with one complete gene, that could not specify the UTR. 
//         In this case, the state from 0 to the first exon in the GFF 
//         is in InterGen and the state after the last exon in the GFF  
//         is not set Introns are all in the state IntronF1.
// -----------------
void SensorTester :: ReadGFFAnnotation(char name[FILENAME_MAX+1], DNASeq *x)
{
  
  FILE *fpCoord;
  char line[MAX_LINE];
  int  i;
  char *feature = new char[FILENAME_MAX];
  int  start,  end;
  char strand, frame;

  std::vector <int>         vPos;
  std::vector <signed char> vState;

  
  if (!(fpCoord = fopen(name, "r"))) 
    {std::cerr<<"Cannot open gff file " << name <<"\n"; exit(2);}

  int j=0;
  while(fgets (line, MAX_LINE, fpCoord) != NULL) {
    if (line[0] != '#') {
      j++;
      i = sscanf(line,"%*s %*s %s %d %d %*s %c %c",
		 feature, &start, &end, &strand, &frame);
      if (i < 5) {
	if (i==-1) {
	  if(j==1) std::cerr<<"WARNING: empty gff file !...";
	} else 
	  {std::cerr<<"\nError in gff file "<<name<<" line "<<j<<".\n";exit(2);}
      } else if (strcmp(feature,"Intron") != 0) {
	if (j==1) {
	  vPos.push_back  ( start-1  );
	  vState.push_back( InterGen );
	  if (strcmp(feature,"UTR5") == 0 || strcmp(feature,"UTR3") == 0) {
	    if (Todo=="TEST") {
	      vPos.push_back  ( end   );
	      vState.push_back( UTR5F );
	    }
	  }
	  else if (strcmp(feature, "E.Init") == 0) {
	    vPos.push_back  ( end    );
	    vState.push_back( InitF1 );
	  }
	  else if (strcmp(feature, "E.Term") == 0) {
	    vPos.push_back  ( end    );
	    vState.push_back( InitR1 );
	  }
	  else if (strcmp(feature, "E.Sngl") == 0)
	    if (strand == '+') {
	      vPos.push_back  ( end    );
	      vState.push_back( InitF1 );
	    }
	    else {
	      vPos.push_back  ( end    );
	      vState.push_back( InitR1 );
	    }
	  else {
	    std::cerr <<"\n Error in gff file "<<name
		      <<" line "<<j<<".\n"
		      <<" WARNING :\n"
		      <<"   - Complete genes only in gff file.\n"
		      <<"   - Feature must be UTR5, UTR3, E.Init,"
		      <<" E.Intr, E.Term or E.Sngl.\n";
	    exit(2);
	  }
	}
	else {
	  if (strcmp(feature,"UTR5") == 0 || strcmp(feature,"UTR3") == 0) {
	    if (Todo=="TEST") {
	      vPos.push_back  ( end   );
	      vState.push_back( UTR5F );
	    }
	  }
	  else if (strcmp(feature,"E.Init") == 0)
	    if (strand == '+') {
	      vPos.push_back  ( end    );
	      vState.push_back( InitF1 );
	    }
	    else {
	      vPos.push_back  ( start-1  );
	      vState.push_back( IntronR1 );
	      vPos.push_back  ( end      );
	      vState.push_back( InitR1   );
	    }
	  else if (strcmp(feature,"E.Term") == 0)
	    if (strand == '-') {
	      vPos.push_back  ( end    );
	      vState.push_back( InitR1 );
	    }
	    else {
	      vPos.push_back  ( start-1  );
	      vState.push_back( IntronF1 );
	      vPos.push_back  ( end      );
	      vState.push_back( InitF1   );
	    }
	  else if (strcmp(feature,"E.Sngl") == 0)
	    if (strand == '+') {
	      vPos.push_back  ( end      );
	      vState.push_back( InitF1   );
	    }
	    else {
	      vPos.push_back  ( end    );
	      vState.push_back( InitR1 );
	    }
	  else if (strcmp(feature,"E.Intr") == 0) {
	    if (strand == '+') {
	      vPos.push_back  ( start-1  );
	      vState.push_back( IntronF1 );
	      vPos.push_back  ( end      );
	      vState.push_back( IntrF1   );
	    }
	    else {
	      vPos.push_back  ( start-1  );
	      vState.push_back( IntronR1 );
	      vPos.push_back  ( end      );
	      vState.push_back( InitR1   );
	    }
	  }
	  else {
	    std::cerr <<"\n Error in gff file "<<name
		      <<" line "<<j<<".\n"
		      <<" "<<feature<<" : unknown feature (UTR5, UTR3, E.Init,"
		      <<" E.Intr, E.Term or E.Sngl).\n";
	    exit(2);
	  }
	}
      }
    }
  }
  fclose(fpCoord);

  // Complete gene ?
  if ((strcmp(feature, "E.Intr")  == 0  ||
      (strcmp(feature, "E.Init") == 0  &&  strand == '+')  ||
      (strcmp(feature, "E.Term") == 0  &&  strand == '-') ) && !(strcmp(feature,"UTR5") == 0 || strcmp(feature,"UTR3") == 0))
    {
      std::cerr <<"\n Error in gff file "<< name <<" line "<<j<<".\n";
      std::cerr <<" WARNING : complete genes only in gff file.\n";
      exit(2);
    }
 
  std::cerr <<"done\n";

  delete [] feature;

  gene = new Prediction(0,x->SeqLen,vPos, vState);
  gene->TrimAndUpdate(x);

  vPos.clear();
  vState.clear();
}

void SensorTester :: ReadGFF3Annotation(GeneFeatureSet & geneFeatureSet , DNASeq *x)
{
  std::string gff_file_name;
  FILE *fpCoord;
  char line[MAX_LINE];
  int  i;
  char *feature = new char[FILENAME_MAX];
  int  start,  end;
  char strand, frame;
  string idSo, ontology_term;
  std::vector <int>         vPos;
  std::vector <signed char> vState;

 
  vector<GeneFeature *>::iterator it = geneFeatureSet.getIterator();
  int nbGeneFeature=geneFeatureSet.getNbFeature();
  
  for ( int j=0 ; j < nbGeneFeature ; j++, it++ )
  {
     cout << (*it)->getString() <<endl;	
     strcpy (feature, (*it)->getType().c_str());
     start = (*it)->getLocus()->getStart();
     end = (*it)->getLocus()->getEnd();
     strand = (*it)->getLocus()->getStrand();
     frame = (*it)->getPhase();
     idSo=(*it)->getType();
     ontology_term=(*it)->getAttributes()->getOntologyTerm();
     if ( idSo.find("SO:") == string::npos )
     {
      string tmp=GeneFeatureSet::soTerms_->getIdFromName(idSo);
      idSo=tmp;
     }
     if ( idSo!="SO:0000316" && idSo!="SO:0000204" && idSo!="SO:0000205" )
     {
	continue;
     }
     if (j==0) {
	  vPos.push_back  ( start-1  );
	  vState.push_back( InterGen );
	  if (idSo=="SO:0000204" || idSo=="SO:0000205") { //if UTR
	    if (Todo=="TEST") {
	      vPos.push_back  ( end   );
	      vState.push_back( UTR5F );
	    }
	  }
	  else if (idSo=="SO:0000316" && ontology_term=="SO:0000196") { //E.Init
	    vPos.push_back  ( end    );
	    vState.push_back( InitF1 );
	  }
	  else if (idSo=="SO:0000316" && ontology_term=="SO:0000197") { //E.Term
	    vPos.push_back  ( end    );
	    vState.push_back( InitR1 );
	  }
	  else if (idSo=="SO:0000316" && ontology_term=="SO:0005845")  { //E.Sngl
	    if (strand == '+') {
	      vPos.push_back  ( end    );
	      vState.push_back( InitF1 );
	    }
	    else {
	      vPos.push_back  ( end    );
	      vState.push_back( InitR1 );
	    }
          }
	  else {
	    std::cerr <<"\n Error in gff file : SOFA : "<< idSo<< "; ontology_term : "<<ontology_term
		      <<" feature "<<j<<".\n"
		      <<" WARNING :\n"
		      <<"   - Complete genes only in gff file.\n"
		      <<"   - Feature must be SO:0000204, SO:0000205 or SO:0000316 with ontology term :"
		      <<" SO:0000196, SO:0000197, SO:0005845.\n";
	    exit(2);
	  }
	}
	else {
	  if (idSo=="SO:0000204" || idSo=="SO:0000205") { //if UTR
	    if (Todo=="TEST") {
	      vPos.push_back  ( end   );
	      vState.push_back( UTR5F );
	    }
	  }
	  else if (idSo=="SO:0000316" && ontology_term=="SO:0000196") { //E.Init
	    if (strand == '+') {
	      vPos.push_back  ( end    );
	      vState.push_back( InitF1 );
	    }
	    else {
	      vPos.push_back  ( start-1  );
	      vState.push_back( IntronR1 );
	      vPos.push_back  ( end      );
	      vState.push_back( InitR1   );
	    }
	  }
	  else if (idSo=="SO:0000316" && ontology_term=="SO:0000197")  { //E.Term
	    if (strand == '-') {
	      vPos.push_back  ( end    );
	      vState.push_back( InitR1 );
	    }
	    else {
	      vPos.push_back  ( start-1  );
	      vState.push_back( IntronF1 );
	      vPos.push_back  ( end      );
	      vState.push_back( InitF1   );
	   	}
	  }
	  else if (idSo=="SO:0000316" && ontology_term=="SO:0005845") { //E.Sngl
	    if (strand == '+') {
	      vPos.push_back  ( end      );
	      vState.push_back( InitF1   );
	    }
	    else {
	      vPos.push_back  ( end    );
	      vState.push_back( InitR1 );
	    }
	  }
	  else if (idSo=="SO:0000316" && ontology_term=="SO:0000004")  { //E.Intr
	    if (strand == '+') {
	      vPos.push_back  ( start-1  );
	      vState.push_back( IntronF1 );
	      vPos.push_back  ( end      );
	      vState.push_back( IntrF1   );
	    }
	    else {
	      vPos.push_back  ( start-1  );
	      vState.push_back( IntronR1 );
	      vPos.push_back  ( end      );
	      vState.push_back( InitR1   );
	    }
	  }
	  else {
	    std::cerr <<"\n Error in gff3 file "
		      <<" line "<<j<<".\n"
		      <<" "<<feature<<"("<< idSo<<")"<< ": unknown gff3 feature (SO:0000204, SO:0000205 or SO:0000316 with ontology term :"
		      <<" SO:0000196, SO:0000197, SO:0005845.\n";
	    exit(2);
	  }
	
    }
 }
  // Complete gene ?
  // in Gff3 : Must be sorted by coordinates !
  if (( (idSo=="SO:0000316" && ontology_term=="SO:0000004") ||
      ((idSo=="SO:0000316" && ontology_term=="SO:0000196") &&  strand == '+')  ||
      ((idSo=="SO:0000316" && ontology_term=="SO:0000197") &&  strand == '-')) && !(idSo=="SO:0000204"||idSo=="SO:0000205"))
    {
      std::cerr <<"\n Error in gff file\n";
      std::cerr <<" WARNING : complete genes only in gff file.\n";
      exit(2);
    }
 
  std::cerr <<"done\n";

  delete [] feature;

  gene = new Prediction(0, x->SeqLen,vPos, vState);
  gene->TrimAndUpdate(x);

  vPos.clear();
  vState.clear();
}


// -------------------------------------------
//  SigType_TF (Cf:enum SigType in SensorIF.h)
// -------------------------------------------
const char* SensorTester :: SigType_TF(int i, int pos, const char **sType)
{
  switch (i) {
  case 2:
    *sType = "Start";
    return gene->IsStart(pos);
  case 3:
    *sType = "Stop";
    return gene->IsStop(pos);
  case 4:
    *sType = "Acc";
    return gene->IsAcc(pos);
  case 5:
    *sType = "Don";
    return gene->IsDon(pos);
  default:
    *sType = "Unknown";
    return "- ? -";
  }
}

// -------
//  State
// -------
const char* SensorTester :: State(int pos)
{
  switch (gene->GetStateAtPos(pos)->GetState()) {
  case 0:
    return "InitF";
  case 3:
    return "InitR";
  case 6:
    return "SnglF";
  case 9:
    return "SnglR";
  case 12:
    return "IntrF";
  case 15:
    return "IntrR";
  case 18:
    return "TermF";
  case 21:
    return "TermR";
  case 24:
    return "IntronF";
  case 27:
    return "IntronR";
  case 30:
    return "IG";
  case 31:
    return "UTR";
  case 35:
    return "IntronUTR";
  case -1:        // From the last element of "prediction"
    return "IG";  // to the end of the sequence...
  default:
    return ".";
  }
}

// -------------------------------------------
// Count the number of TP, FP, TN, FN for each threshold
// Just considering canonical sites : ATG for start, 
// TAG, TAA, TGA for stop, AG for Acc, GT, GC for Don
// -------------------------------------------
void SensorTester :: AnalyzeSPSN(void)
{
  int n, k, i, j, l, m, p;
  int nbj, nbl, nb;
  double d, sn, sp, tp, tn, fp, fn;
  bool is_don;
  std::ofstream f, facc, fdon;

  if (EvaluatedType == EVALUATED_SPLICE) {
    nbj = 4; // 4 columns for acc/don, F/R (Scores, IsAPositons, IsAnnotateds)
    nbl = 2; // 2 cumns for acc/don (TP,FP,TN,FN,Nb)
  } else {
    nbj = 2; // 2 column for signal F/R
    nbl = 1; // 1 column for signal
  }

  // prepare vectors to be used as arrays
  for (n=0; n<(int)Thresholds.size(); n++) {
    std::vector<int> v1; TP.push_back(v1); 
    std::vector<int> v2; FP.push_back(v2); 
    std::vector<int> v3; TN.push_back(v3); 
    std::vector<int> v4; FN.push_back(v4); 
    std::vector<int> v5; Nb.push_back(v5);
    for (l=0; l<nbl; l++) {
      TP[n].push_back(0); FP[n].push_back(0);
      TN[n].push_back(0); FN[n].push_back(0);
      Nb[n].push_back(0);
    }
  }

  // sort thresholds from min to max
  sort( Thresholds.begin(), Thresholds.end() );
 
  // Update number of elements for each threshold
  for (k=0; k<(int)Scores.size(); k++)       // for all sequences
    for (i=0; i<(int)Scores[k].size(); i++)  // for all scores of a sequence
      for (j=0; j<nbj; j++) {                // for all type of information
	d = Scores[k][i][j];
	if (d!=0) {
	  m = std::find(Thresholds.begin(),Thresholds.end(),d) - Thresholds.begin();
	  p = ( (j==DON_F || j==DON_R ) ? 1 : 0 );
	  Nb[m][p] += 1; 
	}
      }	

  // for each threshold compute TP FP TN FN considering all the sequences
  for (n=0; n<(int)Thresholds.size(); n++)
    for (k=0; k<(int)Scores.size(); k++) 
      for (i=0; i<(int)Scores[k].size(); i++)
	for(j=0; j<nbj; j++) 
	  if (IsAPositions[k][i][j]) {
	    d = Scores[k][i][j];
	    is_don = ( (j==DON_F || j==DON_R ) ? true : false );
	    UpdateTP_FP_TN_FN( n, is_don, IsAnnotateds[k][i][j], 
			       d != 0 && d >= Thresholds[n]);
	  }

  // output specificity and sensibility 
  f.open((SensorName + ".SpSn").c_str(), std::ios::out);
  if (EvaluatedType == EVALUATED_SPLICE) {
    facc.open((SensorName + ".Acc").c_str(), std::ios::out);
    fdon.open((SensorName + ".Don").c_str(), std::ios::out);
  }

  std::cout << "Thres.\t\tNb\tTP\tFP\tTN\tFN\tSpec.\tSens.\n";
  for (n=0; n<(int)Thresholds.size(); n++) {
    nb = Nb[n][0]; tp = TP[n][0]; fn = FN[n][0]; fp = FP[n][0]; tn = TN[n][0];
    if (EvaluatedType == EVALUATED_SPLICE) 
      {nb += Nb[n][1]; tp += TP[n][1]; fn += FN[n][1]; fp += FP[n][1]; tn += TN[n][1];}
    sn = 100* (tp/(tp+fn));
    sp = 100* (tp/(tp+fp));
    std::cout <<Thresholds[n]<<"\t"<<nb<<"\t"
	      <<tp<<"\t"<<fp<<"\t"<<tn<<"\t"<<fn<<"\t"
	      <<sp<<"\t"<<sn<<"\n";
    if ((tp+fp > MinNumbers) && (tp+fn > MinNumbers))  
      f <<sp<<"\t"<<sn<<"\n";
  }

  if (EvaluatedType == EVALUATED_SPLICE) {
    std::cout << "\nThres.\t\tNbacc\tTPacc\tFPacc\tTNacc\tFNacc\tSpec.\tSens.\n";
    for (n=0; n<(int)Thresholds.size(); n++) {
      sn = 100* (TP[n][0]/(double)(TP[n][0]+FN[n][0]));
      sp = 100* (TP[n][0]/(double)(TP[n][0]+FP[n][0]));
      std::cout <<Thresholds[n]<<"\t"<<Nb[n][0]<<"\t"
	        << TP[n][0]<<"\t"<<FP[n][0]<<"\t"<<TN[n][0]<<"\t"<<FN[n][0]<<"\t"
		<<sp<<"\t"<<sn<<"\n";
      if ((TP[n][0]+FP[n][0] > MinNumbers) && (TP[n][0]+FN[n][0] > MinNumbers)) 
	facc <<sp<<"\t"<<sn<<"\n";
    }
    std::cout << "\nThres.\t\tNbdon\tTPdon\tFPdon\tTNdon\tFNdon\tSpec.\tSens.\n";
    for (n=0; n<(int)Thresholds.size(); n++) {
      sn = 100* (TP[n][1]/(double)(TP[n][1]+FN[n][1]));
      sp = 100* (TP[n][1]/(double)(TP[n][1]+FP[n][1]));
      std::cout <<Thresholds[n]<<"\t"<<Nb[n][1]<<"\t"
		<< TP[n][1]<<"\t"<<FP[n][1]<<"\t"<<TN[n][1]<<"\t"<<FN[n][1]<<"\t"
		<< sp<<"\t"<<sn<<"\n";
      if ((TP[n][1]+FP[n][1] > MinNumbers) && (TP[n][1]+FN[n][1] > MinNumbers)) 
	fdon <<sp<<"\t"<<sn<<"\n";
    }
  }

  f.close();
  if (EvaluatedType == EVALUATED_SPLICE) {
    facc.close();
    fdon.close();
  }
}


// -------------------------------------------
// Update TP, FP, TN, FN and if required TPacc, FPacc, TNacc, FNacc, 
// TPdon, FPdon, TNdon, FNdon
// -------------------------------------------
void SensorTester :: UpdateTP_FP_TN_FN( int no_threshold, 
      bool is_don, bool is_annotated, bool is_detected )
{
  if (!is_don) {
    if (is_detected) {
      if (is_annotated) TP[no_threshold][0] += 1;
      else FP[no_threshold][0] += 1;
    } else {
      if (is_annotated) FN[no_threshold][0] += 1;
      else TN[no_threshold][0] += 1;
    }
  } else {
    if (is_detected) {
      if (is_annotated) TP[no_threshold][1] += 1;
      else FP[no_threshold][1] += 1;
    } else {
      if (is_annotated) FN[no_threshold][1] += 1;
      else TN[no_threshold][1] += 1;
    }
  }
}


// ----------------------------------------------
// Update the list of threshold
// -------------------------------------------
void SensorTester :: UpdateThreshold( double t )
{
  if (isnan(t))
    {std::cerr<<"ERROR: a score is a NaN in SensorTester::UpdateThreshold.\n"; exit(2);}

  if (t != 0) 
    if( find( Thresholds.begin(), Thresholds.end(), t ) == Thresholds.end() )
      Thresholds.push_back( t );
}
