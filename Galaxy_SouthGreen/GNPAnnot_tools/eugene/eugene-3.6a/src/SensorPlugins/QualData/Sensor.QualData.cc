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
// $Id: Sensor.QualData.cc,v 1.1 2008-09-22 12:54:20 dleroux Exp $
// ------------------------------------------------------------------
// File:     Sensor.FrameShift.cc
// Contents: Sensor FrameShift 
// ------------------------------------------------------------------

#include <cassert>
#include <cstring>

#include "Sensor.QualData.h"

extern Parameters PAR;

#define NORM(x,n) (((n)+(Max(-(n),x)))/(n))

/*************************************************************
 **                      SensorQualData                   **
 *************************************************************/

// ----------------------
//  Default constructor.
// ----------------------
SensorQualData :: SensorQualData (int n, DNASeq *X) : Sensor(n), qualPerPos()
{
  type = Type_FS;
}

// ----------------------
//  Default destructor.
// ----------------------
SensorQualData :: ~SensorQualData ()
{
}

/* return probability of wrong base call associated with phred score */
double phred2prob(double Q) {
	static double C = -log(10.)/10.;
	double X = exp(Q*C);
	return (X/(1.+X));
}


void read_qual(std::vector<double>& vec, double A, char*filename, size_t dnasz) {
	std::string header;
	int Q;
	double prob;
	double weight;
	std::ifstream ifs(filename);
	std::getline(ifs, header);	/* skip header */
	vec.clear();
	/*std::cerr << " =====debug===== " << filename << " factor=" << A << std::endl;*/
	while(!ifs.eof()) {
		ifs >> Q;
		prob = phred2prob((double)Q);
		weight = -log(prob)*A;
		/*std::cerr << Q << '\t' << (100.-prob*100.) << "%\t" << weight << std::endl;*/
		vec.push_back(weight);
	}
	if(dnasz+1 == vec.size()) {
		vec.pop_back();	/* because of last space in file ? */
	}
	/*std::cerr << " ===end=debug=== " << vec.size() << std::endl;*/
}

// ----------------------
//  Init start.
// ----------------------
void SensorQualData :: Init (DNASeq *X)
{
	char tempname[FILENAME_MAX+1];
	sprintf(tempname, "%s.qual", PAR.getC("fstname"));
	char* key = (char*)"Qual.factor*";
	if(PAR.probeKey(key)) {
		factor = PAR.getD(key, GetNumber());
	} else {
		factor = 1.0;
	}
	read_qual(qualPerPos, factor, tempname, X->SeqLen);
	/*std::cerr << " ===DNA=debug=== " << X->SeqLen << std::endl;*/
	assert(qualPerPos.size() == X->SeqLen);

	/* FIXME : Plot ? */
	/*if (PAR.getI("Output.graph")) Plot(X);*/
}

// -----------------------
//  GiveInfo frameshift
// -----------------------
void SensorQualData :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
	/* FIXME : offset pos by -1, actual behaviour should DEFINITELY be spec'd and fixed */
	/*assert(pos>0);*/
	double w = -qualPerPos[pos];
	d->sig[DATA::Ins].weight[Signal::Forward] += w;
	d->sig[DATA::Ins].weight[Signal::Reverse] += w;
	d->sig[DATA::Del].weight[Signal::Forward] += w;
	d->sig[DATA::Del].weight[Signal::Reverse] += w;
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorQualData :: Plot(DNASeq *X)
{
	/* TODO ? */
	/* FIXME : ask T.S. about usefulness and contents */
}

// ------------------
//  Post analyse
// ------------------
void SensorQualData :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
	/* TODO ? */
	/* FIXME : ask T.S. about usefulness and contents */
}
