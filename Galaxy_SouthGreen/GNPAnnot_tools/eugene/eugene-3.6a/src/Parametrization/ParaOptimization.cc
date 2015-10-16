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
// $Id: ParaOptimization.cc,v 1.23 2009-07-07 16:15:26 tschiex Exp $
// ------------------------------------------------------------------
// File:     ParaOptimization.cc
// Contents: The ParaOptimization class optimizes eugene parameters
// ------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <time.h>
#include <unistd.h>
#include <stdio.h>

#include "ParaOptimization.h"

#include "../Param.h"
#include "../Prediction.h"
#include "../DNASeq.h"
#include "../MSensor.h"
#include "LineSearch.h"
#include "Genetic.h"

extern Parameters PAR;
extern MasterSensor*    MS;

extern Prediction* Predict (DNASeq* TheSeq, MasterSensor* MSensor);


//-------------------------------------------------------
// Destructor
//-------------------------------------------------------
ParaOptimization::~ParaOptimization(void)
{
    unsigned int i;
    // delete created instances
    for (i=0; i<Algorithms.size(); i++) delete Algorithms[i];
    for (i=0; i<Sequences.size(); i++) delete Sequences[i];
    for (i=0; i<MSensors.size(); i++) delete MSensors[i];
}


//-------------------------------------------------------
// ParaOptimize : Parameters optimization
//-------------------------------------------------------
void ParaOptimization::ParaOptimize (int argc, char * argv [])
{
    std::string filename;
    OptiAlgorithm* algo;
    bool is_chaining = false;
    Regularizer = 0.0;

    // Initialisation
    Init(argc, argv);

    // Optimize parameters
    for (unsigned int i=0; i<Algorithms.size(); i++)
    {
        AlgoIndex = i;
        if (i>0) is_chaining = true;
        Algorithms[i]->Optimize(is_chaining);
    }

    // Write new parameters file
    algo = Algorithms.back();
    filename = PAR.WriteParam(algo->ParaName, algo->Para);
    std::cerr <<std::endl << "A new parameter file " << filename << " is written." << std::endl;
}


//-------------------------------------------------------
// Initialisation
//-------------------------------------------------------
void ParaOptimization::Init(int argc, char * argv [])
{
    std::string algo_name;
    std::string para_optimization_test;

    ExecutableName = argv[0];
    para_optimization_test = (std::string) PAR.getC("ParaOptimization.Test");
    IsTest = ( ( (para_optimization_test == "1") || (para_optimization_test == "TRUE") )  ? true : false);
    // Inhibit graphic mode
    PAR.set("Output.graph", "0");

    // Creation of required instances of optimization algorithms
    algo_name = PAR.getC("ParaOptimization.Algorithm");
    if (algo_name == "GENETIC")
        Algorithms.push_back( new Genetic() );
    else
        if (algo_name == "LINESEARCH")
            Algorithms.push_back( new LineSearch() );
        else
            if (algo_name == "GENETIC+LINESEARCH")
            {
                Algorithms.push_back( new Genetic() );
                Algorithms.push_back( new LineSearch() );
            }
            else
            {
                std::cerr <<"ERROR: Bad optimization algorithm "<<algo_name<<" in the parameter file"<<std::endl;
                exit(100);
            }

    // Get the Margin Penalty Factor
    Regularizer = PAR.getD("ParaOptimization.Regularizer");

    if (!IsTest)
    {
        int sequence;
        TrueCoordFile = PAR.getC("ParaOptimization.TrueCoordFile");

        // Update the sequences list
        std::cout << "Loading sequence(s) file(s) ...";
        for (sequence = optind; sequence < argc ; sequence++)
        {
            Sequences.push_back(new DNASeq(argv[sequence]) );
            SeqNames.push_back(argv[sequence]);
        }
        std::cout << "done (" << Sequences.size() << " sequence(s))" << std::endl;

        // Load the references
        std::cout << "Loading reference(s) ...";
        this->ReadReferences();
        std::cout << "done (" << References.size() << " reference prediction(s)" << std::endl;

        // Update the master sensors list
        for (sequence = 0; sequence<(int)Sequences.size(); sequence++)
        {
            PAR.set("fstname", SeqNames[sequence].c_str());
            MSensors.push_back(new MasterSensor);
            MS = MSensors[sequence];
            MSensors[sequence]->InitMaster(Sequences[sequence]);
        }
    }
}


//-------------------------------------------------------
// ParaEvaluate : Evaluate the parameters Algorithms[AlgoIndex]->Para
//-------------------------------------------------------
double ParaOptimization::ParaEvaluate (bool is_detail_required)
{
    OptiAlgorithm* algo = Algorithms[AlgoIndex];
    double fitness = 0;
    double spg, sng, spe, sne, spn, snn; // specificity/sensibility gene/exon/nt
    Prediction* pred, *ref;
    int TPg = 0, TPe = 0, TPn = 0; // nb of True positive in gene/exon/nucleotide
    int RGnb = 0, PGnb = 0; //number of real genes/predicted genes
    int REnb = 0, PEnb = 0; //number of real exons/predicted exons
    int RNnb = 0, PNnb = 0; //number of nt included in real gene/nt included in predicted gene
    DetailedEvaluation = "";
    // Init the offset value
    int evalOffset;
    char* offsetKey = (char *)"Eval.offset";
    if ( PAR.probeKey(offsetKey) == false)
    {
        evalOffset = INT_MAX/2;
    }
    else
    {
        evalOffset = PAR.getI("Eval.offset");
    }

    if (algo->Para.size() > 0)
    {
        if (!IsTest)
        {
            // Update new value for parameters
            for (unsigned int i=0; i<Algorithms[AlgoIndex]->ParaName.size(); i++)
            {
                PAR.setD(Algorithms[AlgoIndex]->ParaName[i].c_str(), Algorithms[AlgoIndex]->Para[i]);
            }

            // To later Penalize parameter magnitude.
            // Brute code: all optimized parameters are pooled  a la L1
            double mag_penalty = 0.0;
            for (unsigned int i=0; i<Algorithms[AlgoIndex]->ParaName.size(); i++)
            {
                mag_penalty += fabs(Algorithms[AlgoIndex]->Para[i]);
            }

            // update sensors
            for (unsigned int i=0; i<Sequences.size(); i++)
            {
                PAR.set("fstname", SeqNames[i].c_str());
                MS = MSensors[i];
                MSensors[i]->InitSensors(Sequences[i]);
            }
            // Launch prediction of each sequence and evaluate it
            for (unsigned int i=0; i<Sequences.size(); i++)
            {
                PAR.set("fstname", SeqNames[i].c_str());
                pred = Predict(Sequences[i], MSensors[i]);
                ref  = this->References[i];
                // Evaluate the prediction and compute the number of True Positive
                std::vector<int> vEvalPred = pred->Eval(ref, evalOffset);
                // Inc the number of TP/predicted/real gene/exon/nt
                TPg  += vEvalPred[0];
                PGnb += vEvalPred[1];
                RGnb += vEvalPred[2];
                TPe  += vEvalPred[3];
                PEnb += vEvalPred[4];
                REnb += vEvalPred[5];
                TPn  += vEvalPred[6];
                PNnb += vEvalPred[7];
                RNnb += vEvalPred[8];


//         cout << "\n-------------------\n";
//         cout << "PRED:\n";
//         pred->Print();
//         cout << "REF:\n";
//         ref->Print();
// //         cout << "TPg : "<< evalTP[0] <<  " NB gene predit : " << pred->GetGenes(regionStart, regionStop).size();
// //         cout << " NB gene reel : " << ref->nbGene;
//         cout << "\nTPe : "<< evalTP[1] <<  " NB exon predit : " << pred->GetExons(regionStart, regionStop).size();
//         cout << " NB exons reel : " << ref->GetExonNumber();
//         cout <<"\n-------------------\n";

//        if (i!=Sequences.size()-1) fprintf(fp,"\n");
                delete pred;
            }

            if (RGnb > 0) sng = (double)(TPg*100)/(double)RGnb; else sng = 0;
            if (REnb > 0) sne = (double)(TPe*100)/(double)REnb; else sne = 0;
            if (RNnb > 0) snn = (double)(TPn*100)/(double)RNnb; else snn = 0;
            if (PGnb > 0) spg = (double)(TPg*100)/(double)PGnb; else spg = 0;
            if (PEnb > 0) spe = (double)(TPe*100)/(double)PEnb; else spe = 0;
            if (PNnb > 0) spn = (double)(TPn*100)/(double)PNnb; else spn = 0;

            if (is_detail_required)
            {
                std::stringstream DetailedEvaluationStream;
                DetailedEvaluationStream << "** "  << TPg << " genes bien detectes sur " << RGnb <<" avec " << PGnb << " predictions\n** "<< TPe <<" exons bien detectes sur "<< REnb << " avec " << PEnb << " predictions\n** SNG: " << sng << " SPG: "<< spg <<" SNE: "<< sne <<" SPE: "<< spe << "\n";
                DetailedEvaluation = DetailedEvaluationStream.str();
            }

  		if (!is_detail_required)
		{
		    double wsng = PAR.getD("Fitness.wsng");
		    double wsne = PAR.getD("Fitness.wsne");
		    double wsnn = PAR.getD("Fitness.wsnn");
		    double wspg = PAR.getD("Fitness.wspg");
		    double wspe = PAR.getD("Fitness.wspe");
		    double wspn = PAR.getD("Fitness.wspn");
		    if ((wsng + wspg + wsne + wspe + wsnn + wspn) <= 0.0 )
		    {
		        cerr << "Incorrect fitness weights (zero sum) in parameter file\n";
		        exit(2);
		
		    }
		
		    // Compute fitness
		    fitness = pow(pow(sng,wsng)*pow(spg,wspg)*pow(sne,wsne)*pow(spe,wspe)*pow(snn,wsnn)*pow(spn,wspn), 1.0L/(wsng + wspg + wsne + wspe + wsnn + wspn)) -(mag_penalty*Regularizer);
		}
		else
		{
                fitness = 0;
            }
        }
        else
        {
            fitness = 1;
            for (unsigned int i=0; i<algo->Para.size(); i++)
                fitness *= NormalLaw( algo->Para[i] );
        }
    }

    return fitness;
}


//-------------------------------------------------------
// Normal Law
//-------------------------------------------------------
double ParaOptimization::NormalLaw(double x)
{
    return exp(-pow(x-0.5,2)/2) / sqrt(2*M_PI);
}


// -------------------------------------------------------------------------
// Read the TrueCoordFile file containing the real gene structure of the sequences.
// It is composed of many blocks of data, one block describing the genes of one sequence.
// One line describes one gene.
// first field: the name of the sequence, next fields: the exon begin and exon end, etc.
// CAREFUL: File has to finish by an empty line

// Example of file:

// SEQ1 429 545 665 750 835 1125 1255 1591
// SEQ1 -2001 -2342 -2424 -2522 -2614 -2718 -2802 -2894 -3006 -3098 -3201 -3254 -336
// 0 -3461 -3538 -3612 -3746 -3865
// SEQ1 8382 8471 8576 8667 9006 9061 9363 9436
//
// SEQ2 800 894 1200 1312
//
// -------------------------------------------------------------------------
void ParaOptimization :: ReadReferences()
{
    std::string line;
    std::string space(" \t\n\r\f\v");
    std::string buffer;
    std::string seqName;
    std::string fileName;
    int iSeq = -1; // index to scan Sequences and SeqNames vectors

    // Read the true coordinate File
    std::ifstream ifs(this->TrueCoordFile.c_str());

    // read the coordinate file block per block (one block per sequence)
    std::getline(ifs, line);
    while ( !ifs.eof() )
    {
	// Empty line (all space), end of the gene description of a sequence
        if ( line.find_first_not_of(space) == line.npos )
        {
            iSeq++;
            if ( iSeq >= this->SeqNames.size() )
            {
                std::cerr << "\nERROR: true coordinate file doesn't contain the same number of sequences than loaded sequences.\n";
                exit(100);
            }
            // Check the sequence name is identical to the fasta file name
            stringstream ssBuffer(buffer);
            ssBuffer >> seqName;
            fileName = std::string( BaseName( (char*)this->SeqNames[iSeq].c_str() ) );
            if (seqName != fileName.substr(0, seqName.size()) )
            {
                std::cerr << "\nERROR: can't compare the sequence " << seqName << " with the sequence in the file "<< fileName <<".\n";
                exit(100);
            }

            // Create a new Prediction object for the next sequence
            Prediction* p = new Prediction(buffer, this->Sequences[iSeq]);
            this->References.push_back(p); // Save the prediction
            buffer.clear();
        }
        else
        {
            // add line to the buffer
            buffer.append(line);
            buffer.append("\n");
        }

        std::getline(ifs, line);
    }
}

