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
// $Id: Prediction.cc,v 1.63 2010-02-01 16:18:42 sallet Exp $
// ------------------------------------------------------------------
// File:     Prediction.cc
// Contents: class Prediction
// ------------------------------------------------------------------

#include <iostream>
#include <climits>
#include "Prediction.h"
#include "MSensor.h"
#include "DNASeq.h"

extern Parameters   PAR;

// -----------------------
//  Prediction: IsOriginal checks is the 1st gene of the current
//  alternative prediction does not appear either as a gene in the
//  optimal prediction or as an already predicted alternative gene
// -----------------------
bool Prediction :: IsOriginal ( Prediction* optPred, std::vector <Prediction*>& altPreds, int seuil )
{
	Gene *thisGene = this->vGene[0];
	Gene *otherGene;
	int mismatch;

	// First check in the optimal prediction.
	int idx;

	for ( idx = 0; idx < optPred->vGene.size(); idx++ )
	{
		mismatch = thisGene->isDifferent ( * ( optPred->vGene[idx] ), seuil );
		if ( ( mismatch >= 0 ) && ( mismatch <= seuil ) )
			return false;
	}
	// Now alternative predictions. They contain just one gene

	for ( idx = 0; idx < altPreds.size(); idx++ )
	{
		mismatch = thisGene->isDifferent ( * ( altPreds[idx]->vGene[0] ), seuil );
		if ( ( mismatch >= 0 ) && ( mismatch <= seuil ) )
			return false;
	}

	return true;
}

// -----------------------
//  ESTScan: fills ESTMatch with just Gap and Hit. Used to detect supported UTR
// -----------------------
void Prediction :: ESTScan()
{
	int i, NumEST, EstInstanceNb;
	Hits *ThisEST = NULL, *AllEST = NULL;
	Block *ThisBlock = NULL;
	char tempname[FILENAME_MAX+1];
	FILE *fEST;
	char  *fileExt;           // extension of an EST file
	std::string inputFormat;  // format of an EST file
	ESTMatch = NULL;

	// If no EST data, do  nothing
	if ( ! ( PAR.getI ( "Sensor.Est.use", 0, PAR.getI ( "EuGene.sloppy" ) ) ) ) return;

	ESTMatch = new unsigned char[X->SeqLen+1];
	for ( i = 0; i <= X->SeqLen; i++ ) ESTMatch[i] = 0;

	NumEST        = 0;
	EstInstanceNb = PAR.getI ( "Sensor.Est.use" ); // number of EST instances

	for ( int i=0; i<EstInstanceNb; i++ )
	{
		fileExt = PAR.getC ( "Est.FileExtension", i );
		strcpy ( tempname, PAR.getC ( "fstname" ) );
		strcat ( tempname, fileExt );
		inputFormat = to_string ( PAR.getC ( "Est.format", i, 1 ) );
		if ( inputFormat == "GFF3" )
		{
			strcat ( tempname,".gff3" );
			GeneFeatureSet * geneFeatureSet = new GeneFeatureSet ( tempname );
			AllEST = AllEST->ReadFromGeneFeatureSet ( *geneFeatureSet, &NumEST, -1, 0, X );
			delete geneFeatureSet;
		}
		else
		{
			fEST = FileOpen ( NULL, tempname, "r", PAR.getI ( "EuGene.sloppy" ) );
			if ( fEST )
			{
				AllEST = AllEST->ReadFromFile ( fEST, &NumEST, -1, 0,X->SeqLen );
				fclose ( fEST );
			}
		}
	}

	ThisEST = AllEST; // start from first one.
	while ( ThisEST )
	{
		ThisBlock = ThisEST->Match;
		while ( ThisBlock )
		{
			//Mark Hit
			for ( i = ThisBlock->Start; i <= ThisBlock->End; i++ ) ESTMatch[i] |= Hit;
			// A Gap ?
			if ( ThisBlock->Prev != NULL ) // Mark Gap
				for ( i = ThisBlock->Prev->End; i <= ThisBlock->Start; i++ ) ESTMatch[i] |= Gap;

			ThisBlock = ThisBlock->Next;
		}
		ThisEST = ThisEST->Next;
	}
	delete AllEST;
}

/*************************************************************
 **                       Feature Object                    **
 *************************************************************/
// ------------------------
//  Default constructor.
// ------------------------
Feature :: Feature ()
{
	number   = 0;
	start    = -1;
	end      = -1;
	strand   = '.';
	phase    = 9;
	frame    = 9;
	featureState = new State ( -1 );
}

// ------------------------
//  Constructor.
// ------------------------
Feature::Feature ( signed char n, int s, int e )
{
	number   = 0;
	featureState = new State ( n );
	start    = s;
	end      = e;
	strand   = featureState->GetStrand();
	frame    = featureState->GetFrame();
	phase    = 9;
}

// ------------------------
//  Default destructor.
// ------------------------
Feature :: ~Feature ()
{
	delete featureState;
}

// ------------------------
//  Return true if the feature overlaps the feature f in the same strand
// ------------------------
bool Feature :: Overlap ( const Feature& f )
{
	if ( ( this->strand == f.strand ) &&
	        ( this->end    >= f.start )  &&
	        ( this->start  <= f.end ) )
	{
		return true;
	}
	return false;
}

// ------------------------
//  Return true if the feature is a coding exon
// ------------------------
bool Feature :: IsCodingExon()
{
	return this->featureState->IsCodingExon();
}

// ------------------------
//  Return true if the feature is a non coding rna
// ------------------------
bool Feature :: IsNpcRna()
{
	return this->featureState->IsNpcRna();
}

// ------------------------
//  Return true if the feature is an intron (incuding intron in utr region)
// ------------------------
bool Feature :: IsIntron()
{
	return this->featureState->IsIntron();
}
// ------------------------
//  Return true if the feature is an UTR
// ------------------------
bool Feature :: IsUTR()
{
	return this->featureState->IsUTR();
}

// ------------------------
//  Return true if the feature is an intron in a utr region
// ------------------------
bool Feature :: IsUTRIntron()
{
	return this->featureState->IsUTRIntron();
}

// ------------------------
//  Return true if the feature is transcribed but no spliced  (UTR or Exon)
// ------------------------
bool Feature::IsTranscribedAndUnspliced()
{
	return this->featureState->IsTranscribedAndUnspliced();
}

// ------------------------
//  Return true if the feature is an intergenic region
// ------------------------
bool Feature::IsIntergenic()
{
	return this->featureState->IsIntergenic();
}


// ------------------------
//  Deduce the phase of the feature according to the frame value
// ------------------------
void Feature::ComputePhase ( int seqLength )
{
	if ( this->featureState->InStartStopRegion() )
	{
		int startFeat = ( this->strand == '+' ) ? this->start : ( seqLength - this->end + 1 );
		this->phase   = ( abs ( this->frame ) - ( startFeat%3 ) + 3 ) %3;
	}
	else
	{
		std::cerr << "Careful: you try to compute the phase of a feature which is not an exon or an intron!\n";
	}
}


/*************************************************************
 **                         Gene Object                     **
 *************************************************************/
// ------------------------
//  Default constructor.
// ------------------------
Gene :: Gene ()
{
	complete   = 0;
	isvariant  = false;
	isNpcRna    = false;
	hasvariant = 0;
	tuStart    = 0;
	tuEnd      = 0;
	clear();
}

// ------------------------
// Constructor
// line contains the exon structure of the gene. First field is the name of the sequence
// line example:  SEQ 429 545 665 750 835 1125 1255 1591
// ------------------------
Gene :: Gene ( std::string line, DNASeq& seq )
{
	std::vector<int> exons; // vector containing the exon positions
	int lpos, rpos;         // left and right positions of an exon
	int intrLPos, intrRPos; // left and right positions of an intron
	signed char state;      // nature of the exon/intron
	int rest_codon_lg = 0;  // number of nt in the exon belonging to the previous codon
	int exon_length   = 0;
	int seqLength     = seq.SeqLen;

	//read the first field, the name of the sequence
	std::string field;
	istringstream iss ( line );
	iss >> field;

	// save the exon positions in the vector
	while ( iss >> lpos )
	{
		iss >> rpos; // read the end position of the exon
		exons.push_back ( lpos );
		exons.push_back ( rpos );
	}

	if ( exons.size() > 0 )
	{
		// Gene on the forward strand
		if ( exons[0] >= 0 )
		{
			// Add one exon and the next intron
			for ( int i=0; i < exons.size()-1; i+=2 )
			{
				lpos = exons[i];
				rpos = exons[i+1];
				// Compute the state of the exon
				if ( exons.size() == 2 )            state = SnglF1; // Case single exon
				else if ( i == 0 )                  state = InitF1; // Case first exon
				else if ( ( i+1 ) == exons.size()-1 ) state = TermF1; // Case last exon
				else                                state = IntrF1; // Internal exon
				// Compute the frame
				state += ( lpos + rest_codon_lg + 2 ) % 3; // [Init|Term|Intr]F[1|2|3]]
				this->AddFeature ( state, lpos, rpos ); // Add the exon

				// Use to compute the frame of the intron and of the next exon
				exon_length   = rpos - lpos + 1 - rest_codon_lg;
				rest_codon_lg = ( 3 - exon_length%3 ) %3;

				// Case there is a next intron
				if ( exons.size() != 2 && ( ( i+1 ) < exons.size()-1 ) )
				{
					intrLPos    = rpos+1;
					intrRPos    = exons[i+2]-1;
					// compute the state of the intron
					if ( rest_codon_lg == 0 )         state = IntronF1;
					else if ( rest_codon_lg == 2 )
					{
						if ( seq[intrLPos-2] == 't' ) state = IntronF2T;
						else                        state = IntronF2;
					}
					else
					{
						if ( seq[intrLPos-3] == 't' && seq[intrLPos-2] == 'g' ) state = IntronF3TG;
						else if ( seq[intrLPos-3] == 't' && seq[intrLPos-2] == 'a' ) state = IntronF3TA;
						else  state = IntronF3;
					}
					this->AddFeature ( state, intrLPos, intrRPos );
				}
			}
		}
		else // Strand -
		{
			std::vector<signed char> vStates; // vector use to save states of the features
			std::vector<int> vStart; // vector use to save start pos of the features
			std::vector<int> vStop; // vector use to save stop pos of the features
			// scan the exons starting by the end
			for ( int i=exons.size()-1; i > 0; i-=2 )
			{
				lpos = abs ( exons[i-1] );
				rpos = abs ( exons[i] );
				// Compute the nature of the exons
				if ( exons.size() == 2 )        state = SnglR1; // Case single exon
				else if ( ( i-1 ) == 0 )          state = TermR1; // Case last exon
				else if ( i == exons.size()-1 ) state = InitR1; // Case first exon
				else                            state = IntrR1; // Case internal exon
				// Compute the frame
				state += ( seqLength - rpos + rest_codon_lg ) %3; //[Init|Term|Intr]R[1|2|3]]
				vStates.push_back ( state );
				vStart.push_back ( lpos );
				vStop.push_back ( rpos );
				// will be use to compute the frame of the next exon
				exon_length   = rpos - lpos + 1 - rest_codon_lg;
				rest_codon_lg = ( 3 - exon_length%3 ) %3;

				// Case there is a next intron
				if ( exons.size() != 2 && ( ( i-1 ) > 0 ) )
				{
					intrLPos    = abs ( exons[i-2] ) +1;
					intrRPos    = lpos-1;
					// compute the state of the intron
					if ( rest_codon_lg == 0 )       state = IntronR1;
					else if ( rest_codon_lg == 2 )
					{
						if ( seq[intrLPos-3] == 't' && seq[intrLPos-2] == 'c' ) state = IntronR2AG;
						else                                                  state = IntronR2;
					}
					else
					{
						if ( seq[intrLPos-2] == 'c' )  state = IntronR3G;
						else if ( seq[intrLPos-2] == 't' ) state = IntronR3A;
						else                        state = IntronR3;
					}
					vStates.push_back ( state );
					vStart.push_back ( intrLPos );
					vStop.push_back ( intrRPos );
				}
			}
			// Add the features
			while ( !vStates.empty() )
			{
				this->AddFeature ( vStates.back(), vStart.back(), vStop.back() );
				vStates.pop_back();
				vStart.pop_back();
				vStop.pop_back();
			}
		}
	}

	this->Update ( seqLength );
}

// ------------------------
//  Default destructor.
// ------------------------
Gene :: ~Gene ()
{
	for ( int i=0; i<nbFea(); i++ )
	{
		delete  vFea[i];
		vFea[i] = NULL;
	}
	vFea.clear();
}

//
// Return the variant code when alternative forms are predicted
//
char Gene :: GetVariantCode ()
{
	char code_variant = '\0';

	if ( this->hasvariant )
	{
		code_variant = ( this->isvariant ) ? ( char ) ( 'b'+ ( this->hasvariant-1 ) ) : 'a';
	}
	return code_variant;
}

// ------------------------
//  Default cleaner.
// ------------------------
void Gene :: clear()
{
	exNumber   = exLength   = 0;
	inNumber   = inLength   = 0;
	utrLength  = mrnaLength = geneLength = 0;
	cdsStart   = cdsEnd     = trStart    = trEnd = -1;
}

// ------------------------
//  comparison: same CDS ?
// Assumes complete genes with similar stucture (number of exons, positions)
// returns the cumulated diffence in position or a negative number if different structure.
// ------------------------
int Gene :: isDifferent ( const Gene& o, int threshold )
{
	int idxo = 0;
	int idxt = 0;
	int nDiff = 0;

	//first coding exon in each
	while ( ! o.vFea[idxo]->IsCodingExon() )     idxo++;
	while ( ! this->vFea[idxt]->IsCodingExon() ) idxt++;
	int shift = 0;

	while ( ( idxo+shift < o.vFea.size() ) &&
	        ( idxt+shift < this->nbFea() ) &&
	        ( !o.vFea[idxo+shift]->IsUTR() ) && // not UTR
	        ( !this->vFea[idxt]->IsUTR() )    // not UTR
	      )
	{
		if ( ( o.vFea[idxo+shift]->featureState->GetState() == this->vFea[idxt+shift]->featureState->GetState() ) &&
		        ( abs ( o.vFea[idxo+shift]->start - this->vFea[idxt+shift]->start ) <= threshold ) &&
		        ( abs ( o.vFea[idxo+shift]->end - this->vFea[idxt+shift]->end )  <= threshold ) )
		{
			nDiff += abs ( o.vFea[idxo+shift]->start - this->vFea[idxt+shift]->start );
			nDiff += abs ( o.vFea[idxo+shift]->end - this->vFea[idxt+shift]->end );
			shift++;
		}
		else return -1;
	}
	// we assume the gene structures satisfy the structural constraints
	// of complete genes. If the final exon is identical, it is a term
	// exon and the 2 CDS are identical.
	return nDiff;
}

// ------------------------
//  comparison: same CDS ?
// Assumes complete genes
// ------------------------
bool Gene :: operator== ( const Gene& o )
{
	int idxo = 0;
	int idxt = 0;

	//first coding exon in each
	while ( !o.vFea[idxo]->IsCodingExon() ) idxo++;
	while ( !this->vFea[idxt]->IsCodingExon() ) idxt++;
	int shift = 0;

	while ( ( idxo+shift < o.vFea.size() ) &&
	        ( idxt+shift < this->nbFea() ) &&
	        ( !o.vFea[idxo+shift]->IsUTR() ) &&  // not UTR
	        ( !this->vFea[idxt+shift]->IsUTR() ) // not UTR
	      )
	{
		if ( ( o.vFea[idxo+shift]->start != this->vFea[idxt+shift]->start ) ||
		        ( o.vFea[idxo+shift]->end   != this->vFea[idxt+shift]->end )   ||
		        ( o.vFea[idxo+shift]->featureState->GetState() != this->vFea[idxt+shift]->featureState->GetState() ) ) return false;
		shift++;
	}
	// we assume the gene structures satisfy the structural constraints
	// of complete genes. If the final exon is identical, it is a term
	// exon and the 2 CDS are identical.

	return true;
}

// ------------------------
//  Compare the exons of the two genes
// ------------------------
bool Gene :: HasSameExons ( const Gene& o )
{
	int idxo = 0; // index of the o gene
	int idxt = 0; // index of the current gene

	while ( ( idxo < o.vFea.size() ) && ( idxt < this->nbFea() ) )
	{
		if ( !o.vFea[idxo]->IsCodingExon() )
		{
			idxo++;
		}
		else if ( !this->vFea[idxt]->IsCodingExon() )
		{
			idxt++;
		}
		else
		{
			if ( ( o.vFea[idxo]->featureState->GetState() != this->vFea[idxt]->featureState->GetState() ) ||
			        ( o.vFea[idxo]->start != this->vFea[idxt]->start ) ||
			        ( o.vFea[idxo]->end   != this->vFea[idxt]->end ) )
			{
				return false;
			}
			idxo++;
			idxt++;
		}
	}
	return true;
}


// ------------------------
//  Add Feature.
// ------------------------
void Gene :: AddFeature ( signed char state, int start, int end )
{
	State featState ( state );
	vFea.push_back ( new Feature ( state, start, end ) );

	// Complete : +1 pour un exon init|sngl et +2 pour term|sngl
	if ( featState.IsCodingExon() )
	{
		if ( featState.IsInitExon() ) complete += 1;
		else if ( featState.IsSnglExon() ) complete += 3;
		else if ( featState.IsTermExon() ) complete += 2;
	}
	if ( featState.IsNpcRna() )
	{
		this->isNpcRna = true;
		complete += 3;
	}
}

// ------------------------
//  Update.
// ------------------------
void Gene :: Update ( int seqLength )
{
	int forward   = ( vFea[0]->strand == '+' ) ? 1 : 0;

	clear();

	// Compute Gene informations: [ex|in]Number, [ex|in|utr]Length, ..
	for ( int i=0; i<nbFea(); i++ )
	{
		if ( vFea[i]->IsCodingExon()  || vFea[i]->IsNpcRna() )
		{
			exNumber++;
			exLength += vFea[i]->end - vFea[i]->start + 1;
			if ( cdsStart == -1 ) cdsStart = vFea[i]->start - 1;
			cdsEnd = vFea[i]->end - 1;
		}
		else if ( vFea[i]->IsIntron() )
		{
			inNumber++;
			inLength += vFea[i]->end - vFea[i]->start + 1;
		}
		else if ( vFea[i]->IsUTR() )
		{
			utrLength += vFea[i]->end - vFea[i]->start + 1;
		}
	}
	trStart = vFea[0]->start - 1;
	trEnd   = vFea[nbFea()-1]->end - 1;
	if ( vFea[nbFea()-1]->featureState->GetState() < UTR5F ) cdsEnd = trEnd;

	mrnaLength = exLength   + utrLength;
	geneLength = mrnaLength + inLength;
	if ( trStart == -1 )
	{
		cdsStart = trStart = 0;
	}
	if ( trEnd   == -1 )
	{
		cdsEnd   = trEnd   = 0;
	}

	// The number of the current exon
	int nb = ( forward ) ? 1 : exNumber;

	for ( int i=0; i<nbFea(); i++ )
	{
		// keep the exons
		if ( !vFea[i]->IsCodingExon() ) continue;

		// Update the exon number
		vFea[i]->number = nb;
		( forward ) ? nb++ : nb--;

		// Compute the phase of the exon
		vFea[i]->ComputePhase ( seqLength );
		// If the previous feature is an intron, update the value of its phase
		if ( forward && ( i > 0 ) && vFea[i-1]->featureState->IsIntronInStartStopRegion() )
		{
			vFea[i-1]->phase = vFea[i]->phase;
		}
		else if ( !forward && ( i < nbFea()-1 ) && vFea[i+1]->featureState->IsIntronInStartStopRegion() )
		{
			vFea[i+1]->phase = vFea[i]->phase;
		}
	}
}

// ------------------------
//  PrintGene.
// ------------------------
void Gene :: PrintInfo ( FILE *F, int nb, char *seqName )
{
	bool nofirst = false;
	int i, state;
	char comp[10];
	switch ( complete )
	{
		case 1:
			strcpy ( comp, "Cfrag" );
			break;
		case 2:
			strcpy ( comp, "Nfrag" );
			break;
		case 3:
			strcpy ( comp, "Full" );
			break;
		default:
			strcpy ( comp, "?" );
	}

	if ( this->isNpcRna ) // Just display start and stop of the rna
	{
		fprintf ( F,"%s.%d  \tEuGene_misc\tncRNA\t%d\t%d\t%d\t%c\t.\t",
		          seqName, nb, cdsStart+1, cdsEnd+1, exLength, vFea[0]->strand );
		fprintf ( F,"%s\t", comp );
	}
	else
	{
		// CDS
		fprintf ( F,"%s.%d  \tEuGene_misc\tCDS\t%d\t%d\t%d\t%c\t.\t",
		          seqName, nb, cdsStart+1, cdsEnd+1, exLength, vFea[0]->strand );
		fprintf ( F,"%s\t", comp );
		for ( i=0; i<nbFea(); i++ )
		{
			if ( vFea[i]->IsCodingExon() )
			{
				if ( nofirst ) fprintf ( F,"," );
				else         nofirst = true;
				fprintf ( F,"%d..%d",vFea[i]->start,vFea[i]->end );
			}
		}
		fprintf ( F,"\n" );

		// Gene
		fprintf ( F,"%s.%d  \tEuGene_misc\tGene\t%d\t%d\t%d\t%c\t.\t",
		          seqName, nb, trStart+1, trEnd+1,
		          geneLength, vFea[0]->strand );
		fprintf ( F,"%s\t", comp );
		nofirst = false;
		for ( i=0; i<nbFea(); i++ )
		{
			if ( vFea[i]->IsCodingExon() || vFea[i]->IsUTR() )
			{
				char sep = ',';
				// Si exon et etat d'avant utr OU
				// si utr  et etat d'avant exon
				if ( ( vFea[i]->IsCodingExon()  && ( i != 0 && ( vFea[i-1]->IsUTR() || vFea[i-1]->IsUTRIntron() ) ) ) ||
				        ( ( vFea[i]->IsUTR() ||vFea[i]->IsUTRIntron() )  && ( i != 0 && vFea[i-1]->IsCodingExon() ) ) )
					sep = ':';
				if ( nofirst ) fprintf ( F,"%c",sep );
				else         nofirst = true;
				fprintf ( F,"%d..%d",vFea[i]->start,vFea[i]->end );
			}
		}
	}


	fprintf ( F,"\n" );
}

// ------------------------
//  Return the number of exons
// ------------------------
int Gene :: GetExonNumber()
{
	return exNumber;
}

// ------------------------
//  Print Gene data
// ------------------------
void Gene::Print()
{
	signed char state;
	int start, end, phase, frame;

	for ( int i=0; i < this->nbFea(); i++ )
	{
		state = this->vFea[i]->featureState->GetState();
		start = this->vFea[i]->start;
		end   = this->vFea[i]->end;
		phase = this->vFea[i]->phase;
		frame = this->vFea[i]->frame;
		printf ( "state %d, start %d, end %d, phase %d, frame %d\n", state, start, end, phase, frame );
	}
}

/*************************************************************
 **                        Prediction                       **
 ************************************************************/
// ------------------------
//  Default constructor.
// ------------------------
Prediction :: Prediction ()
{
	clear();
}

// ------------------------
//  constructor.
// ------------------------
Prediction :: Prediction ( int From, int To, std::vector <int> vPos,
                           std::vector <signed char> vState )
{
	clear();
	int start = 1+From;

	for ( int i=0; i< ( int ) vPos.size(); i++ )
	{
		if ( i!=0 ) start = vPos[i-1] + 1;

		// Si state=intergenique OU start > seqlen on passe
		if ( vState[i] == InterGen || start > vPos[vPos.size()-1] ) continue;

		// Si premier item OU nvx gene avec ou sans utr
		// (=> etat=(UTR5F|UTR3R) && etat precedent est intergenique)
		if ( i==0 ||
		        ( ( vState[i] == UTR5F || vState[i] == UTR3R || vState[i] == RnaF || vState[i] == RnaR ) && vState[i-1] == InterGen ) ||
		        ( vState[i] <= TermR3 && vState[i-1] == InterGen ) )
		{
			nbGene++;
			vGene.push_back ( new Gene() );
		}
		vGene[nbGene-1]->AddFeature ( vState[i], start, vPos[i] );
	}
}

// --------------------------
// constructor reading description text
// Each description line describes one gene of the sequence:
// first field the name of the sequence, next fields the exon begin and exon end, etc.
// Example of desc:

// SEQ1 429 545 665 750 835 1125 1255 1591
// SEQ1 -2001 -2342 -2424 -2522 -2614 -2718 -2802 -2894 -3006 -3098 -3201 -3254 -336
// 0 -3461 -3538 -3612 -3746 -3865
// SEQ1 8382 8471 8576 8667 9006 9061 9363 9436
// --------------------------
Prediction :: Prediction ( const std::string & desc, DNASeq* sequence )
{
	std::string line;
	std::stringstream ss ( desc );

	clear();

	// Init the dna sequence
	X = sequence;

	// Create new Gene objects and save them in the vGene vector
	while ( std::getline ( ss, line, '\n' ) )
	{
		// SEQ1 429 545 665 750 835 1125 1255 1591
		Gene* gene = new Gene ( line, * ( this->X ) );
		this->nbGene++;
		this->vGene.push_back ( gene );
	}
}


// ------------------------
//  Default destructor.
// ------------------------
Prediction :: ~Prediction ()
{
	for ( int i=0; i<nbGene; i++ )
	{
		delete vGene[i];
		vGene[i] = NULL;
	}
	vGene.clear();
	delete [] ESTMatch;

}
// ------------------------
//  Initializer
// ------------------------
void Prediction :: clear()
{
	MS = NULL;
	X  = NULL;
	ESTMatch = NULL;
	seqName[0]  = '\000';
	nbGene      = 0;
	optimalPath = 0;
}
// --------------------------
// Sanity check
// --------------------------
void Prediction :: SanityCheck()
{
	for ( int i=0; i<nbGene; i++ )
	{
		std::vector <Feature*>::iterator featindex;
		std::vector <Feature*>::reverse_iterator rfeatindex;

		for ( featindex = vGene[i]->vFea.begin(); featindex != vGene[i]->vFea.end(); featindex++ )
		{
			//printf("start %i end %i trStart %i trEnd %i\n",(*featindex)->start,(*featindex)->end,vGene[i]->trStart,vGene[i]->trEnd);

			assert ( ( *featindex )->start >= vGene[i]->trStart+1 );
			assert ( ( *featindex )->end <= vGene[i]->trEnd+1 );
		}
	}
	printf ( "Prediction Sanity Check passed\n" );
}
// --------------------------
//  TrimAndUpdate:
//  - modifies the prediction to Trim UTR according to EST support
//    - computes all statistics (length and intron numbers)
//    - deletes short or partial genes if needed
// --------------------------
void Prediction :: TrimAndUpdate ( DNASeq* x )
{
	int  state,start,end;
	State* featState;
	X = x;

	if ( PAR.getI ( "Output.UTRtrim" ) )
	{

		// Fill in ESTMatch data if any
		ESTScan();

		for ( int i=0; i<nbGene; i++ )
		{
			// delete unsupported or trim
			std::vector <Feature*>::iterator featindex;
			std::vector <Feature*>::reverse_iterator rfeatindex;

			for ( featindex = vGene[i]->vFea.begin(); featindex != vGene[i]->vFea.end(); )
			{
				state = ( *featindex )->featureState->GetState();
				featState = ( *featindex )->featureState;
				//if (state <= InterGen) break; // Not UTR or UTRIntron
				if ( !featState->IsUTRIntron() && !featState->IsUTR() ) break;// Not UTR or UTRIntron
				start = ( *featindex )->start;
				end = ( *featindex )->end;
				if ( UTRCheckAndTrim ( &start,&end,featState->GetState() ) )
				{
					( *featindex )->start = start;
					break;
				}
				else
				{
					//          printf("Erasing %d-%d (state %d)\n",(*featindex)->start,(*featindex)->end,(*featindex)->state);
					featindex = vGene[i]->vFea.erase ( featindex );
				}
			}

			// Same on right
			for ( rfeatindex = vGene[i]->vFea.rbegin(); rfeatindex != vGene[i]->vFea.rend(); )
			{
				//state = (*rfeatindex)->featureState->GetState();
				featState = ( *rfeatindex )->featureState;
				// if (state <= InterGen) break; // Not UTR or UTRIntron
				if ( !featState->IsUTRIntron() && !featState->IsUTR() ) break;
				start = ( *rfeatindex )->start;
				end   = ( *rfeatindex )->end;
				if ( UTRCheckAndTrim ( &start,&end,featState->GetState() ) )
				{
					//          printf("Triming %d-%d (state %d) to E%d\n",(*featindex)->start,(*featindex)->end,
					//            (*featindex)->state,end);
					( *rfeatindex )->end = end;
					break;
				}
				else
				{
					//          printf("Erasing %d-%d (state %d)\n",(*rfeatindex)->start,(*rfeatindex)->end,(*rfeatindex)->state);
					rfeatindex++;
					vGene[i]->vFea.pop_back();
				}
			}
		}
	}

	if ( ESTMatch ) delete [] ESTMatch; // ESTScan() delete
	ESTMatch = NULL;

	UpdateAndDelete();

}

// --------------------------
//  Gene update and deletion of short or fragmentary genes
// --------------------------
void Prediction :: UpdateAndDelete()
{
	std::vector <Gene*>::iterator geneindex;
	int MinCDSLen = PAR.getI ( "Output.MinCDSLen" );
	int RemoveFrags = PAR.getI ( "Output.RemoveFrags" );
	int gIdx = 0;

	for ( geneindex = vGene.begin(); geneindex != vGene.end(); )
	{
		int empty_gene = ( ( *geneindex )->nbFea() < 1 );
		if ( ! empty_gene )
			( *geneindex )->Update ( this->X->SeqLen );
		if ( ( *geneindex )->isNpcRna == false &&
		        ( ( *geneindex )->exLength <= MinCDSLen ||  empty_gene || ( RemoveFrags && ( ( *geneindex )->complete != 3 ) ) ) )
		{
			nbGene--;
			geneindex = vGene.erase ( geneindex );
		}
		else
		{
			( *geneindex )->geneNumber = gIdx;
			gIdx++;
			geneindex++;
		}
	}
}
// --------------------------
//  Delete genes outside of range:
//  TS modified. Keeps only one gene with maximum overlap with the region.
// --------------------------
void Prediction :: DeleteOutOfRange ( int s,int e )
{
	std::vector <Gene*>::iterator geneindex;
	int overlap = 0;
	int left,right;

// 1st compute largest Olap
	for ( geneindex = vGene.begin(); geneindex != vGene.end(); )
	{
		left =  Max ( s, ( *geneindex )->trStart ); // left of overlap region (if any)
		right = Min ( e, ( *geneindex )->trEnd );   // right of overlap region (if any)
		if ( ( right - left ) > overlap ) overlap = right - left;

		if ( ( right - left ) <= 0 )  // delete non overlapping genes
		{
			nbGene--;
			geneindex = vGene.erase ( geneindex );
		}
		else geneindex++;
	}

// then keep only one olap'ing gene
	for ( geneindex = vGene.begin(); geneindex != vGene.end(); )
	{
		left =  Max ( s, ( *geneindex )->trStart ); // left of overlap region (if any)
		right = Min ( e, ( *geneindex )->trEnd );   // right of overlap region (if any)
		if ( ( right - left ) < overlap ) // throw away
		{
			nbGene--;
			geneindex = vGene.erase ( geneindex );
		}
		else
		{
			overlap++; // will never be found again. Only one gene kept
			geneindex++;
		}
	}
}


// --------------------------
//  print prediction (master)
// --------------------------
void Prediction :: Print ( DNASeq* x, MasterSensor *ms, FILE *OPTIM_OUT, const char append )
{

	MS = ms;
	FILE *OUT;
	char nameformat[20];
	char outputFormat[20];
	int  trunclen =      PAR.getI ( "Output.truncate" );
	strcpy ( outputFormat, PAR.getC ( "Output.format" ) );

	const char *mode;
	std::ios_base::openmode ccmode;

	if ( append )
	{
		mode = "ab";
		ccmode = std::ofstream::out |std:: ofstream::app;
	}
	else
	{
		mode = "wb";
		ccmode = std::ofstream::out | std::ofstream::trunc;
	}

	// SeqName
	if ( trunclen ) sprintf ( nameformat,"%%%d.%ds",trunclen,trunclen );
	else strcpy ( nameformat,"%s" );
	sprintf ( seqName, nameformat, X->Name );

	// If optimisation mode
	if ( OPTIM_OUT != NULL )
	{
		PrintEgnL ( OPTIM_OUT, seqName );
	}
	else
	{
		for ( int i=0; i<strlen ( outputFormat ); i++ )
		{
			char filename[FILENAME_MAX];
			strcpy ( filename,PAR.getC ( "prefixName" ) );
			switch ( outputFormat[i] )
			{
				case 'a':
					OUT = FileOpen ( NULL, strcat ( filename, ".egn.ara" ), mode );
					PrintEgnL ( OUT, seqName, 1 );
					fclose ( OUT );
					break;
				case 'd':
					OUT = FileOpen ( NULL, strcat ( filename, ".egn.debug" ), mode );
					PrintEgnD ( OUT );
					fclose ( OUT );
					break;
				case 'g':
					OUT = FileOpen ( NULL, strcat ( filename,".gff1" ), mode );
					PrintGff ( OUT, seqName );
					fclose ( OUT );
					{
						std::string filename_gff3 ( PAR.getC ( "prefixName" ) );
						filename_gff3 += ".gff3";
						std::ofstream out ( filename_gff3.c_str(),ccmode );
						//  , std::ios_base::binary);
						//le 1 est a verifier
						if ( ! append )
						{
							out << Gff3Line::header ( X->Name, 1, X->SeqLen );
						}
						PrintGff3 ( out,seqName,append );
						out.close();
					}
					break;
				case 'h':
					OUT = FileOpen ( NULL, strcat ( filename, ".html" ), mode );
					PrintHtml ( OUT, seqName );
					fclose ( OUT );
					break;
				case 'l':
					OUT = FileOpen ( NULL, strcat ( filename, ".egn" ), mode );
					PrintEgnL ( OUT, seqName );
					fclose ( OUT );
					break;
				case 's':
					OUT = FileOpen ( NULL, strcat ( filename, ".egn.short" ), mode );
					PrintEgnS ( OUT );
					fclose ( OUT );
					break;
				case 'o':
					fprintf ( stderr,"\n    Seq         Type    S       Lend    Rend   Length  Phase   Frame      Ac      Do      Pr\n\n" );
					PrintEgnL ( stdout, seqName );
					fprintf ( stdout, "\n" );
					break;
			}
		}
	}
}
// ------------------------
//  print prediction (GFF)
// ------------------------
void Prediction :: PrintGff ( FILE *OUT, char *seqName )
{
	int offset = PAR.getI ( "Output.offset" );
	int stepid = PAR.getI ( "Output.stepid" );
	int estopt = PAR.getI ( "Sensor.Est.use", 0, PAR.getI ( "EuGene.sloppy" ) );
	int start, end;
	State* featState;
	int incons = 0, cons = 0;

	/* name source feature start end score strand frame */
	for ( int i=0; i<nbGene; i++ )
	{
		for ( int j=0; j<vGene[i]->nbFea(); j++ )
		{
			featState = vGene[i]->vFea[j]->featureState;
			// Print introns ?
//             if ( (PAR.getI("Output.intron") == 0) &&
//                     (state >= IntronF1 && state <= InterGen) ||
//                     (state >= IntronU5F) )
			if ( ( PAR.getI ( "Output.intron" ) == 0 && featState->IsIntronInStartStopRegion() ) ||
			        featState->IsUTRIntron() )
				continue;

			start = vGene[i]->vFea[j]->start;
			end   = vGene[i]->vFea[j]->end;
			if ( estopt ) CheckConsistency ( start, end, featState->GetState(), &cons, &incons );

			fprintf ( OUT, "%s.%d",seqName, ( ( vGene[i]->geneNumber ) *stepid ) +1 );
			if ( vGene[i]->hasvariant )
				fprintf ( OUT, "%c",vGene[i]->GetVariantCode() );
			fprintf ( OUT, ".%d\tEuGene\t%s\t%d\t%d\t%.0f.%.0f\t%c\t", vGene[i]->vFea[j]->number,
			          featState->State2GFFString(), start+offset, end+offset,
			          100.0* ( double ) cons/ ( end-start+1 ), 100.0* ( double ) incons/ ( end-start+1 ),
			          vGene[i]->vFea[j]->strand );
			if ( vGene[i]->vFea[j]->phase == 9 )
				fprintf ( OUT, ".\n" );
			else
				fprintf ( OUT, "%d\n",abs ( vGene[i]->vFea[j]->phase ) );
		}
	}
	//ajoute une ligne vide en fin de fichier
	fprintf ( OUT, "\n" );
}

//SEB ------------------------
//  print prediction (GFF3)
// ------------------------
void Prediction :: PrintGff3 ( std::ofstream& out, char *seqName, char append )
{
	int offset = PAR.getI ( "Output.offset" );
	int stepid = PAR.getI ( "Output.stepid" );
	int estopt = PAR.getI ( "Sensor.Est.use", 0, PAR.getI ( "EuGene.sloppy" ) );
	int start, end;
	State* featState;
	int incons = 0, cons = 0;

	/* name source feature start end score strand frame */
	for ( int i=0; i<nbGene; i++ )
	{
		Gff3Line gene_line ( seqName );
		std::string gene_name = to_string ( seqName ) +'.'+to_string ( ( vGene[i]->geneNumber*stepid+1 ) );
		//je pars du principe que la premiere structure me donne le brin du gene
		// done also on variants to have a corrrect strand on variant mRNA
		if ( vGene[i]->nbFea() > 0 )
			gene_line.setStrand ( vGene[i]->vFea[0]->strand );
		std::string gene_id;
		if ( vGene[i]->isNpcRna )
		{
			gene_id = Sofa::getName ( SOFA_NCRNA ) + ":" + gene_name;
			gene_line.setStart ( vGene[i]->trStart+1+offset );
			gene_line.setEnd ( vGene[i]->trEnd+1+offset );
			gene_line.setType ( SOFA_NCRNA );
			gene_line.setAttribute ( "ID=" + gene_id );
			gene_line.addAttribute ( "Name=" + gene_name );
			gene_line.addAttribute ( "length=" + to_string ( vGene[i]->geneLength ) );
			gene_line.print ( out );
		}
		else
		{
			//  --  --  la ligne du gene  --  --

			gene_id = Sofa::getName ( SOFA_GENE ) + ":" + gene_name;

			if ( append &&  ! vGene[i]->isvariant )
			{
				fprintf ( stderr,"Warning: in append mode, only variant can be printed out gene=%d\n",i );
				// je laisse l'erreur dans .gff1 pour aider au debugage
				continue;
			}

			if ( ! vGene[i]->isvariant ) // only for the optimal prediction
			{
				int tr_min = ( vGene[i]->tuStart ) ?  vGene[i]->tuStart : vGene[i]->trStart;
				int tr_max = ( vGene[i]->tuEnd )   ?  vGene[i]->tuEnd   : vGene[i]->trEnd;

				//la ligne decrivant le gene
				gene_line.setType ( SOFA_GENE );
				gene_line.setStart ( tr_min+1+offset );
				gene_line.setEnd ( tr_max+1+offset );

				//l'identifiant du gene = type_sofa:nom_eugene

				gene_line.setAttribute ( "ID=" + gene_id );
				gene_line.addAttribute ( "Name=" + gene_name );
				gene_line.addAttribute ( "length=" + to_string ( vGene[i]->geneLength ) );
				gene_line.print ( out );
			}
//  --  --  fin de la ligne du gene  --  --

//  --  --  la ligne de l'ARNm  -- <=> basee sur le gene --
			//En fait je me contente d'ecraser ce que je ne veux plus

			gene_line.setStart ( vGene[i]->trStart+1+offset );
			gene_line.setEnd ( vGene[i]->trEnd+1+offset );
			std::string mrna_name = gene_name;
			if ( vGene[i]->hasvariant )
			{
				mrna_name += to_string ( vGene[i]->GetVariantCode() );
			}

			gene_line.setType ( SOFA_MRNA );
			std::string rna_id = Sofa::getName ( SOFA_MRNA ) + ":"+ mrna_name;
			gene_line.setAttribute ( "ID=" + rna_id );
			gene_line.addAttribute ( "Name=" + mrna_name );
			gene_line.addAttribute ( "Parent=" + gene_id );
			gene_line.addAttribute ( "nb_exon=" + to_string ( vGene[i]->exNumber ) );

//----- fin de la ligne du mRNA

			Gff3Line *pre_arn_line, *arn_line;
			std::vector<Gff3Line*> pre_arn_lines, arn_lines;
			std::string fea_name ( gene_name + "." );
			std::string fea_id;
			//postulat : ici on a un seul ARNm par gene // JER je change de postulat.
			int taille_mRNA = 0;

			for ( int j=0; j<vGene[i]->nbFea(); j++ )
			{
				featState = vGene[i]->vFea[j]->featureState;
				// Print introns ?
//             if ( (PAR.getI("Output.intron") == 0) &&
//                     (state >= IntronF1 && state <= InterGen) ||
//                     (state >= IntronU5F) )
				if ( ( PAR.getI ( "Output.intron" ) == 0 && featState->IsIntronInStartStopRegion() ) ||
				        featState->IsUTRIntron() )
					continue;

				start = vGene[i]->vFea[j]->start;
				end   = vGene[i]->vFea[j]->end;
				char code_variant = ( vGene[i]->hasvariant ) ? vGene[i]->GetVariantCode() : 0;
				if ( estopt ) CheckConsistency ( start, end, featState->GetState(), &cons, &incons );

				int type_sofa_cds = featState->getTypeSofa ( true ); // CDS
				int type_sofa     = featState->getTypeSofa ( false ); // exon
				//selon le type on renseigne une ligne existante ou nouvelle
				// d'arn, de cds ou des 2
				switch ( type_sofa )
				{
						//cas d'une region intergenique, ou d'un intron
						//on ne le signale que dans le pre_arn
					case SOFA_INTERGEN  :
					case SOFA_INTRON    :
					case SO_UTR_INTRON  :
						pre_arn_line = fillGff3Line ( type_sofa, start+offset, end+offset,
						                              vGene[i]->vFea[j]->strand,
						                              vGene[i]->vFea[j]->phase );
						//les attributs
						setGff3Attributes ( pre_arn_line, featState, type_sofa, fea_name, j,code_variant, rna_id, false );
						//on l'ajoute au vecteur
						pre_arn_lines.push_back ( pre_arn_line );
						break;
						//cas d'un UTR 3' ou 5'
						//c'est [la fin d']un exon du pre_arn, et un utr de l'arn
					case SOFA_5_UTR :
					case SOFA_3_UTR :
						if ( ( !pre_arn_lines.empty() ) &&
						        ( previousExonMustBeUpdated ( pre_arn_lines.back(), start+offset ) ) )
							//fusion des exons
							pre_arn_lines.back()->setEnd ( end+offset );
						else
						{
							pre_arn_line = fillGff3Line ( SOFA_EXON, start+offset, end+offset,
							                              vGene[i]->vFea[j]->strand,
							                              vGene[i]->vFea[j]->phase );
							//attributs
							setGff3Attributes ( pre_arn_line, featState, SOFA_EXON,
							                    fea_name, j, code_variant, rna_id, false );
							pre_arn_lines.push_back ( pre_arn_line );
						}

						arn_line = fillGff3Line ( type_sofa_cds, start+offset, end+offset,
						                          vGene[i]->vFea[j]->strand,
						                          vGene[i]->vFea[j]->phase );
						//attributs
						setGff3Attributes ( arn_line, featState, type_sofa_cds,
						                    fea_name, j, code_variant, rna_id, true );
						arn_lines.push_back ( arn_line );
						//je mets � jour la taille de l'ARNm
						taille_mRNA += ( end-start+1 );
						break;
						//cas d'un exon
						//il est peut etre la fin du precedent
					case SOFA_EXON :
//               fin = -99999;
						if ( ( !pre_arn_lines.empty() ) &&//c'est le premier
						        ( previousExonMustBeUpdated ( pre_arn_lines.back(), start+offset ) ) )
							//fusion des exons
							pre_arn_lines.back()->setEnd ( end+offset );
						else
						{
							pre_arn_line = fillGff3Line ( SOFA_EXON, start+offset, end+offset,
							                              vGene[i]->vFea[j]->strand,
							                              vGene[i]->vFea[j]->phase );
							pre_arn_lines.push_back ( pre_arn_line );
						}
						//les attributs
//               je le mets ici pour permettre la MAJ d'ontology term.
//               Pour bien faire il faudrait une fonction SetOntologyTerm()
//                 que l'on placerait au niveau de la MAJ
//                 et le setAttributes resterait dans le else
						//si je suis sur un single exon, le code de l'ontology doit etre SOFA_EXON
						if ( featState->IsSnglExon() )
							//je passe nbTracks pour etre sur d'etre hors zone
							setGff3Attributes ( pre_arn_lines.back(), new State ( -1 ), type_sofa, fea_name,
							                    vGene[i]->vFea[j]->number, code_variant, rna_id, false );
						else
							setGff3Attributes ( pre_arn_lines.back(), featState, type_sofa, fea_name,
							                    vGene[i]->vFea[j]->number, code_variant, rna_id, false );

						arn_line = fillGff3Line ( SOFA_CDS, start+offset, end+offset,
						                          vGene[i]->vFea[j]->strand,
						                          vGene[i]->vFea[j]->phase );
						//les attributs
						setGff3Attributes ( arn_line, featState, SOFA_CDS,
						                    fea_name, vGene[i]->vFea[j]->number, code_variant, rna_id, true );
						arn_lines.push_back ( arn_line );
						//je mets � jour la taille de l'ARNm
						taille_mRNA += ( end-start+1 );
						break;
						//par defaut je n'ajoute que dans le pre-arn
					default :
						pre_arn_line = fillGff3Line ( type_sofa, start+offset, end+offset,
						                              vGene[i]->vFea[j]->strand,
						                              vGene[i]->vFea[j]->phase );
						//les attributs
						setGff3Attributes ( pre_arn_line, featState, SOFA_CDS, fea_name, j, code_variant, rna_id, false );
						pre_arn_lines.push_back ( pre_arn_line );
				}//fin switch
				//si on vient d'ajouter une ligne dans les ARNs
				//Je lui ajoute le score
				if ( estopt )
				{
					if ( ( type_sofa==SOFA_5_UTR ) ||
					        ( type_sofa==SOFA_3_UTR ) ||
					        ( type_sofa==SOFA_EXON ) )
					{
						double score = 100.0* ( double ) cons/ ( end-start+1 );
						char buff[6];
						snprintf ( buff, 6, "%.1f", score );
						arn_lines.back()->addAttribute ( "est_cons=" + std::string ( buff ) );
						score  = 100.0* ( double ) incons/ ( end-start+1 );
						snprintf ( buff, 6, "%.1f", score );
						arn_lines.back()->addAttribute ( "est_incons=" + to_string ( buff ) );
					}
				}
			}//fin des structures du gene

			//J'ajoute la taille du mRNA � sa ligne
			gene_line.addAttribute ( "length=" + to_string ( taille_mRNA ) );
			//j'affiche la ligne du mRNA juste avant les UTRs et CDS
			gene_line.print ( out );

			//on affiche les lignes, j'en profite pour liberer les pointeurs
			std::vector<Gff3Line*>::iterator itr;
			//le pre-arn
			for ( itr = pre_arn_lines.begin();
			        itr != pre_arn_lines.end();
			        ++itr )
			{
				( *itr )->print ( out );
				delete ( *itr );
			}
			//l'arn
			for ( itr = arn_lines.begin();
			        itr != arn_lines.end();
			        ++itr )
			{
				( *itr )->print ( out );
				delete ( *itr );
			}
		} // end of else (not ncrna)
		// fin du gene
		out << Gff3Line::endComplexFeature();

	}//fin des genes
}
//SEB
// ----------------------------
//  print prediction (egn long)
// ----------------------------
void Prediction :: PrintEgnL ( FILE *OUT, char *seqName, int a )
{
	int offset = PAR.getI ( "Output.offset" );
	int stepid = PAR.getI ( "Output.stepid" );
	int estopt = PAR.getI ( "Sensor.Est.use", 0, PAR.getI ( "EuGene.sloppy" ) );
	int forward, start, end, don, acc, phase, frame;
	int incons = 0, cons = 0;
	State* featState;

	/* Seq Type S Lend Rend Length Phase Frame Ac Do Pr */
	for ( int i=0; i<nbGene; i++ )
	{
		forward = ( vGene[i]->vFea[0]->strand == '+' ) ? 1 : 0;
		for ( int j=0; j<vGene[i]->nbFea(); j++ )
		{
			featState = vGene[i]->vFea[j]->featureState;

			// Print introns ?
			/*if ( (PAR.getI("Output.intron") == 0) &&
			        (state >= IntronF1 && state <= InterGen) ||
			        (state >= IntronU5F) )*/
			if ( ( PAR.getI ( "Output.intron" ) == 0 && featState->IsIntronInStartStopRegion() ) ||
			        featState->IsUTRIntron() )
				continue;
			phase = vGene[i]->vFea[j]->phase;
			frame = vGene[i]->vFea[j]->frame;
			start = vGene[i]->vFea[j]->start;
			end   = vGene[i]->vFea[j]->end;

			// Frameshift detection: we just concatenate all the parts
			// of the exon in one. The phase will be the "initial" phase.
			while ( ( vGene[i]->vFea[j]->IsCodingExon() ) &&
			        ( j < vGene[i]->nbFea()-1 ) && ( vGene[i]->vFea[j+1]->IsCodingExon() ) )
			{
				j++;
				end = vGene[i]->vFea[j]->end;
			}

			if ( forward )
			{
				don = offset + start - 1;
				acc = offset + end   + 1;
			}
			else
			{
				acc = offset + start - 1;
				don = offset + end   + 1;
			}

			// araset ?
			if ( a ) fprintf ( OUT, "%s ",seqName );

			else
			{
				fprintf ( OUT, "%s.%d",seqName, ( ( vGene[i]->geneNumber ) *stepid ) +1 );
				if ( vGene[i]->hasvariant )
					fprintf ( OUT, "%c",vGene[i]->GetVariantCode() );
				fprintf ( OUT, ".%d\t", vGene[i]->vFea[j]->number );
			}

			fprintf ( OUT, "%s    %c    %7d %7d     %4d  ",
			          featState->State2EGNString(),
			          vGene[i]->vFea[j]->strand,
			          start+offset, end+offset, end - start +1 );

			if ( featState->IsUTR() || featState->IsNpcRna() )
			{
				fprintf ( OUT,"   NA      NA" ); // No phase No frame
				fprintf ( OUT,"      NA      NA " );  // No Don No Acc
			}
			else
			{
				// print the phase if it is defined
				if ( phase == 9 )
					fprintf ( OUT," Unk." );
				else
					fprintf ( OUT,"   %d", phase );
				// print the frame if it is defined and if its is an exon
				if ( frame == 9 )
					fprintf ( OUT,"    Unk." );
				else if ( featState->IsIntron() || frame == 0 )
					fprintf ( OUT,"      NA" );
				else
					fprintf ( OUT,"      %+2d", frame );
				// print donor and acceptor
				fprintf ( OUT," %7d %7d ", don, acc );
			}

			if ( estopt ) CheckConsistency ( start, end, featState->GetState(), &cons, &incons );
			fprintf ( OUT,"  %3.0f.%-3.0f\n",100.0* ( double ) cons/ ( end-start+1 ),
			          100.0* ( double ) incons/ ( end-start+1 ) );
		}
	}
}


// -----------------------------
//  print prediction (egn short)
// -----------------------------
void Prediction :: PrintEgnS ( FILE *OUT )
{
	int offset = PAR.getI ( "Output.offset" );
	int start, end;

	for ( int i=0; i<nbGene; i++ )
	{
		if ( i!=0 )
			fprintf ( OUT,"\n" );
		for ( int j=0; j<vGene[i]->nbFea(); j++ )
		{
			if ( vGene[i]->vFea[j]->IsCodingExon() || vGene[i]->vFea[j]->IsNpcRna() )
			{
				start = offset + vGene[i]->vFea[j]->start;
				end   = offset + vGene[i]->vFea[j]->end;
				fprintf ( OUT,"%d %d ", start, end );
			}
		}
	}
}

// -----------------------------
//  print prediction (egn debug)
// -----------------------------
void Prediction :: PrintEgnD ( FILE *OUT )
{
	DATA Data;

	fprintf ( OUT, "   pos nt  EF1   EF2   EF3   ER1   ER2   ER3    IF    IR    IG   U5F   U5R   U3F   U3R   IUF   IUR   RF   RR FWD tSta  tSto   Sta   Sto   Acc   Don   Ins   Del   tStaNcp   tStopNcp REV tSta  tSto   Sta   Sto   Acc   Don   Ins   Del   tStaNcp   tStopNcp noF tSta  tSto   Sta   Sto   Acc   Don   Ins   Del   tStaNcp   tStopNcp noR tSta  tSto   Sta   Sto   Acc   Don   Ins   Del   tStaNcp   tStopNcp\n" );
	for ( int i=0; i<X->SeqLen; i++ )
	{
		MS->GetInfoAt ( X, i, &Data );
		MS->PrintDataAt ( X, i, &Data, OUT );
	}
	fprintf ( OUT, "   pos nt  EF1   EF2   EF3   ER1   ER2   ER3    IF    IR    IG   U5F   U5R   U3F   U3R   IUF   IUR   RF   RR FWD tSta  tSto   Sta   Sto   Acc   Don   Ins   Del   tStaNcp   tStopNcp REV tSta  tSto   Sta   Sto   Acc   Don   Ins   Del   tStaNcp   tStopNcp noF tSta  tSto   Sta   Sto   Acc   Don   Ins   Del   tStaNcp   tStopNcp noR tSta  tSto   Sta   Sto   Acc   Don   Ins   Del   tStaNcp   tStopNcp\n" );
}

// --------------------------
//  print prediction (Html)
// --------------------------
void Prediction :: PrintHtml ( FILE *OUT, char *seqName )
{
	int start, end;
	int offset     = PAR.getI ( "Output.offset" );
	int stepid     = PAR.getI ( "Output.stepid" );
	char *html_dir = new char[FILENAME_MAX+1];
	strcpy ( html_dir, PAR.getC ( "web_dir" ) );
	State * featState;
	StartHTML ( html_dir, OUT );
	OutputHTMLFileNames ( 1, html_dir, OUT );
	fprintf ( OUT,
	          "         </select>\n"
	          "       </td>\n"
	          "       <td width=\"100\" "
	          "align=\"center\"\n"
	          "           bgcolor=\"#c0dbe2\">\n"
	          "         <img src=\"%s"
	          "/Images/next.jpg\"\n"
	          "              onclick=\"next();\"\n"
	          "              title=\"Next\" "
	          "align=\"middle\">\n"
	          "         &nbsp; &nbsp; &nbsp;\n"
	          "         <img src=\"%s"
	          "/Images/last.jpg\"\n"
	          "              onclick=\"last();\"\n"
	          "              title=\"Jump to end\" "
	          "align=\"middle\">\n", html_dir, html_dir );
	fprintf ( OUT,
	          "       </td>\n"
	          "           </tr>\n"
	          "           <tr>\n"
	          "       <td bgcolor=\"#9bc7d0\"></td>\n"
	          "       <td align=\"center\" "
	          "bgcolor=\"white\">\n"
	          "         <div class=\"conteneur\">\n"
	          "           <div class=\"frame\">\n"
	          "             <table width=\"100%%\" "
	          "cellpadding=\"1\" cellspacing=\"1\">\n"
	          "           <tr>\n"
	          "             <td>\n"
	          "               <table "
	          "width=\"100%%\""
	          " cellpadding=\"2\" cellspacing=\"1\">\n" );
	fprintf ( OUT,
	          "<tr align=\"center\" class=\"fonce\">\n"
	          " <td><font color=\"white\" face=\"monospace\">\n"
	          "  <b>SeqName_GeneN.ExonN</b></font></td>\n"
	          " <td><font color=\"white\" face=\"monospace\">\n"
	          "  <b>Source</b></font></td>\n"
	          " <td><font color=\"white\" face=\"monospace\">\n"
	          "  <b>Feature</b></font></td>\n"
	          " <td><font color=\"white\" face=\"monospace\">\n"
	          "  <b>Start</b></font></td>\n"
	          " <td><font color=\"white\" face=\"monospace\">\n"
	          "  <b>End</b></font></td>\n"
	          " <td><font color=\"white\" face=\"monospace\">\n"
	          "  <b>Score</b></font></td>\n"
	          " <td><font color=\"white\" face=\"monospace\">\n"
	          "  <b>Strand</b></font></td>\n"
	          " <td><font color=\"white\" face=\"monospace\">\n"
	          "  <b>Frame</b></font></td>\n"
	          "</tr>\n" );
	if ( nbGene == 0 )
		fprintf ( OUT,
		          "<tr class=\"A0\">\n <td colspan=\"8\" "
		          "align=\"center\">\n  <font face=\"monospace\">\n  "
		          "<b>No exons/genes predicted in your submitted "
		          "sequence !</b></td>\n</tr>\n" );

	int lc = 0;
	for ( int i=0; i<nbGene; i++ )
	{
		for ( int j=0; j<vGene[i]->nbFea(); j++ )
		{
			featState = vGene[i]->vFea[j]->featureState;
			// Print introns ?
//             if ( (PAR.getI("Output.intron") == 0) &&
//                     (state >= IntronF1 && state <= InterGen) ||
//                     (state >= IntronU5F) )
			if ( ( PAR.getI ( "Output.intron" ) == 0 && featState->IsIntronInStartStopRegion() ) ||
			        featState->IsUTRIntron() )
				continue;

			start = offset + vGene[i]->vFea[j]->start;
			end   = offset + vGene[i]->vFea[j]->end;
			fprintf ( OUT,
			          "<tr class=\"A%d\">\n"
			          " <td align=\"center\">\n"
			          "  <font face=\"monospace\">%s_%d.%d</font></td>\n"
			          "<td align=\"center\">\n"
			          "  <font face=\"monospace\">EuGene_%s</font></td>\n",
			          lc++%2, seqName, ( i*stepid ) +1, vGene[i]->vFea[j]->number,
			          PAR.getC ( "EuGene.organism" ) );
			fprintf ( OUT,
			          " <td align=\"center\">\n  <font face="
			          "\"monospace\">%s</font></td>\n", featState->State2GFFString() );
			fprintf ( OUT,
			          " <td align=\"right\">\n"
			          "  <font face=\"monospace\">%d</font></td>\n"
			          " <td align=\"right\">\n"
			          "  <font face=\"monospace\">%d</font></td>\n"
			          " <td align=\"center\">\n"
			          "  <font face=\"monospace\">."
			          "</font></td>\n"
			          " <td align=\"center\">\n"
			          "  <font face=\"monospace\">%c</font></td>\n"
			          " <td align=\"center\">\n"
			          "  <font face=\"monospace\">",
			          start, end, vGene[i]->vFea[j]->strand );

			if ( vGene[i]->vFea[j]->phase == 9 )
				fprintf ( OUT, "." );
			else
				fprintf ( OUT, "%d",abs ( vGene[i]->vFea[j]->phase ) );
			fprintf ( OUT, "</font></td>\n</tr>\n" );
		}
	}
	EndHTML ( OUT );
	delete [] html_dir;
}


// ----------------------------------------------------------------
// If the feature checked is a UTR, then the UTR will be checked for
// support (bool returned) and start/end positions may be trimmed to
// match with the first experimental (EST) evidence. All EST evidence
// is used.
// ----------------------------------------------------------------
bool Prediction :: UTRCheckAndTrim ( int* debut, int* fin, int etat )
{
	if ( ( etat == UTR5F ) || ( etat == UTR3R ) )
		for ( int k = *debut; k <= *fin; k++ )
			if ( GetESTMatch ( k-1 ) & Hit )
			{
				*debut = k;
				return true;
			}

	if ( ( etat == UTR3F ) || ( etat == UTR5R ) )
		for ( int k = *fin; k >= *debut; k-- )
			if ( GetESTMatch ( k-1 ) & Hit )
			{
				*fin = k;
				return true;
			}

	return false;
}
// ----------------------------------------------------------------
// Verif coherence EST: calcul le nombre de nuc. coherents et
// incoherents avec les match est debut/fin/etat: debut et fin de la
// seq. dont l'etat est etat cons/incons: retour des valeurs.  Seuls
// les EST informant ESTmatch (non filtrs) sont pris en compte.
// ----------------------------------------------------------------
void Prediction :: CheckConsistency ( int debut, int fin, int etat,
                                      int* cons, int* incons )
{
	int i, con = 0, inc = 0;
	DATA dTMP;

	// les valeurs qui sont coherentes avec chaque etat
	const unsigned char Consistent[NbTracks] =
	{
		HitForward|MForward,    HitForward|MForward,    HitForward|MForward,
		HitReverse|MReverse,    HitReverse|MReverse,    HitReverse|MReverse,
		HitForward|MForward,    HitForward|MForward,    HitForward|MForward,
		HitReverse|MReverse,    HitReverse|MReverse,    HitReverse|MReverse,
		HitForward|MForward,    HitForward|MForward,    HitForward|MForward,
		HitReverse|MReverse,    HitReverse|MReverse,    HitReverse|MReverse,
		HitForward|MForward,    HitForward|MForward,    HitForward|MForward,
		HitReverse|MReverse,    HitReverse|MReverse,    HitReverse|MReverse,
		GapForward|MForward,    GapForward|MForward,    GapForward|MForward,
		GapReverse|MReverse,    GapReverse|MReverse,    GapReverse|MReverse,
		GapForward|MForward,    GapForward|MForward,    GapForward|MForward,
		GapReverse|MReverse,    GapReverse|MReverse,    GapReverse|MReverse,
		0,
		HitForward|MForward, HitForward|MForward,
		HitReverse|MReverse, HitReverse|MReverse,
		GapForward|MForward, GapReverse|MReverse,
		GapForward|MForward, GapReverse|MReverse
	};

	const unsigned char MaskConsistent[NbTracks] =
	{
		Hit|Margin,    Hit|Margin,    Hit|Margin,
		Hit|Margin,    Hit|Margin,    Hit|Margin,
		Hit|Margin,    Hit|Margin,    Hit|Margin,
		Hit|Margin,    Hit|Margin,    Hit|Margin,
		Hit|Margin,    Hit|Margin,    Hit|Margin,
		Hit|Margin,    Hit|Margin,    Hit|Margin,
		Hit|Margin,    Hit|Margin,    Hit|Margin,
		Hit|Margin,    Hit|Margin,    Hit|Margin,
		Gap|Margin,    Gap|Margin,    Gap|Margin,
		Gap|Margin,    Gap|Margin,    Gap|Margin,
		Gap|Margin,    Gap|Margin,    Gap|Margin,
		Gap|Margin,    Gap|Margin,    Gap|Margin,
		0,
		Hit|Margin,    Hit|Margin,
		Hit|Margin,    Hit|Margin,
		Gap|Margin,    Gap|Margin,
		Gap|Margin,    Gap|Margin
	};

	for ( i=Max ( debut-1,0 ); i<fin; i++ )
	{

		MS->GetInfoSpAt ( Type_Content, X, i, &dTMP );

		// y a t'il de l'info
		if ( dTMP.EstMatch )
		{
			// y a t'il une info incoherente avec l'etat
			if ( dTMP.EstMatch & ~MaskConsistent[etat] )
				inc++;
			else if ( dTMP.EstMatch & Consistent[etat] )
				con++;
			else
				inc++;
		}
	}
	*cons = con;
	*incons = inc;
}

//---------------------------------------------------
// PrintHtml : -ph print the begin of the HTML output
//---------------------------------------------------
void Prediction :: StartHTML ( char* html_dir, FILE *OUT )
{
	char *d = new char[MAX_LINE];
	GetStrDate ( d );

	fprintf ( OUT,
	          "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n"
	          "<html>\n"
	          "  <head>\n"
	          "    <title>EuGene : Prediction</title>\n"
	          "    <link rel=\"STYLESHEET\" type=\"text/css\" "
	          "href=\"%s/Style/eugene.css\">\n\n"
	          "    <script language=\"JavaScript1.2\" "
	          "src=\"%s/Javascripts/diap.js\">\n"
	          "    </script>\n"
	          "  </head>\n\n"
	          "  <body leftmargin=\"0\" topmargin=\"0\" marginwidth=\"0\" "
	          "marginheight=\"0\"\n"
	          " bgcolor=\"white\">\n"
	          "    <table height=\"100%%\" width=\"100%%\" "
	          "CELLSPACING=\"0\" CELLPADDING=\"0\"\n"
	          "    border=\"0\" rows=\"2\">\n"
	          " <tr height=\"140\">\n"
	          "   <td BACKGROUND=\"%s/Images/top.jpg\" "
	          "colspan=\"2\" height=\"140\"\n"
	          "       valign=\"top\">&nbsp;</td>\n"
	          " </tr>\n"
	          " <tr>\n"
	          "   <td width=\"145\" valign=\"top\" align=\"left\"\n"
	          "       BACKGROUND=\"%s/Images/left.jpg\">\n"
	          "     <img src=\"%s/Images/left.jpg\" "
	          "width=\"145\" height=\"10\"></td>\n"
	          "   <td width=\"100%%\" valign=\"top\">\n\n"
	          "     <!-- DEBUT PAGE... -->\n"
	          "     <table CELLSPACING=\"0\" CELLPADDING=\"6\" "
	          "border=\"0\" width=\"100%%\">\n"
	          "   <tr>\n"
	          "                  <td colspan=\"2\"> </td>\n"
	          "                </tr>\n"
	          "                <tr>\n"
	          "                  <td colspan=\"2\"> "
	          "<img src=\"%s/Images/euGpred_on.jpg\">\n",
	          html_dir, html_dir, html_dir,
	          html_dir, html_dir, html_dir );

	fprintf ( OUT,
	          "<br><font size=-1> - "
	          "<a href=\"http://www.inra.fr/bia/T/EuGene/index.html\">EuGene</a> "
	          "%s for %s, %s -",
	          PAR.getC ( "EuGene.version" ), PAR.getC ( "EuGene.organism" ), d );

	if ( rindex ( PAR.getC ( "Output.format" ), 'g' ) != NULL )
		fprintf ( OUT," <a href=\"./%s.gff\">Gff output</a> -", PAR.getC ( "prefixName" ) );

	if ( PAR.getI ( "Sensor.Est.use" )    && PAR.getI ( "Est.PostProcess" )  ||
	        PAR.getI ( "Sensor.BlastX.use" ) && PAR.getI ( "BlastX.PostProcess" ) )
		fprintf ( OUT," <a href=\"./%s.misc_info\">Evidences</a> -", PAR.getC ( "prefixName" ) );

	fprintf ( OUT,
	          "<br></font></td>\n"
	          "                </tr>\n"
	          "                <tr>\n"
	          "      <td>\n"
	          "       <form name=\"formname\">\n"
	          "         <div align=\"center\">\n"
	          "     <center>\n"
	          "       <table cellspacing=\"3\" cellpadding=\"5\" "
	          "border=\"0\" bgcolor=\"#c0dbe2\">\n"
	          "           <tr>\n"
	          "       <td colspan=\"3\" align=\"center\" "
	          "bgcolor=\"white\">\n",
	          PAR.getC ( "EuGene.version" ), PAR.getC ( "EuGene.organism" ), d );

	delete [] d;
}

//-------------------------------------------------
// PrintHtml : -ph print the end of the HTML output
//-------------------------------------------------
void Prediction :: EndHTML ( FILE *OUT )
{
	fprintf ( OUT,
	          "                                              </table>\n"
	          "             </td>\n"
	          "           </tr>\n"
	          "             </table>\n"
	          "           </div>\n"
	          "         </div>\n"
	          "       </td>\n"
	          "       <td  bgcolor=\"#9bc7d0\"></td>\n"
	          "           </tr>\n"
	          "       </table>\n"
	          "     </center>\n"
	          "         </div>\n"
	          "       </form>\n"
	          "     </td>\n"
	          "   </tr>\n"
	          "     </table>\n"
	          "   </td>\n"
	          " </tr>\n"
	          "    </table>\n"
	          "    <!-- FIN PAGE -->\n"
	          "  </body>\n"
	          "</html>\n" );
}

// ------------------------
//  Print gene info.
// ------------------------
void Prediction :: PrintGeneInfo ( FILE* F )
{
	fprintf ( F, "#=========================================#\n" );
	fprintf ( F, "#=           Gene informations           =#\n" );
	fprintf ( F, "#=========================================#\n" );
	for ( int i=0; i<nbGene; i++ )
	{
		vGene[i]->PrintInfo ( F, i+1, seqName );
	}
}

// ------------------------------------
//  Find the index of the 1st gene maximally overlapping a segment.
//  NULL if not found. Assumes a sorted prediction.
// ------------------------------------
Gene *Prediction :: FindGene ( int start, int end )
{
	int left, right;
	int overlap = 0;
	int idx = -1;

	for ( int i=0; i<nbGene; i++ )
	{
		left =  Max ( start, vGene[i]->trStart );
		right = Min ( end, vGene[i]->trEnd );

		if ( ( ( right - left ) <= 0 ) && overlap != 0 ) return vGene[idx]; // cannot be -1
		if ( ( right - left ) > overlap )
		{
			idx = i;
			overlap = right-end;
		}
	}
	if ( idx == -1 ) return NULL;
	else return vGene[idx];
}

// ------------------------
//  plot a prediction.
// ------------------------
void Prediction :: PlotPred ()
{
	const int predWidth = 2;
	int start, end;
	short int featureFrame;
	int endBack = 0;

	for ( int i=0; i<nbGene; i++ )
	{
		for ( int j=0; j<vGene[i]->nbFea(); j++ )
		{
			featureFrame = vGene[i]->vFea[j]->frame;
			start        = vGene[i]->vFea[j]->start;
			end          = vGene[i]->vFea[j]->end;
			// Plot interg
			for ( int k=endBack+1; k<start; k++ )
				PlotBarI ( k, 0, 0.4, predWidth, 4 );
			endBack = end;

			for ( int k=start; k<end; k++ )
			{
				if ( vGene[i]->vFea[j]->IsNpcRna() ) // change the color if its a npcrna
					PlotBarI ( k, featureFrame, 0.4, predWidth, 11 );
				else
					PlotBarI ( k, featureFrame, 0.4, predWidth, 1 );
			}
		}
	}
}

// -------------------------------
//  get the state for a given pos.
// -------------------------------
State* Prediction :: GetStateAtPos ( int pos )
{
	for ( int i=0; i<nbGene; i++ )
	{
		if ( pos < vGene[i]->vFea[0]->start ) return new State ( InterGen ); // IG
		else
			if ( pos <= vGene[i]->vFea[vGene[i]->nbFea()-1]->end )
				for ( int j=0; j<vGene[i]->nbFea(); j++ )
					if ( pos <= vGene[i]->vFea[j]->end )
						return vGene[i]->vFea[j]->featureState;
	}
	return new State ( -1 );
}


// ------------------------
//  isStart.
// True if a start codon starts at the position (p+1) in the strand forward
// True if a stard codon starts at the position p in the reverse strand
// ------------------------
const char* Prediction :: IsStart ( int p )
{
	State* currentState = GetStateAtPos ( p );
	State* nextState    = GetStateAtPos ( p+1 );

	// If the next position is coding exon on strand forward AND that position is IG or UTR or UTR Intron
	if ( nextState->IsForwardCodingExon() && ( !currentState->InStartStopRegion() ) )
		return "True";

	if ( currentState->IsReverseCodingExon()  &&
	        ( ( !nextState->InStartStopRegion() ) || ( p+1 >= vGene[0]->cdsEnd + 1 ) ) )
		return "True";

	return "False";
}

// ------------------------
//  isStop.
// True if a stop codon stops at the position p on strand forward
// True if a stop codon stops at the position p+1 on reverse strand
// ------------------------
const char* Prediction :: IsStop ( int p )
{
	State *currentState = GetStateAtPos ( p );
	State *nextState    = GetStateAtPos ( p+1 );

	if ( currentState->IsDefined()  &&  currentState->IsForwardCodingExon()  &&
	        ( ( !nextState->InStartStopRegion() ) || ( p+1 >= vGene[0]->cdsEnd + 1 ) ) )
		return "True";

	if ( ( !currentState->InStartStopRegion() ) && nextState->IsReverseCodingExon() )
		return "True";
	return "False";
}

// ------------------------
//  isDon.
// ------------------------
const char* Prediction :: IsDon ( int p )
{
	State *currentState = GetStateAtPos ( p );
	State *nextState    = GetStateAtPos ( p+1 );

	if ( currentState->IsForwardCodingExon() &&  nextState->GetState() == IntronF1 )
		return "True";
	if ( nextState->IsReverseCodingExon() &&  currentState->GetState() == IntronR1 )
		return "True";
	return "False";
}

// ------------------------
//  isAcc.
// ------------------------
const char* Prediction :: IsAcc ( int p )
{
	State *previousState = GetStateAtPos ( p );
	State *currentState  = GetStateAtPos ( p+1 );

	if ( currentState->IsForwardCodingExon() &&  previousState->GetState() == IntronF1 )
		return "True";
	if ( previousState->IsReverseCodingExon()  &&  currentState->GetState() == IntronR1 )
		return "True";
	return "False";
}

// ------------------------
// IsState: Is a nucleotide in a given state ?
// BE CAREFUL:
// this method could be used only for a prediction with one complete gene.
//         the prediction could be the representation of an external gff annotation
//         that could not specify the UTR. In this case, the state from 0 to the first exon
//         in the annotation (the first element of vPos and vState)
//         is InterGen and the state after the last exon
//         in the annotation (the last element of vPos and vState)
//         is not set (getStateForPos returns -1).
//         To identify Start Reverse and Stop Forward, an other condition is used
//         pos == vPos[0] that compares the position with
//                             in Forward the end of the last exon of the gene (vPos[0])
//                             in Reverse the start of the first exon of the gene (vPos[0])
// ------------------------
bool Prediction :: IsState ( DATA::SigType sig_type, int pos, char strand )
{
	bool is_state = false;
	bool bad_strand = false;
	//int state, nState, pState;
	State *currentState;
	State *nextState;
	State *previousState;
	if ( sig_type == DATA::Start )
	{
		currentState = GetStateAtPos ( pos );
		nextState    = GetStateAtPos ( pos+1 );

		if ( strand == '+' )
		{
			// Intergenic then coding exon
			if ( ( currentState->IsIntergenic() ) && ( nextState->IsForwardCodingExon() ) )
				is_state = true;
		}
		else if ( strand == '-' )
		{
			if ( currentState->IsReverseCodingExon() &&
			        ( nextState->GetState() == UTR5R || pos == vGene[0]->cdsEnd + 1 ) )
				is_state = true;
		}
		else bad_strand = true;

	}
	else if ( sig_type == DATA::Stop )
	{
		currentState = GetStateAtPos ( pos );
		nextState    = GetStateAtPos ( pos+1 );

		if ( strand == '+' )
		{
			if ( currentState->IsForwardCodingExon() &&
			        ( ( nextState->GetState() == UTR3F ) || ( pos == vGene[0]->cdsEnd + 1 ) ) )
				is_state = true;
		}
		else if ( strand == '-' )
		{
			if ( nextState->IsReverseCodingExon() &&
			        ( currentState->GetState() == UTR3R || currentState->IsIntergenic() ) )
				is_state = true;
		}
		else bad_strand = true;

	}
	else if ( sig_type == DATA::Acc )
	{
		previousState = GetStateAtPos ( pos );
		currentState  = GetStateAtPos ( pos+1 );

		if ( strand == '+' )
		{
			// If the the current position is in an exon and the previous position is in an intron
			if ( currentState->IsForwardCodingExon() && previousState->IsForwardIntron() ) is_state = true;
		}
		else if ( strand == '-' )
		{
			// If the current position (pos+1) is in an intron and the previous is in coding on the reverse strand
			if ( previousState->IsReverseCodingExon() && currentState->IsReverseIntron() ) is_state = true;
		}
		else bad_strand = true;

	}
	else if ( sig_type == DATA::Don )
	{
		currentState = GetStateAtPos ( pos );
		nextState    = GetStateAtPos ( pos+1 );

		if ( strand == '+' )
		{
			// If the current position is coding and the next position is an intron on the forward strand
			if ( currentState->IsForwardCodingExon() && nextState->IsForwardIntron() ) is_state = true;
		}
		else if ( strand == '-' )
		{
			// If the current position is an intron and the next position an exon in the reverse strand
			if ( currentState->IsReverseIntron() && nextState->IsReverseCodingExon() ) is_state = true;
		}
		else bad_strand = true;

	}
	else
	{
		std::cerr<<"ERROR: bad state "<<sig_type<<" given in argument in Prediction::IsState.\n";
		exit ( 2 );
	}

	if ( bad_strand )
	{
		std::cerr<<"ERROR: bad strand "<<strand<<"  given in argument in Prediction::IsState.\n";
		exit ( 2 );
	}

	return is_state;
}


//SEB--  --  --  --  --  --  --
//  --  --  delocalisation pour eclaircir le code
Gff3Line*
Prediction::fillGff3Line ( int type_sofa, int start, int end,
                           char strand, int phasegff )
{
	Gff3Line* line = new Gff3Line ( seqName );
	line->setType ( type_sofa );
	line->setStart ( start );
	line->setEnd ( end );
	line->setStrand ( strand );
	if ( phasegff != 9 )
		line->setPhase ( abs ( phasegff ) );
	return line;
}

bool
Prediction::previousExonMustBeUpdated ( Gff3Line* line, int start )
{
	if ( !line )
		std::cerr << "bad previous !" << std::endl;
//     return false;
	if ( line->getType() ==SOFA_EXON )
		return ( ( line->getEnd() +1 ) == start );
	return false;
}

// -----------------------------------------
// Fill the 'line' attribute.
// -----------------------------------------
void Prediction::setGff3Attributes ( Gff3Line* line, State* featState,
                                     int type_sofa, std::string fea_name,
                                     int j, char code_variant, std::string gene_id, bool coding )
{
	std::string fea_id = Sofa::getName ( type_sofa ) + ":"
	                     + fea_name + to_string ( j );
	if ( code_variant )
		fea_id += to_string ( code_variant );
	line->setAttribute ( "ID=" + fea_id );
	line->addAttribute ( "Parent="+gene_id );
	line->addAttribute ( "Ontology_term="
	                     + Sofa::codeToString ( featState->getTypeSofa ( coding, false ) ) );
}


// ------------------------
//  Print: print on stdout each feature of each gene of the prediction
// ------------------------
void Prediction::Print()
{
	int state, start, end, phase, frame;
	char strand;
	for ( int i=0; i<nbGene; i++ )
	{
		printf ( "Gene %d :\n", ( i+1 ) );

		for ( int j=0; j<vGene[i]->nbFea(); j++ )
		{
			state = vGene[i]->vFea[j]->featureState->GetState();
			start = vGene[i]->vFea[j]->start;
			end   = vGene[i]->vFea[j]->end;
			phase = vGene[i]->vFea[j]->phase;
			frame = vGene[i]->vFea[j]->frame;
			strand = vGene[i]->vFea[j]->strand;
			printf ( "  Feat %d : state %d=%s, start %d, end %d, phase %d, frame %d strand %c\n", j, state, vGene[i]->vFea[j]->featureState->State2EGNString(), start, end, phase, frame, strand );
		}
	}
}

// ------------------------
// Eval the prediction in comparaison with the reference in gene/exon/nucleotide level
// TP: true positive, P* predicted *, R* real/reference *, N nuccleotide, E exon, G gene
// Return a vector [TPg, PGnb, RGnb, TPe, PEnb, REnb, TPn, PNnb, RNnb]
// offset is the size of the region added before/after the reference
// ------------------------
std::vector<int> Prediction :: Eval ( Prediction * ref, int offset )
{
	std::vector<int> vEval;    //[TPg, PGnb, RGnb, TPe, PEnb, REnb, TPn, PNnb, RNnb]
	std::vector<int> vExonEval;
	// Compute positions of the region around the reference +/- offset
	int regionStart, regionStop;

	if ( ref->nbGene > 0 )
	{
		int refCdsStart = ref->vGene[0]->cdsStart;           // cds start of the first gene
		int refCdsEnd   = ref->vGene[ref->nbGene-1]->cdsEnd; // cds end of the last gene
		regionStart     = max ( 1, refCdsStart - offset );
		regionStop      = min ( this->X->SeqLen, refCdsEnd + offset );
	}
	else
	{
		regionStart = 1;
		regionStop  = this->X->SeqLen;
	}

	// Compute the nb of TP genes, of predicted genes in the region and of ref genes
	vEval     = this->EvalGene ( ref, regionStart, regionStop );
	// Compute the nb of TP exon/nt, of predicted exon/nt in the region and of ref exon/nt
	vExonEval = this->EvalExon ( ref, regionStart, regionStop );
	vEval.insert ( vEval.end(), vExonEval.begin(), vExonEval.end() );
	return vEval;
}

// ------------------------
//  Return true if the CDS of the gene overlaps the CDS of the gene g
// ------------------------
bool Gene :: Overlap ( const Gene& g )
{
	return ( ( this->cdsStart <= g.cdsEnd ) && ( this->cdsEnd >= g.cdsStart ) );
}

// ------------------------
// Eval the predicted genes in comparaison with the reference genes
// Return a vector of 3 int:
//<TP gene nb, the nb of predicted genes in the region [start-end], the nb of real genes>
// NOTE: THE NUMBER OF REAL GENES IS EQUAL TO THE NUMBER OF GENES OF THE REFERENCE BECAUSE WE SUPPOSE THAT REF IS COMPLETELY INCLUDED BETWEEN START AND END POSITIONS
// ------------------------
std::vector<int> Prediction :: EvalGene ( Prediction* ref, int start, int end )
{
	std::vector<int> vResult ( 3,0 ); // [TPg, PGnb, RGnb]

	Gene* predGene; // Predicted gene
	Gene* refGene;  // Reference gene
	int iPredGene  = 0; // index used to scan predicted genes
	int iRefGene   = 0; // index used to scan genes of the reference
	int TPGene     = 0; // number of true positive genes
	int realGeneNb = ref->nbGene;  // number of real genes
	// Get the predicted genes in the region of reference
	std::vector<Gene*> vPredGenes = this->GetGenes ( start, end );
	int predGeneNb = vPredGenes.size(); // nb of predicted genes between start and end

	while ( ( iRefGene < realGeneNb ) && ( iPredGene < predGeneNb ) )
	{
		predGene = vPredGenes[iPredGene];
		refGene  = ref->vGene[iRefGene];
		if ( !predGene->Overlap ( *refGene ) )
		{
			if ( predGene->cdsStart < refGene->cdsStart ) iPredGene++;
			else                                        iRefGene++;
		}
		else // overlapping genes
		{
			if ( *predGene == *refGene ) // if the gene are identical
			{
				TPGene += 1; // gene TP ++
				iPredGene++;
				iRefGene++;
			}
			else // overlapping and different genes
			{
				// if the next predicted gene overlaps the reference, go to the next predicted gene
				if ( ( ( iPredGene+1 ) < predGeneNb ) &&
				        ( vPredGenes[iPredGene+1]->Overlap ( *refGene ) ) )
				{
					iPredGene++; // next predicted gene
				}
				// if the next ref gene overlaps the current pred gene, move to the next ref gene
				else if ( ( ( iRefGene+1 ) < realGeneNb ) &&
				          ( ref->vGene[iRefGene+1]->Overlap ( *predGene ) ) )
				{
					iRefGene++; // next reference gene
				}
				else
				{
					iPredGene++;
					iRefGene++;
				}
			}
		}
	}
	vResult[0] = TPGene;
	vResult[1] = predGeneNb;
	vResult[2] = realGeneNb;
	return vResult;
}

// ------------------------
// Eval the predicted exons in comparaison with the reference exons
// Return a vector of 6 int:
//<TP exons nb, the nb of predicted exons in the region [start-end], the nb of real exons,
// TP nucleotides nb, the nb of nt predicted as coding in the region, the nb of nt really coding>
// NOTE: THE NUMBER OF REAL EXONS IS EQUAL TO THE NUMBER OF EXONS OF THE REFERENCE BECAUSE WE SUPPOSE THAT REF IS COMPLETELY INCLUDED BETWEEN START AND END POSITIONS
// ------------------------
std::vector<int>  Prediction :: EvalExon ( Prediction * ref, int start, int end )
{
	std::vector<int> vResult ( 6,0 ); // //[TPe, PEnb, REnb, TPn, PNnb, RNnb]

	int TPexons  = 0; // number of true positive exons
	int TPnt     = 0; // number of true positive nucleotides
	int PredNtNb = 0; // number of coding nt between start and end
	int overlapStart, overlapEnd;

	int RealExonNb = ref->GetExonNumber();
	int RealNtNb   = ref->GetExonLength();
	std::vector<Feature*> vPredExons = this->GetExons ( start, end );
	std::vector<Feature*> vRefExons  = ref->GetExons ( start, end );
	int PredExonNb = vPredExons.size(); // number of exons predicted between start and stop

	// If some exons were predicted in this region
	if ( PredExonNb > 0 && vRefExons.size() > 0 )
	{
		int iPredExon = 0;
		int iRefExon  = 0;
		Feature* refExon;  // first real exon
		Feature* predExon; // first predicted exon

		// Compute the number of predicted coding nt (PredNtNb)
		// A Ameliorer!!! EK
		for ( int i = 0; i < PredExonNb; i++ )
		{
			predExon = vPredExons[i];
			PredNtNb += min ( predExon->end, end ) - max ( predExon->start, start ) + 1;
		}

		while ( iPredExon < PredExonNb && iRefExon < RealExonNb )
		{
			refExon  = vRefExons[iRefExon];
			predExon = vPredExons[iPredExon];

			// if the predictions are not overlapping
			if ( !predExon->Overlap ( *refExon ) )
			{
				if ( predExon->start < refExon->start )  iPredExon++;
				else                                   iRefExon++;
			}
			else // prediction overlapping
			{
				// If the exons are identical, inc TP exon
				if ( ( predExon->start == refExon->start ) && ( predExon->end == refExon->end ) )
				{
					TPexons += 1;
					iRefExon++;// next real gene
					iPredExon++; // next predicted exon
				}
				else
				{
					// if the next pred exon overlaps the current reference, go to the next pred exon
					if ( ( ( iPredExon+1 ) < PredExonNb ) &&
					        ( vPredExons[iPredExon+1]->Overlap ( *refExon ) ) )
					{
						iPredExon++; // next prediction
					}
					// if the next ref exon overlaps the current pred exon, move to the next ref exon
					else if ( ( ( iRefExon+1 ) < RealExonNb ) &&
					          ( vRefExons[iRefExon+1]->Overlap ( *predExon ) ) )
					{
						iRefExon++;
					}
					else
					{
						iPredExon++;
						iRefExon++;
					}
				}

				// Compute the number of overlapping nucleotides and inc nucleotide TP
				overlapStart = max ( predExon->start, refExon->start );
				overlapEnd   = min ( predExon->end, refExon->end );
				TPnt += overlapEnd - overlapStart + 1; // nucleotide TP += overlap length
			}
		}
	}

	vResult[0] = TPexons;
	vResult[1] = PredExonNb;
	vResult[2] = RealExonNb;
	vResult[3] = TPnt;
	vResult[4] = PredNtNb;
	vResult[5] = RealNtNb;

	return vResult;
}

// ------------------------
//  Extract the exons including or overlapping the region [ begin,end].
// ------------------------
std::vector<Feature*> Prediction :: GetExons ( int begin, int end )
{
	std::vector<Feature *> vSelectedExons; // vector of exons overlapping the region
	Feature* feat;

	if ( begin > end )
	{
		int saveBegin = begin;
		begin         = end;
		end           = saveBegin;
	}

	// for each gene
	for ( int i=0; i < this->nbGene; i++ )
	{
		// for each feature of the gene, check its an exon overlapping the region
		for ( int j = 0; j < this->vGene[i]->nbFea(); j++ )
		{
			feat = this->vGene[i]->vFea[j];
			if ( ( feat->IsCodingExon() == true ) && ( feat->start <= end ) && ( feat->end >= begin ) )
			{
				vSelectedExons.push_back ( feat );
			}
		}
	}
	return vSelectedExons;
}

// ------------------------
//  Extract the genes which are completly or partially included between begin and end positions.
// ------------------------
std::vector<Gene*> Prediction :: GetGenes ( int begin, int end )
{
	std::vector<Gene *> vSelectedGenes; // vector of genes including or overlapping in the region

	if ( begin > end )
	{
		int saveBegin = begin;
		begin         = end;
		end           = saveBegin;
	}

	// for each gene of the prediction
	for ( int i=0; i < this->nbGene; i++ )
	{
		if ( ( this->vGene[i]->cdsStart <= end ) &&
		        ( this->vGene[i]->cdsEnd   >= begin ) )
			vSelectedGenes.push_back ( vGene[i] );
	}

	return vSelectedGenes;
}

// ------------------------
// Compute the number of exons in all the predicted genes
// ------------------------
int Prediction:: GetExonNumber()
{
	int ExonNumber = 0;
	// for each gene, get the exon number
	for ( int i = 0; i < nbGene; i++ )
	{
		ExonNumber += vGene[i]->GetExonNumber();
	}
	return ExonNumber;
}

// ------------------------
// Compute the length of all the exons
// ------------------------
int Prediction :: GetExonLength()
{
	int codingNtNumber = 0;

	// for each gene, get the coding nt number
	for ( int i = 0; i < this->nbGene; i++ )
	{
		codingNtNumber += this->vGene[i]->exLength;
	}

	return codingNtNumber;
}
