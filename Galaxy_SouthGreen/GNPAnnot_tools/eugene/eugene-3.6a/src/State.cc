// ------------------------------------------------------------------
// Copyright (C) 2009 INRA <eugene@ossau.toulouse.inra.fr>
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
// $Id
// ------------------------------------------------------------------
// File:     State.h
// Contents: Class State

#include "State.h"

/*************************************************************
 **                      State                              **
 *************************************************************/

// ----------------------
//  Default constructor.
// ----------------------
State::State()
{
    this->state = -1;
}

// ----------------------
//  Constructor.with a parameter: the state in char format
// ----------------------
State::State(char state)
{
    this->state = state;
}

// ----------------------
//  Return the char value of the state
// ----------------------
char State::GetState()
{
    return this->state;
}

// ----------------------
//  Return true if it is intergenic
// ----------------------
bool State::IsIntergenic(void)
{
    if (State2Status[this->state] == UNTRANSCRIBED)
    {
        return true;
    }
    return false;
}


// ----------------------
//  Return true if it is an intron (including intron in utr region)
// ----------------------
bool State::IsIntron(void)
{
    if (State2Status[this->state] == SPLICED_TRANSCRIBED)
    {
        return true;
    }
    return false;
}

// ----------------------
//  Return true if it is an intron of a region between a start and stop (excuding intron in utr region)
// ----------------------
bool State::IsIntronInStartStopRegion(void)
{
	if (this->state >= IntronF1 && this->state <= IntronR2AG)
		return true;
	return false;	
}

// ----------------------
//  Return true if it is an intron of an UTR region
// ----------------------
bool State::IsUTRIntron(void)
{
	if ( (this->state >= IntronU5F) && (this->state <= IntronU3R))
		return true;
	return false;
}
// ----------------------
//  Return true if it is an UTR
// ----------------------
bool State::IsUTR(void)
{
    if (State2Status[this->state] == UNTRANSLATED)
    {
        return true;
    }
    return false;
}

// ----------------------
//  Return true if it is an UTR in 5'
// ----------------------
bool State::IsUTR5(void)
{
    if (this->state == UTR5F  ||  this->state == UTR5R)
    {
        return true;
    }
    return false;
}

// ----------------------
//  Return true if it is an UTR in 3'
// ----------------------
bool State::IsUTR3(void)
{
    if (this->state == UTR3F  ||  this->state == UTR3R)
    {
        return true;
    }
    return false;
}

// ----------------------
//  Return true if it is a coding exon
// ----------------------
bool State::IsCodingExon(void)
{
    if (State2Status[this->state] == TRANSLATED)
    {
        return true;
    }
    return false;
}

// ----------------------
//  Return true if coding exon on the strand forward
// NOTE  A modifier! Return (IsCodingExon() && IsForward())
// ----------------------
bool State::IsForwardCodingExon()
{
    if (1 <= State2Frame[this->state] && State2Frame[this->state] <= 3)
    {
        return true;
    }
    return false;
}

// ----------------------
//  Return true if coding exon on the reverse forward
// NOTE  A modifier! Return (IsCodingExon() && IsReverse())
// ----------------------
bool State::IsReverseCodingExon()
{
    if (-3 <= State2Frame[this->state] && State2Frame[this->state] <= -1)
    {
        return true;
    }
    return false;
}

// ----------------------
//  Return true if it is transcribed and unspliced: UTR or coding exon
// ----------------------
bool State::IsTranscribedAndUnspliced(void)
{
    if (this->IsUTR() || this->IsCodingExon())
    {
        return true;
    }
    return false;
}

// -----------------------
// Return true if its an intron in the forward strand
// -----------------------
bool State::IsForwardIntron(void)
{
    if (this->IsIntron() && this->GetStrand() == '+')
    {
        return true;
    }
    return false;
}


// -----------------------
// Return true if its an intron in the reverse strand
// NOTE: remplacer ce test par IsIntron && IsReverse
// -----------------------
bool State::IsReverseIntron(void)
{
    //if (State2Frame[this->state] == FrameIntronR)
    if (this->IsIntron() && this->GetStrand() == '-')
    {
        return true;
    }
    return false;
}


// -----------------------
// Return true if its an element including between a start and a stop codon
// Not IG, UTR and UTR Intron, RNA
// -----------------------
bool State::InStartStopRegion(void)
{
    if (this->state < InterGen)
        return true;
    return false;
}

// ------------------------
//  Return the frame - See the definition on the EuGene trac
// ------------------------
short int State::GetFrame()
{
    return State2Frame[this->state];
}

// ------------------------
//  Return true if the state is defined
// NOTE: for the moment, test if the value of state is positive; this test can be improved
// ------------------------
bool State::IsDefined()
{
    if (this->state >= 0)
        return true;
    return false;
}

// ------------------------
//  Return the strand 
// ------------------------
char State::GetStrand()
{
	char strand = (State2Frame[this->state] > 0) ? '+' : '-';
	    // If state == utr State2Frame returns 0 so
    	if (this->state == UTR5F  ||  this->state == UTR3F) strand = '+';
	return strand;
}

// ------------------------
//  Return true if it is an initial exon
// ------------------------
bool State::IsInitExon()
{
	if (this->state <= InitR3)
		return true;
	return false;
}

// ------------------------
//  Return true if it is an single exon
// ------------------------
bool State::IsSnglExon()
{
	if (this->state >= SnglF1 && this->state <= SnglR3)
		return true;
	return false;
}

// ------------------------
//  Return true if it is a terminal exon
// ------------------------
bool State::IsTermExon()
{
	if (this->state >= TermF1 && this->state <= TermR3)
		return true;
	return false;
}

// ------------------------
//  Return true if it is a non protein coding rna
// ------------------------
bool State::IsNpcRna()
{
	if (this->state == RnaF || this->state == RnaR)
	{
		return true;
	}
	return false;
}



// -----------------------------------------
//  Convert the state in string
// -----------------------------------------
const char* State::State2EGNString()
{
    if (this->IsInitExon())                  return "Init";
    if (this->state <= SnglR3)                  return "Sngl";
    if (this->state <= IntrR3)                  return "Intr";
    if (this->state <= TermR3)                  return "Term";
    if (this->state <= IntronR2AG)              return "Intron";
    if (this->state == InterGen)                return "InterG";
    if (this->state == UTR5F || state == UTR5R) return "Utr5";
    if (this->state == UTR3F || state == UTR3R) return "Utr3";
    if (this->state == RnaR || this->state == RnaF) return "ncRNA ";
    if (this->state >= IntronU5F)               return "Intron";

}

// -----------------------------------------
//  State2GFFString (convert state to char*)
// -----------------------------------------
const char* State :: State2GFFString ()
{
    if (this->IsInitExon())                  return "E.Init";
    if (this->state <= SnglR3)                  return "E.Sngl";
    if (this->state <= IntrR3)                  return "E.Intr";
    if (this->state <= TermR3)                  return "E.Term";
    if (this->state <= IntronR2AG)              return "Intron";
    if (this->state == InterGen)                return "InterG";
    if (this->state == UTR5F || state == UTR5R) return "UTR5";
    if (this->state == UTR3F || state == UTR3R) return "UTR3";
    if (this->state == RnaR || this->state == RnaF) return "ncRNA";
    if (this->state >= IntronU5F)               return "Intron";
}

// -----------------------------------------
//  Return the SOFA type. If sofa==false, return the SO code 
// TODO:  Modifier ce test if (!this->IsDefined()) return SOFA_EXON; 
// -----------------------------------------
int State::getTypeSofa(bool coding, bool sofa)
{
	if (!this->IsDefined()) return SOFA_EXON; 

    if(this->state <= InitR3)
//        str = "five_prime_coding_exon";//"E.Init";
        return (sofa) ? SOFA_EXON : (coding ? SO_5_CODING_EXON : SO_5_EXON);
    if(this->state <= SnglR3)
//        str = "single_exon";//"E.Sngl";
        return (sofa) ? SOFA_EXON : SO_SINGLE_EXON;
    if(this->state <= IntrR3)
        //interior_exon et interior_coding_exon ne sont pas dans SOFA
      return (sofa) ? SOFA_EXON : SO_INTERIOR_EXON;  //"E.Intr";
//        str = "SO:0000201";  //"E.Intr";
    if(this->state <= TermR3)
//        str = "three_prime_coding_exon";//"E.Term";
      return (sofa) ? SOFA_EXON : (coding ? SO_3_CODING_EXON : SO_3_EXON);
    if(this->state <= IntronR2AG)
      return SOFA_INTRON;    //"Intron";
    if(this->state == InterGen)
      return SOFA_INTERGEN;  //"InterG";
    if(this->state == UTR5F || this->state == UTR5R)
      return SOFA_5_UTR;  //"UTR5";
    if(this->state == UTR3F || this->state == UTR3R)
      return SOFA_3_UTR;  //"UTR3";
//    if(state_egn >= IntronU5F)
    if(this->state <= IntronU3R)
      return SO_UTR_INTRON;    //"Intron";
//     si ce n'est rien de tout ca, je renvoies SOFA_EXON
//     a cause des attributs,
//     je dois pouvoir recuperer un SOFA_EXON avec un sofa a false
//     a changer quand la fonction setOntology sera ecrite
    return SOFA_EXON;
}
