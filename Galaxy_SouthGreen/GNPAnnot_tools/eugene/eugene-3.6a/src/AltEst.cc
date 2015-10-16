// ------------------------------------------------------------------
// Copyright (C) 2005 INRA <eugene@ossau.toulouse.inra.fr>
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
// $Id: AltEst.cc,v 1.16 2009-09-07 15:32:40 sallet Exp $
// ------------------------------------------------------------------
// File:     AltEst.cc
// Contents: class Alternative Est
// ------------------------------------------------------------------

#include "AltEst.h"
#include <algorithm>
#include "./Hits.h"
extern Parameters   PAR;

/*************************************************************
 **                      OneAltEst
 *************************************************************/
// ----------------------
//  Default constructor.
// ----------------------
OneAltEst :: OneAltEst()
{
    id[0] = 0;
    start = end = index = exonsNumber = 0;
    totalLength = altSplicingEvidence = 0;
}

OneAltEst :: OneAltEst(char *ID, int i, int j)
{
    strcpy(id, ID);
    this->start               = i;
    this->index               = -1;
    this->exonsNumber         = 0;
    this->altSplicingEvidence = 0;
    this->totalLength         = 0;
    this->AddExon(i, j);
}

// ----------------------
//  Default destructor.
// ----------------------
OneAltEst :: ~OneAltEst ()
{
}

// ------------------------
//  ResetEst
// ------------------------
void OneAltEst :: Reset ()
{
    this->start               = -1;
    this->end                 = -1;
    this->index               = -1;
    this->exonsNumber         = 0;
    this->altSplicingEvidence = 0;
    this->totalLength         = 0;
    vi_ExonStart.clear();
    vi_ExonEnd.clear();
}

// ----------------------
//  Push an exon
// ----------------------
void OneAltEst :: AddExon(int i, int j)
{
    vi_ExonStart.push_back(i);
    vi_ExonEnd.push_back(j);
    this->end = j;
    this->exonsNumber++;
    this->totalLength += (j-i+1);
}

// ----------------------
//  Remove an exon
// ----------------------
void OneAltEst :: RemoveExon(int i)
{
    totalLength -= (vi_ExonEnd[i] - vi_ExonStart[i] + 1);
    vi_ExonStart.erase(vi_ExonStart.begin() + i);
    vi_ExonEnd.erase(vi_ExonEnd.begin() + i);
    exonsNumber--;
    if ((i==0) || (i==exonsNumber))
        UpdateBoundaries();
}

// ------------------------------------------
//  ExtremitiesTrim
//   To avoid alignment errors, clean the est
//   by "trimming" the extremities
//   like a pac-man or an exonuclease...
// ------------------------------------------
void OneAltEst :: ExtremitiesTrim(int exonuclen)
{
    if (exonuclen == 0) return;
    // shortening first exon
    //print();
    if (vi_ExonEnd[0] - vi_ExonStart[0] > exonuclen)
    {
        vi_ExonStart[0] += exonuclen;
        totalLength     -= exonuclen;
    }
    else   // tiny first exon, it has to be removed
    {
        //fprintf(stdout,"\n1 %s %d %d\n", id, vi_ExonStart[0], vi_ExonEnd[0]);
        if (exonsNumber>1) // security test
            RemoveExon(0);
    }

    // same shortening for the last exon
    if (vi_ExonEnd[exonsNumber-1]-vi_ExonStart[exonsNumber-1] > exonuclen)
    {
        vi_ExonEnd[exonsNumber-1] -= exonuclen;
        totalLength               -= exonuclen;
    }
    else   // tiny last exon, it has to be removed
    {
        //fprintf(stdout,"2 %s\n", id);
        if (exonsNumber>1) // security test
            RemoveExon(exonsNumber-1);
    }
    UpdateBoundaries();
    //Print();
}

// -------------------------------------------
//  IsFiltered
//   Return 0 if no filtered
//          1 if filtered because of unspliced
//          2 if filtered because of exon len
// -------------------------------------------
int OneAltEst :: IsFiltered(bool unspliced, bool extremelen, bool verbose,
                            int  minIn,     int  maxIn,      int  maxEx, int minEx,
                            int  minEstLen, int  maxEstLen)
{
    // remove if unspliced
    if ((unspliced) && (exonsNumber == 1))
    {
        if (verbose) fprintf(stderr,"\n%s removed (unspliced) ...", id);
        return 1;
    }
    if (extremelen)
    {
        // remove if est too short or too long (after trimming!)
        if (totalLength < minEstLen)
        {
            if (verbose) fprintf(stderr,"\n%s removed (too short) ...", id);
            return 2;
        }
        if (totalLength > maxEstLen)
        {
            if (verbose) fprintf(stderr,"\n%s removed (too long) ...", id);
            return 2;
        }
        if ((exonsNumber>1) && (vi_ExonEnd[0] - vi_ExonStart[0]-1 > maxEx))
        {
            if (verbose) fprintf(stderr,"\n%s removed (exon too long) ...", id);
            return 2;
        }

        // remove if one intron or exon is too short or too long,
        for (int j=1; j<exonsNumber;j++)
        {
            if ((vi_ExonStart[j] - vi_ExonEnd[j-1] -1) < minIn)
            {
                if (verbose) fprintf(stderr,"\n%s removed (intron too short) ...", id);
                fprintf(stderr,"j : %d vi_ExonStart[j] %d vi_ExonEnd[j-1] %d", j, vi_ExonStart[j], vi_ExonEnd[j-1]);
                return 2;
            }
            if ((vi_ExonStart[j] - vi_ExonEnd[j-1] -1) > maxIn)
            {
                if (verbose) fprintf(stderr,"\n%s removed (intron too long) ...", id);
                return 2;
            }
            if ((vi_ExonEnd[j] - vi_ExonStart[j] -1) > maxEx)
            {
                if (verbose) fprintf(stderr,"\n%s removed (exon too long) ...", id);
                return 2;
            }
            if ((vi_ExonEnd[j] - vi_ExonStart[j] -1 < minEx) && (j<exonsNumber-1))
            {
                // internal exon too short
                if (verbose) fprintf(stderr,"\n%s removed (int. exon too short) ...", id);
                return 2;
            }
        }
    }
    return 0;
}

// -----------------------------------------------------
//  Test splice sites identity of 2 Est
//  they have to be sorted: "this" starts before "other"
// -----------------------------------------------------
bool OneAltEst :: IsInconsistentWith(OneAltEst *other)
{
    int i, j = 0, k;
    //printf("%s IsInconsistentWith %s ? ... ",this->GetId(),other->GetId());

    // if this and other are not overlaping they can't be inconsistent
    if (other->start >= this->end) return false;

    // i reaches the first exon that ends after the begin of other,
    // that is the first exon overlaping other
    for (k=0; k<this->exonsNumber; k++)
    {
        if (vi_ExonEnd[k] >= other->start) break;
    }
    // control if other starts in an intron of this
    if (other->start < this->vi_ExonStart[k]) return true;

    for (i=k; i<this->exonsNumber; i++)
    {
        // for each exon i of this (other starts at exon j=0, because this is before other)
        if ((i < (this->exonsNumber-1)) && (j < (other->exonsNumber-1)))
        {
            // they have each at least one remaining exon
            if ( (this->vi_ExonEnd[i]     != other->vi_ExonEnd[j]) ||
                    (this->vi_ExonStart[i+1] != other->vi_ExonStart[j+1]))
                return true; // a difference in the pairs donor - acceptor
            // same coordinates :
            j++;
            continue; // OK, increment the next exon of this (loop "for" upper)
        }
        else   // at least one est is ending
        {
            if (i == this->exonsNumber-1)   // end of this
            {
                if ((other->vi_ExonEnd[j] >= this->end) || (j == other->exonsNumber-1)) return false;
                // current exon of other ends after the last exon of this, so there can't be any inconsistency
                // if not, it has to be the last exon of other
                return true;
                // an intron begin in other while i is still in its last exon
            }
            else   // only other is in its last exon
            {
                if (this->vi_ExonEnd[i] >= other->end) return false;
                // idem than upper, the last exon of other ends in an exon of this,
                // so no inconsistency is possible.
                // If not, the last exon of other continue in an intron of this
                return true;
            }
        }
    }
    // We're not supposed to reach this line!...
    fprintf(stderr, "\n\n INTERNAL ERROR in AlternativeEst::IsInconsistentWith\n");
    return false; // just to avoid a compilation warning
}

// -----------------------------------------------------------------
// Check compatibility of a prediction with an AltEst
// -----------------------------------------------------------------
bool OneAltEst :: CompatibleWith(Prediction *pred)
{
    int idxf, idxe=0;
    Gene *g;

    //locate gene
    g = pred->FindGene(start,end);

    // if no such gene then it is incompatible.
    if (g == NULL) return false;

    // check first exon start is in transcribed matured region
    int  nbFeature = 0;
    nbFeature=g->nbFea();

    for (idxf = 0; idxf < nbFeature; idxf++)
    {
        if ((g->vFea[idxf]->start-1 <= vi_ExonStart[idxe]) &&
                (g->vFea[idxf]->end-1   >= vi_ExonStart[idxe]))
        {
	    if (! g->vFea[idxf]->IsTranscribedAndUnspliced() )
                return false;
            else break;
        }
    }
    //   if (idxf == g->nbFea()) return false;
    idxf++;
    bool firstOk = false;
    // test des fronti�res suivantes
    while ((idxe < exonsNumber-1) && (idxf < nbFeature))
    {
        if (!firstOk && (g->vFea[idxf]->start-1 == vi_ExonEnd[idxe]+1))
            firstOk = true;
	if (!firstOk && ( ! g->vFea[idxf]->IsTranscribedAndUnspliced()) ) // IG or intron: broken
            return false;

        if (firstOk && (g->vFea[idxf]->end == vi_ExonStart[idxe+1]))
            if ( !g->vFea[idxf]->IsTranscribedAndUnspliced() ) // IF or intron
            {
                idxe++;
                idxf++;
                firstOk = false;
                continue;
            }
            else return false;
        idxf++;
    }

    // test de la fin
    for (; idxf < nbFeature; idxf++)
    {
        if ((g->vFea[idxf]->start-1 <= vi_ExonEnd[idxe]) &&
                (g->vFea[idxf]->end-1   >= vi_ExonEnd[idxe]))
	      return  (g->vFea[idxf]->IsTranscribedAndUnspliced());
    }
    return false;
}

// --------------------------------------
//  Penalize according to the EST
// --------------------------------------
void OneAltEst :: Penalize(int pos, DATA *Data, double altPenalty)
{
    if ((pos < start) | (pos > end)) return;

    int idx, inExon = 1;

    for (idx =0; idx < exonsNumber-1; idx++)
    {
        if ((vi_ExonStart[idx] <= pos) &&
                (vi_ExonEnd[idx] >= pos))
        {
            inExon = 1;
            break;
        }

        if ((vi_ExonEnd[idx] < pos) &&
                (vi_ExonStart[idx+1] > pos))
        {
            inExon = 0;
            break;
        }
    }

    //we don't choose the strand. The splice sites will do this for us
    if (inExon == 0)
    {
        Data->contents[DATA::UTR5F] -= altPenalty;
        Data->contents[DATA::UTR3F] -= altPenalty;
        Data->contents[DATA::UTR5R] -= altPenalty;
        Data->contents[DATA::UTR3R] -= altPenalty;
        for (int i=0; i<6; i++)
            Data->contents[i] -= altPenalty;
    }
    else
    {
        Data->contents[DATA::IntronF] -= altPenalty;
        Data->contents[DATA::IntronUTRF] -= altPenalty;
        Data->contents[DATA::IntronR] -= altPenalty;
        Data->contents[DATA::IntronUTRR] -= altPenalty;
    }
    Data->contents[DATA::InterG] -= altPenalty;

}

// --------------------------------------
//  Print the OneAltEst (nuc coordinates)
// --------------------------------------
void OneAltEst :: Print()
{
    fprintf(stdout," %s:", id);
    for (int i=0; i<(int)vi_ExonStart.size(); i++)
    {
        fprintf(stdout,"\t%d - %d", vi_ExonStart[i]+1, vi_ExonEnd[i]+1);
    }
    fprintf(stdout,"\n");
}

bool StartAtLeft( OneAltEst A,  OneAltEst B)
{
    return (A.GetStart() < B.GetStart());
};
bool EndAtRight( OneAltEst A,  OneAltEst B)
{
    return (A.GetEnd() > B.GetEnd());
};
/*************************************************************
 **                        AltEst
 *************************************************************/
// ----------------------
//  Default constructor.
// ----------------------
AltEst :: AltEst(DNASeq *X)
{
    int i,nbIncomp,nbNoevidence,nbIncluded,nbUnspliced, nbExtremLen;
    char tempname[FILENAME_MAX+1];
    i=nbIncomp=nbNoevidence=nbIncluded=nbUnspliced=nbExtremLen=0;

    altPenalty = PAR.getD("AltEst.Penalty");
    includedEstFilter= PAR.getI("AltEst.includedEstFilter");
    compatibleEstFilter= PAR.getI("AltEst.compatibleEstFilter");
    unsplicedEstFilter= PAR.getI("AltEst.unsplicedEstFilter");
    extremeLengthFilter= PAR.getI("AltEst.extremeLengthFilter");
    maxEstLength= PAR.getI("AltEst.maxEstLength");
    minEstLength= PAR.getI("AltEst.minEstLength");
    maxIn= PAR.getI("AltEst.maxIn");
    minIn= PAR.getI("AltEst.minIn");
    maxEx= PAR.getI("AltEst.maxEx");
    minEx= PAR.getI("AltEst.minEx");
    exonucleasicLength= PAR.getI("AltEst.exonucleasicLength");
    altEstDisplay= PAR.getI("AltEst.altEstDisplay");
    verbose= PAR.getI("AltEst.verbose");
    totalAltEstNumber = 0;

    strcpy(tempname, PAR.getC("fstname"));
    strcat(tempname, ".alt.est");
    if (!ProbeFile(NULL,tempname)) return;

    fprintf(stderr, "Reading alt. spl. evidence...");
    fflush(stderr);
    //BEGIN NEW CN
    std::string inputFormat = to_string(PAR.getC("AltEst.format",0,1));
    int NumEST=0;
    Hits * AllEST= NULL;
    if ( inputFormat == "GFF3" )
    {
        strcat(tempname,".gff3");
        GeneFeatureSet * geneFeatureSet = new GeneFeatureSet (tempname);
        AllEST = AllEST->ReadFromGeneFeatureSet(*geneFeatureSet, &NumEST, -1, 0, X);
        delete geneFeatureSet;
    }
    else
    {
        FILE* fEST = FileOpen(NULL, tempname, "r", PAR.getI("EuGene.sloppy"));
        if (fEST)
        {
            //ReadFromFile (EstFile  EstNumber  Level  Margin)
            AllEST = AllEST->ReadFromFile(fEST, &NumEST, -1, 0,X->SeqLen);
            fclose(fEST);
        }
    }
    i=convertHitsToAltEst(AllEST,nbUnspliced,nbExtremLen,NumEST);
    //END NEW CN
    //i=ReadAltFile (tempname,nbUnspliced,nbExtremLen);

    fprintf(stderr, " %d read, ",i);
    fflush(stderr);

    // sort all Est by their begin coordinates
    sort(voae_AltEst.begin(), voae_AltEst.end(), StartAtLeft);

    Compare(nbIncomp, nbNoevidence, nbIncluded);
    fprintf(stderr,"%d removed (%d incl., %d unsp., %d no alt.spl., %d len.), %d inc. pairs, ",
            nbIncluded+nbNoevidence+nbUnspliced+nbExtremLen,
            nbIncluded,nbUnspliced,nbNoevidence,nbExtremLen,nbIncomp);

    fprintf(stderr, "%d kept ...",totalAltEstNumber);
    fprintf(stderr, " done\n");
    fflush(stderr);

    // and set indexes
    for (i=0; i<totalAltEstNumber; i++)
    {
        voae_AltEst[i].PutIndex(i);
        // So, with the special AltEst "Init",
        // indexes are the same than the DAG[] indexes
        if (altEstDisplay) voae_AltEst[i].Print();
    }

    nextAdd = nextRemove = (totalAltEstNumber > 1);

}

// ----------------------
//  Default destructor.
// ----------------------
AltEst :: ~AltEst ()
{
}

// ----------------------
//  Read from Hits class
// ----------------------
int AltEst :: convertHitsToAltEst (Hits * AllEST, int &nbUnspliced, int &nbExtremLen, int &NumEST)
{
    int  read, deb, fin, poids, brin, EstDeb, EstFin, filtertype, totalread;
    OneAltEst oaetmp;
    bool first = 1;
    totalread = 0;
    int i;
    Hits *ThisEST = NULL;
    Block *ThisBlock = NULL;
    for (i = 0, ThisEST = AllEST; i < NumEST; i++, ThisEST = ThisEST->Next)
    {

        ThisBlock = ThisEST->Match;
        while (ThisBlock)
        {
            if (ThisBlock == ThisEST->Match) //si premier
            {
                oaetmp = OneAltEst(ThisEST->getName(), ThisBlock->Start, ThisBlock->End);
            }
            else
            {
                oaetmp.AddExon(ThisBlock->Start, ThisBlock->End);
            }
            ThisBlock = ThisBlock->Next;
        }
        oaetmp.ExtremitiesTrim(exonucleasicLength);
        filtertype = oaetmp.IsFiltered(unsplicedEstFilter, extremeLengthFilter, verbose,
                                       minIn, maxIn, maxEx, minEx,
                                       minEstLength, maxEstLength);
        if (filtertype == 0)
        {
            voae_AltEst.push_back(oaetmp);
            totalAltEstNumber++;
        }
        else
        {
            if (filtertype == 1) nbUnspliced++;
            if (filtertype == 2) nbExtremLen++;
        }
        oaetmp.Reset();
        totalread++;
    }

    return totalread;
}
// ----------------------
//  Read .alt file.
// ----------------------
int AltEst :: ReadAltFile (char name[FILENAME_MAX+1], int &nbUnspliced, int &nbExtremLen)
{
    int  read, deb, fin, poids, brin, EstDeb, EstFin, filtertype, totalread;
    int  PEstFin;
    char *EstId, *PEstId, *tmp;
    char A[128], B[128];
    FILE *fp;
    OneAltEst oaetmp;
    bool first = 1;

    A[0]   = B[0] = 0;
    EstId  = A;
    PEstId = B;
    PEstFin = -1;

    fp = FileOpen(NULL, name, "r");
    totalread = 0;

    while ((read=fscanf(fp,"%d %d %d %*s %d %s %d %d\n",
                        &deb, &fin, &poids, &brin, EstId, &EstDeb, &EstFin)) == 7)
    {

        if ((strcmp(EstId, PEstId) == 0) && (EstDeb > PEstFin))
        {
            //voae_AltEst[totalAltEstNumber-1].AddExon(deb-1,fin-1);
            oaetmp.AddExon(deb-1, fin-1);
        }
        else
        {
            totalread++;
            if (first) first = 0;
            else
            {
                oaetmp.ExtremitiesTrim(exonucleasicLength);
                filtertype = oaetmp.IsFiltered(unsplicedEstFilter, extremeLengthFilter, verbose,
                                               minIn, maxIn, maxEx, minEx,
                                               minEstLength, maxEstLength);
                if (filtertype == 0)
                {
                    voae_AltEst.push_back(oaetmp);
                    totalAltEstNumber++;
                }
                else
                {
                    if (filtertype == 1) nbUnspliced++;
                    if (filtertype == 2) nbExtremLen++;
                }
                oaetmp.Reset();
            }
            oaetmp = OneAltEst(EstId, deb-1, fin-1);
            tmp    = PEstId;
            PEstId = EstId;
            EstId  = tmp;
        }
        PEstFin = EstFin;
    }
    // last est
    oaetmp.ExtremitiesTrim(exonucleasicLength);
    filtertype = oaetmp.IsFiltered(unsplicedEstFilter, extremeLengthFilter, verbose,
                                   minIn, maxIn, maxEx, minEx,
                                   minEstLength, maxEstLength);
    if (filtertype == 0)
    {
        voae_AltEst.push_back(oaetmp);
        totalAltEstNumber++;
    }
    else
    {
        if (filtertype == 1) nbUnspliced++;
        if (filtertype == 2) nbExtremLen++;
    }

    if (read != EOF) fprintf(stderr,"Incorrect ALTEST file !\n");
    fclose(fp);
    return totalread;
}

// -----------------------------------------------------------------
// Filters the AltEst:
//  Remove Est strictly included in another one,
//  TODO: Set a list of Est clusters,
//   (a cluster can be seen as a connex component in a graph
//   built with est as edges and overlap as vertices ;
//   it allows to reduce the time complexity by restricting the
//   pairwise est-est comparisons to est of a same cluster.)
//  And keep only Est confering an evidence of alternative splicing,
//  (that is every EST that is inconsistent with another one).
// -----------------------------------------------------------------
void AltEst :: Compare(int &nbIncomp, int &nbNoevidence, int &nbIncluded)
{
    int i, j, k = 0;
    // Until now, only one cluster is considered (time max)
    // test all pairwise est comparisons in the cluster
    // NxN comparisons could be tested, but it has been reduced to
    // ((NxN)-N)/2 , only one comparison per pair, without the diag. (cf. k)
    // WARNING : voae_AltEst[0] is the special INIT (counted in totalAltEstNumber)

    for (i=0; i<totalAltEstNumber; i++)
    {
        if (compatibleEstFilter || includedEstFilter)
        {
            // Compare this est with the others to check incompatibility or inclusion
            for (j=1+k; j<totalAltEstNumber; j++)
            {
                char *iID = voae_AltEst[i].GetId();
                char *jID = voae_AltEst[j].GetId();
                if (voae_AltEst[i].IsInconsistentWith(&voae_AltEst[j]))
                {
                    if (verbose) fprintf(stderr,"\nincompatibility: %s vs. %s ...", jID, iID);
                    voae_AltEst[i].PutAltSplE(true);
                    voae_AltEst[j].PutAltSplE(true);
                    nbIncomp++;
                }
                else   // no inconsistency
                {
                    // j strictly included in i
                    if ((includedEstFilter) && (voae_AltEst[j].GetEnd() <= voae_AltEst[i].GetEnd()))
                    {
                        if (verbose) fprintf(stderr,"\n%s removed (included in %s) ...", jID, iID);
                        voae_AltEst.erase(voae_AltEst.begin() + j);
                        totalAltEstNumber--;
                        j--;
                        nbIncluded++;
                    }
                }
            }
            k++;
        }
    }

    if (compatibleEstFilter)
    {
        for (i=0; i<totalAltEstNumber; i++)
        {
            if  (! voae_AltEst[i].GetAltSplE())
            {
                if (verbose) fprintf(stderr,"\n%s removed (no alt.spl. evidence) ...", voae_AltEst[i].GetId());
                voae_AltEst.erase(voae_AltEst.begin() + i);
                totalAltEstNumber--;
                i--;
                nbNoevidence++;
            }
        }
    }
}
// -----------------------------------------------------------------
// Penalize according to EST number i
// -----------------------------------------------------------------
void AltEst :: Penalize(int i, int pos, DATA *Data)
{
    voae_AltEst[i].Penalize(pos,Data,altPenalty);
}
