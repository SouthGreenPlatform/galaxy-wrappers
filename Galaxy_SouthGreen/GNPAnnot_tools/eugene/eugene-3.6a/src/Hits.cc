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
// $Id: Hits.cc,v 1.22 2009-05-12 15:42:23 tschiex Exp $
// ------------------------------------------------------------------
// File:     Hits.cc
// Contents: Definitions for a class representing alignements with genomic seq.
// ------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#include <stdlib.h>
#include <assert.h>
#include <algorithm>
#ifdef STDC_HEADERS
#include <string.h>
#else
#include <strings.h>
#endif

#include "Const.h"
#include "Hits.h"
#include "System.h"

// ---------------------------------------------------------------------
//  Default constructor.
// ---------------------------------------------------------------------
Block :: Block ()
{
    Prev = Next = NULL;
    HitSeq = NULL;
}

// ---------------------------------------------------------------------
//  Constuct from Deb/Fin
// ---------------------------------------------------------------------
Block ::  Block(int start, int end,int lstart,int lend,int Ph,int Scr)
{
    //fprintf (stderr, "New block : %d %d %d %d %d %d %c \n",start,end,lstart,lend,Ph,Scr);
    Start  = start;
    End    = end;
    LStart = lstart;
    LEnd   = lend;
    Phase  = Ph;
    Score  = Scr;
    Prev   = Next = NULL;
    HitSeq = NULL;
}

// ---------------------------------------------------------------------
//  Create a new Block and insert it after
// ---------------------------------------------------------------------
void Block ::  AddBlockAfter(int start,int end,int lstart,int lend,int Ph,int Scr, char *HSP)
{

    Block *ABlock = new Block(start,end,lstart,lend,Ph,Scr);
    this->Next    = ABlock;
    ABlock->Prev  = this;
    HitSeq = HSP;
}

// ---------------------------------------------------------------------
//  Insert newBlock after
// ---------------------------------------------------------------------
void Block::AddBlockAfter(Block* newBlock, char *HSP)
{
    this->Next     = newBlock;
    newBlock->Prev = this;
    HitSeq = HSP;
}

// ---------------------------------------------------------------------
//  Compare the two blocks according to their genomic positions
// ---------------------------------------------------------------------
bool BlockPosCompare(const Block* a, const Block* b)
{
    return a->Start < b->Start;
}

// ---------------------------------------------------------------------
//  Destroy this Block and all blocks after him
// ---------------------------------------------------------------------
Block :: ~ Block  ()
{
    if (HitSeq) free(HitSeq);
    delete Next;
}

// ---------------------------------------------------------------------
//  Default constructor.
// ---------------------------------------------------------------------
Hits :: Hits  ()
{
    Rejected = 0;
    Start    = 0;
    End      = 0;
    Evalue   = 0.0;
    Level    = 0;
    Support  = 0;
    Length   = 0;
    Strand   = 0;
    NGaps    = 0;
    Name     = NULL;
    Match    = NULL;
    Next     = NULL;
}

// ---------------------------------------------------------------------
//  Construct from a char*...
// ---------------------------------------------------------------------
Hits :: Hits  (char* name, int length, char strand, int deb, int fin,
               int ldeb, int lfin, int Ph, int Scr, double Prob, int level,
               int sup)
{
    Rejected = 0;
    Name     = new char[strlen(name)+1];
    strcpy(Name,name);
    Strand = strand;
    Length = length;
    Start  = deb;
    End    = fin;
    Evalue = Prob;
    Level  = level;
    Support= sup;
    Match  = new Block(deb,fin,ldeb,lfin,Ph,Scr);
    NGaps  = 0;
    Next   = NULL;
}

// ---------------------------------------------------------------------
//  Read a table of Hits from a file
// ---------------------------------------------------------------------
Hits* Hits::ReadFromFile(FILE* HitFile, int *NumHits, int level, int margin, int maxPos)
{
    char   *HitId, *PHitId;
    int    deb, fin, phase, Pphase, HSPDeb, HSPFin, poids, read;
    double evalue, Pevalue;
    char   A[512], B[512];
    char *HSP = NULL;
    Block *ThisBlock = NULL;
    Hits  *OneHit    = NULL, *ThisHit = this, *AllHit = this;
    const int MaxHitLen = 15000;

    Pevalue = -1.0;
    Pphase = 0;
    A[0]    = B[0] = 0;
    HitId   = A;
    PHitId  = B;

    if (ThisHit != NULL)
        for (int i=0; i<*NumHits-1; i++) ThisHit = ThisHit->Next;

    while ((read=fscanf(HitFile,"%d %d %d %lf %d %s %d %d %as\n", &deb, &fin,
                        &poids, &evalue, &phase, HitId, &HSPDeb, &HSPFin,HSP)) >= 8)
    {
        if (HSP) fprintf(stderr,HSP);
        if (phase < 0  &&  deb > fin)
        {
            int tmp = deb;
            deb     = fin;
            fin     = tmp;
            tmp     = HSPDeb;
            HSPDeb  = HSPFin;
            HSPFin  = tmp;
        }

        if ((deb>maxPos) || (fin>maxPos))
        {
            fprintf(stderr,"Hit does not map on sequence. Check %s\n",
                    HitId);
            exit(1);
        }

        if (abs(fin-deb) > MaxHitLen)
        {
            fprintf(stderr,"Similarity of extreme length rejected. Check %s\n",
                    HitId);
            continue;
        }
        if ((strcmp(HitId,PHitId) == 0)      && (phase*Pphase >= 0) &&
                (deb + margin > ThisBlock->End)  &&
                (phase >= 0  &&  (HSPDeb + margin > ThisBlock->LEnd)    ||
                 phase < 0   &&  (HSPDeb - margin < ThisBlock->LEnd)))
            // si HitId et PHitId sont egaux, alors il y a un Hit en cours
            // de meme nom on verifie que c'est bien compatible en terme
            // de position (sur l'est et le genomique) et en e-value, en
            // prenant en compte le brin !
        {
            ThisHit->NGaps++;
            ThisHit->End = fin-1;
            ThisBlock->AddBlockAfter(deb-1,fin-1,HSPDeb,HSPFin,phase,poids);
            ThisBlock = ThisBlock->Next;

        }
        else
        {
            (*NumHits)++;
            OneHit = new Hits(HitId, poids, (phase>0 ? '+' : '-'), deb-1, fin-1,
                              HSPDeb, HSPFin, phase, poids, evalue, level, 0);
            ThisBlock = OneHit->Match;
            //CN Correct bug, assignment wasn't on data but on adress : PHitId = HitId
            strcpy (PHitId , HitId );

            Pevalue   = evalue;
            Pphase = phase;

            if (AllHit == NULL)
            {
                AllHit = OneHit;
                ThisHit = OneHit;
            }
            else
            {
                ThisHit->Next = OneHit;
                ThisHit       = OneHit;
            }
        }
    }
    if (read != EOF)
        fprintf(stderr,"\nIncorrect similarity file after seq. %s\n",HitId);


    return AllHit;
}

// ---------------------------------------------------------------------
//  Read a table of Hits from a file
// ---------------------------------------------------------------------
Hits* Hits::ReadFromGeneFeatureSet(GeneFeatureSet & HitSet , int *NumHits, int level, int margin, DNASeq *X )
{
    int maxPos = X->SeqLen;
    int    deb, fin, phase, HSPDeb, HSPFin, poids;
    double hitEvalue;
    string idSo;
    char   HitId[512], strand, hitStrand;
    Block *ThisBlock = NULL;
    Hits  *OneHit    = NULL, *ThisHit = this, *AllHit = this;
    Block * newBlock;
    // Blocks of one parent
    vector<Block*> vBlocks;
    const int MaxHitLen = 15000;
    HitId[0]   = 0;

    //  because this Hit mais be non empty at the method's start
    if (ThisHit != NULL)
        for (int i=0; i<*NumHits-1; i++) ThisHit = ThisHit->Next;

    // Iterator to scan parent after parent
    map<string, vector<GeneFeature *> >::iterator itParent = HitSet.getIteratorParentToChildren();
    int nbParents = HitSet.getNbParentFeature();

    for (int j=0 ; j < nbParents; j++, itParent++)
    {
        vBlocks.clear();
        if (itParent->first == "") continue;

        vector<GeneFeature *>::iterator it = itParent->second.begin();
        int nbGeneFeature = itParent->second.size();
        int i=0;
        // Create one block per child
        for ( i = 0 ; i < nbGeneFeature ; i++, it++ )
        {
            if ( ! (*it)->hasTarget() ) continue;

            deb    = (*it)->getLocus()->getStart();
            fin    = (*it)->getLocus()->getEnd();
            poids  = (*it)->getAttributes()->getTarget()->getScoreHit();
            strand = (*it)->getLocus()->getStrand();
            idSo   = (*it)->getType();
            if ( poids <= 0)
            {
                poids=(*it)->getLength();
            }
            //Get SO code if feature correspond to the name or the synonym.
            if ( idSo.find("SO:") == string::npos )
            {
                string tmp=GeneFeatureSet::soTerms_->getIdFromName(idSo);
                idSo=tmp;
            }
            strcpy (HitId, (*it)->getAttributes()->getTarget()->getName().c_str());
            HSPDeb = (*it)->getAttributes()->getTarget()->getLocus()->getStart();
            HSPFin = (*it)->getAttributes()->getTarget()->getLocus()->getEnd();

            if ( idSo == "SO:0000668" ) //EST_match
            {
                phase  = (strand == '+' ? 0 : 1);
            }
            else //BlastX
            {
                phase  = (*it)->getAttributes()->getTarget()->getFrameHit();
                if (phase == 0) //not specified
                {
                    phase  = X->Pos2Frame(deb,strand);
                }
                else
                {
                    int computedFrame=X->Pos2Frame(deb,strand);
                    if (phase != computedFrame) //check frame computrd and read are the same.
                    {
                        fprintf( stderr, "Computed frame (%d) and input frame (%d) are different. Check : %s %d %d %c %d %d\n",computedFrame, phase, HitId, deb, fin ,strand, HSPDeb, HSPFin);
                        fflush(stderr);
                        continue;
                    }
                }
            }

            if (phase < 0  &&  strand == '-' )
            {
                int tmp = HSPDeb;
                HSPDeb  = HSPFin;
                HSPFin  = tmp;
            }
            if ((deb>maxPos) || (fin>maxPos))
            {
                fprintf(stderr,"Hit does not map on sequence. Check %s\n", HitId);
                exit(1);
            }
            if (abs(fin-deb) > MaxHitLen)
            {
                fprintf(stderr,"Similarity of extreme length rejected. Check %s\n",HitId);
                continue;
            }

            // first child of the hits: save the strand and the evalue
            if (i == 0)
            {
                hitStrand = strand;
                hitEvalue = (*it)->getScore();
            }
            // Add a new block and save it in the vector
            newBlock = new Block(deb-1,fin-1,HSPDeb,HSPFin,phase,poids);
            vBlocks.push_back(newBlock);
        }

        if (vBlocks.size() > 0)
        {
            // sort the block according to the genomic positions
            std::sort(vBlocks.begin(), vBlocks.end(), BlockPosCompare);

            // Get the first block
            ThisBlock = vBlocks[0];
            // create a new hit and the first block of the hit
            OneHit = new Hits(HitId, ThisBlock->Score, hitStrand, ThisBlock->Start, ThisBlock->End, ThisBlock->LStart, ThisBlock->LEnd, ThisBlock->Phase, ThisBlock->Score, hitEvalue, level, 0);
            (*NumHits)++;
            ThisBlock = OneHit->Match;

            // Add other blocks to the hits and inc the number of gaps
            for (int k = 1; k< vBlocks.size(); k++)
            {
                ThisBlock->AddBlockAfter(vBlocks[k]);
                ThisBlock = ThisBlock->Next;
                OneHit->NGaps++;
            }
            // Update the end position of the Hit = end position of the last block
            OneHit->End   = vBlocks[vBlocks.size()-1]->End;

            if (AllHit == NULL)
            {
                AllHit   = OneHit;
                ThisHit  = OneHit;
            }
            else
            {
                ThisHit->Next = OneHit;
                ThisHit       = OneHit;
            }
        }
    }


    return AllHit;
}


// ---------------------------------------------------------------------
//  Destroy this alignement
// ---------------------------------------------------------------------
Hits :: ~ Hits  ()
{
    delete Match;
    delete [] Name;
    delete Next;
}
char * Hits::getName ()
{
    return Name;
}
