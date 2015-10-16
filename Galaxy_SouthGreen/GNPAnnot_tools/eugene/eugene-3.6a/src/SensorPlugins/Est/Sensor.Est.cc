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
// $Id: Sensor.Est.cc,v 1.62 2010-01-25 17:01:22 sallet Exp $
// ------------------------------------------------------------------
// File:     Sensor.Est.cc
// Contents: Sensor
// ------------------------------------------------------------------

#include "Sensor.Est.h"
#include "../../MSensor.h"

#define EXTREMEMARGIN 1
#define Inconsistent(x) (((x) & Hit) && ((x) & Gap))

extern MasterSensor *MS;
extern Parameters   PAR;

/*************************************************************
 **                        SensorEst                        **
 *************************************************************/
int HitsCompare(const void *A, const void *B)
{
    Hits **UA,**UB;

    UA = (Hits **) A;
    UB = (Hits **) B;

    if ((*UA)->NGaps > (*UB)->NGaps)
        return -1;
    if ((*UA)->NGaps < (*UB)->NGaps)
        return  1;
    if ((*UA)->Length > (*UB)->Length)
        return -1;
    if ((*UA)->Length < (*UB)->Length)
        return  1;
    return strcmp((*UA)->Name,(*UB)->Name); //pour un tri sans ambiguit�
}

int HitsCompareLex(const void *A, const void *B)
{
    Hits **UA,**UB;

    UA = (Hits **) A;
    UB = (Hits **) B;

    if ((*UA)->Start < (*UB)->Start)
        return -1;
    if ((*UA)->Start > (*UB)->Start)
        return 1;
    if ((*UA)->NGaps > (*UB)->NGaps)
        return -1;
    if ((*UA)->NGaps < (*UB)->NGaps)
        return 1;
    if ((*UA)->Length > (*UB)->Length)
        return -1;
    if ((*UA)->Length < (*UB)->Length)
        return 1;
    return 0;
}

int HitsCompareSup(const void *A, const void *B)
{
    Hits **UA,**UB;

    UA = (Hits **) A;
    UB = (Hits **) B;

    if ((*UA)->Support > (*UB)->Support)
        return -1;
    if ((*UA)->Support < (*UB)->Support)
        return 1;
    return 0;
}

// ----------------------
//  Default constructor.
// ----------------------
SensorEst :: SensorEst (int n, DNASeq *X) : Sensor(n)
{
    FILE *fEST;
    char tempname[FILENAME_MAX+1];
    int  i;

    type     = Type_Content;
    HitTable = NULL;
    N        = n;

    vPos.clear();
    vESTMatch.clear();
    fileExt        = PAR.getC("Est.FileExtension", N);
    estM           = PAR.getI("Est.estM",N);
    utrM           = PAR.getI("Est.utrM",N);
    ppNumber       = PAR.getI("Est.PPNumber",N);
    stepid         = PAR.getI("Output.stepid");
    MinDangling    = PAR.getI("Est.MinDangling",N);
    MaxIntron      = PAR.getI("Est.MaxIntron",N);
    MaxIntIntron   = PAR.getI("Est.MaxInternalIntron",N);
    mRNAOnly       = PAR.getI("Est.mRNAOnly", N);
    DonorThreshold = PAR.getD("Est.StrongDonor",N);
    DonorThreshold = log(DonorThreshold/(1-DonorThreshold));

    index = 0;

    ESTMatch = new unsigned char[X->SeqLen+1];
    for (i = 0; i <= X->SeqLen; i++)
        ESTMatch[i] = 0;

    fprintf(stderr,"Reading cDNA hits............");
    fflush(stderr);

    strcpy(tempname, PAR.getC("fstname"));
    strcat(tempname, fileExt);
    NumEST = 0;
    Hits * AllEST=NULL;
      
    inputFormat_ = to_string(PAR.getC("Est.format", GetNumber(),1));

    if ( inputFormat_ == "GFF3" )
    {
      strcat(tempname,".gff3");

      GeneFeatureSet * geneFeatureSet = new GeneFeatureSet (tempname);
      
      AllEST = AllEST->ReadFromGeneFeatureSet(*geneFeatureSet, &NumEST, -1, 0, X);
      delete geneFeatureSet;
      HitTable = ESTAnalyzer(AllEST, ESTMatch, estM, &NumEST, X);
    }
    else
    {
      fEST = FileOpen(NULL, tempname, "r", PAR.getI("EuGene.sloppy"));
      if (fEST)
      {   //ReadFromFile (EstFile  EstNumber  Level  Margin)
	  AllEST = AllEST->ReadFromFile(fEST, &NumEST, -1, 0,X->SeqLen);
	  HitTable = ESTAnalyzer(AllEST, ESTMatch, estM, &NumEST, X);
	  fclose(fEST);
      }
    }
    
    for (i = 0; i<= X->SeqLen; i++)
        if(ESTMatch[i] != 0)
        {
            vPos.push_back      ( i );
            vESTMatch.push_back ( ESTMatch[i] );
        }
    //Print(tempname);
    delete [] ESTMatch;
}

// ----------------------
//  Default destructor.
// ----------------------
SensorEst :: ~SensorEst ()
{
    vPos.clear();
    vESTMatch.clear();

    if(HitTable != NULL)
        delete [] HitTable;
    HitTable = NULL;
}

// --------------------------------
//  Init est.
//  Exploiting spliced alignements
//  against EST and complete cDNA.
// --------------------------------
void SensorEst :: Init (DNASeq *X)
{
    estP = PAR.getD("Est.estP*",N);
    utrP = PAR.getD("Est.utrP*",N);
    spliceBoost    = PAR.getD("Est.SpliceBoost*",N);
    //for(int jj=0;jj<(int)vPos.size();jj++)
    //printf("vPos[%d]:%d\tvESTM[%d]:%d\n",jj,vPos[jj]+1,jj,vESTMatch[jj]);
}

// -----------------------
//  GiveInfo signal est.
// -----------------------
void SensorEst :: GiveInfo (DNASeq *X, int pos, DATA *d)
{
    unsigned char cESTMatch = 0; // current ESTMatch
    unsigned char Oriented = 0;

    // Peut on faire un bete acces sequentiel ?
    if((index != 0                &&  vPos[index-1] >= pos) ||
            (index < (int)vPos.size()  &&  vPos[index]   <  pos))
    {
        // Non... on se repositionne en temps logarithmique
        iter = lower_bound(vPos.begin(), vPos.end(), pos);
        index = iter-vPos.begin();
    }
    // On est juste avant ou sur pos

    // Si on est dessus
    if (index < (int)vPos.size()  &&  vPos[index] == pos)
    {

        cESTMatch = vESTMatch[index];
        Oriented = !((cESTMatch & AllForward)  & ((cESTMatch & AllReverse) >> ReverseToForward));

        // Favor splice sites in marginal exon-intron regions
        if ((cESTMatch & MLeftForward) &&
                d->sig[DATA::Don].weight[Signal::Forward] != 0.0)
            d->sig[DATA::Don].weight[Signal::Forward] += spliceBoost;

        if ((cESTMatch & MRightForward) &&
                d->sig[DATA::Acc].weight[Signal::Forward] != 0.0)
            d->sig[DATA::Acc].weight[Signal::Forward] += spliceBoost;

        if ((cESTMatch & MLeftReverse) &&
                d->sig[DATA::Acc].weight[Signal::Reverse])
            d->sig[DATA::Acc].weight[Signal::Reverse] += spliceBoost;

        if ((cESTMatch & MRightReverse) &&
                d->sig[DATA::Don].weight[Signal::Reverse] != 0.0)
            d->sig[DATA::Don].weight[Signal::Reverse] += spliceBoost;

        // Exon ou UTR Forward (ou Rna Forward if not mRNAOnly)
        // Si on a un Gap EST ou si l'on connait le sens du match EST
        if ((cESTMatch & Gap) ||
                ((cESTMatch & Hit) && !(cESTMatch & HitForward)))
        {
            d->contents[DATA::UTR5F] += estP;
            d->contents[DATA::UTR3F] += estP;
            for(int i=0; i<3; i++)
                d->contents[i] += estP;
            if (!this->mRNAOnly) 
		d->contents[DATA::RNAF] += estP;
        }

        // Exon ou UTR ou Reverse (ou Rna reverse if not mRNAOnly)
        // Si on a un Gap EST ou si l'on connait le sens du match EST
        if ((cESTMatch & Gap) ||
                ((cESTMatch & Hit) && !(cESTMatch & HitReverse)))
        {
            d->contents[DATA::UTR5R] += estP;
            d->contents[DATA::UTR3R] += estP;
            for(int i=3; i<6; i++)
                d->contents[i] += estP;
            if (!this->mRNAOnly)
		d->contents[DATA::RNAR] += estP; 
        }

        // Intron Forward
        // Si on a un Hit EST orienté ou si l'on a un Gap de l'autre côté seulement
        if(((cESTMatch & Hit) && Oriented) ||
                ((cESTMatch & Gap) && !(cESTMatch & GapForward)))
        {
            d->contents[DATA::IntronF] += estP;
            d->contents[DATA::IntronUTRF] += estP;
        }

        // Intron Reverse
        // Si on a un Hit EST orienté ou si l'on a un Gap de l'autre côté seulement
        if(((cESTMatch & Hit) && Oriented) ||
                ((cESTMatch & Gap) && !(cESTMatch & GapReverse)))
        {
            d->contents[DATA::IntronR] += estP;
            d->contents[DATA::IntronUTRR] += estP;
        }

        // Intergenique: tout le temps si on a un match (gap ou hit)
        d->contents[DATA::InterG] += ((cESTMatch & (Gap|Hit)) != 0)*estP;

        // Penalize RNA if there is a match and that mRNAonly is activated
        if (this->mRNAOnly && (cESTMatch & (Gap|Hit)))
        {
          d->contents[DATA::RNAF] += estP;
          d->contents[DATA::RNAR] += estP;
        }   

        d->EstMatch = cESTMatch;  // WARNING : EST -> on est dans intron
        index++;
    }

    // Pour que les UTR soient support�s par un EST
    if (cESTMatch == 0  &&  (int)vPos.size() != 0)
    {
        // Left
        for (int k=1; k<=utrM; k++)
        {
            iter = lower_bound(vPos.begin(), vPos.end(), pos+k);
            if(*iter == pos+k)
            {
                cESTMatch = vESTMatch[iter-vPos.begin()];
                // If only Margin (-> EST extremities) then penalize all utr tracks
                if((cESTMatch & Margin) && !(cESTMatch & Gap) && !(cESTMatch & Hit))
                {
                    d->contents[DATA::UTR5F] += log(utrP);
                    d->contents[DATA::UTR5R] += log(utrP);
                    d->contents[DATA::UTR3F] += log(utrP);
                    d->contents[DATA::UTR3R] += log(utrP);
                    break;
                }
            }
        }
        // Right
        for (int k=utrM; k>0; k--)
        {
            iter = lower_bound(vPos.begin(), vPos.end(), pos-k);
            if(*iter == pos-k)
            {
                cESTMatch = vESTMatch[iter-vPos.begin()];
                // If only Margin (-> EST extremities) then penalize all utr tracks
                if((cESTMatch & Margin) && !(cESTMatch & Gap) && !(cESTMatch & Hit))
                {
                    d->contents[DATA::UTR5F] += log(utrP);
                    d->contents[DATA::UTR5R] += log(utrP);
                    d->contents[DATA::UTR3F] += log(utrP);
                    d->contents[DATA::UTR3R] += log(utrP);
                    break;
                }
            }
        }
    }
}

// -----------------------
//  ESTAnalyzer.
// -----------------------
Hits** SensorEst :: ESTAnalyzer(Hits *AllEST, unsigned char *ESTMatch,
                                int EstM, int *NumEST, DNASeq *X)
{
    int i,j,k;
    int Rejected = 0;
    Hits *ThisEST = NULL;
    Block *ThisBlock = NULL;
    
    
    fprintf(stderr,"%d sequences read\n",*NumEST);
    fflush(stderr);

    // on trie les hits sur le nombre de gaps et la
    // longueur. L'idee est d'eliminer les epissages partiels et de
    // favoriser la longueur de toute facon.
    Hits **HitTable = new Hits *[*NumEST+1];
    for (i = 0, ThisEST = AllEST; i < *NumEST; i++, ThisEST = ThisEST->Next)
    {        HitTable[i] = ThisEST;    }
    // pour memoriser le premier (a liberer)
    HitTable[*NumEST] = AllEST;

    qsort((void *)HitTable, *NumEST, sizeof(void *), HitsCompare);

    for (int index = 0; index < *NumEST; index++)
    {
        int Inc;
        int ExonInc,TheStrand;
        double WorstSpliceF, WorstSpliceR;
        double DonF,AccF,DonR,AccR;
        DATA dTmp;

        ThisEST = HitTable[index];
        // Le veritable  brin est a priori indetermine
        TheStrand = HitForward | HitReverse;
        Inc = 0;
        ExonInc = 0;
        WorstSpliceF = WorstSpliceR = -NINFINITY;

        // Look for each match in the current Hit
        ThisBlock = ThisEST->Match;

        // First hit: delete short dangling
        if ((ThisBlock->Next != NULL) &&
                ((abs(ThisBlock->Start - ThisBlock->End) <= MinDangling) ||
                 (abs(ThisBlock->End - ThisBlock->Next->Start) >= MaxIntron)))
        {
            ThisEST->Match = ThisBlock->Next;
            ThisBlock->Next->Prev = NULL;
            fprintf(stderr,
                    "   [%s]: Suspicious dangling match removed (len. %d, gap %d)\n",
                    ThisEST->Name,abs(ThisBlock->Start - ThisBlock->End),
                    abs(ThisBlock->End - ThisBlock->Next->Start));
            ThisBlock->Next = NULL;
            delete ThisBlock;
            ThisEST->NGaps--;
            ThisBlock = ThisEST->Match;
        }

        // First Step: tries to determine strandedness
        while (ThisBlock)
        {

            // Last block: delete short dangling
            if ((ThisBlock->Next == NULL) && (ThisBlock->Prev != NULL) &&
                    ((abs(ThisBlock->Start-ThisBlock->End) <= MinDangling) ||
                     (abs(ThisBlock->Start-ThisBlock->Prev->End) >= MaxIntron)))
            {
                ThisBlock->Prev->Next = NULL;
                fprintf(stderr,
                        "   [%s]: Suspicious dangling match removed (len. %d, gap %d)\n",
                        ThisEST->Name,abs(ThisBlock->Start-ThisBlock->End),
                        abs(ThisBlock->Start-ThisBlock->Prev->End));
                delete  ThisBlock;
                ThisBlock = NULL;
                ThisEST->NGaps--;
                break;
            }

            // si on a un gap ?
            if (ThisBlock->Prev != NULL)
//             {cerr << "coord hsp" << ThisBlock->Prev->Start << "-" <<ThisBlock->Prev->End <<" , " << ThisBlock->Start << "-" << ThisBlock->End << "\n";
//             cerr << "gen  hsp" << ThisBlock->Prev->LStart << "-" <<ThisBlock->Prev->LEnd <<" , " << ThisBlock->LStart << "-" << ThisBlock->LEnd << "\n";}
            // si on a un gap ?
            if ( (ThisBlock->Prev != NULL) &&
               ( ( (ThisBlock->Prev->LStart <=  ThisBlock->LStart) &&
                   (abs(ThisBlock->LStart-ThisBlock->Prev->LEnd) <= estM)) || 
                 ( (ThisBlock->Prev->LStart > ThisBlock->LStart) &&
                   (abs(ThisBlock->Prev->LStart-ThisBlock->LEnd <= estM) ))
               ))
            {
                DonF = NINFINITY;
                DonR = NINFINITY;
                AccF = NINFINITY;
                AccR = NINFINITY;

                //	printf("[%d %d]\n",ThisBlock->Start,ThisBlock->End);

                for (j = -EstM; j <= EstM; j++)
                {
                    k = Min(X->SeqLen,Max(0,ThisBlock->Prev->End+j+1));
                    MS->GetInfoSpAt(Type_Acc|Type_Don, X, k, &dTmp);

                    DonF = Max(DonF, dTmp.sig[DATA::Don].weight[Signal::Forward]-
                               dTmp.sig[DATA::Don].weight[Signal::ForwardNo]);
                    AccR = Max(AccR, dTmp.sig[DATA::Acc].weight[Signal::Reverse]-
                               dTmp.sig[DATA::Acc].weight[Signal::ReverseNo]);

                    k = Min(X->SeqLen,Max(0,ThisBlock->Start+j));
                    if(MS->GetInfoSpAt(Type_Acc|Type_Don, X, k, &dTmp))
                    {
                        DonR = Max(DonR, dTmp.sig[DATA::Don].weight[Signal::Reverse]-
                                   dTmp.sig[DATA::Don].weight[Signal::ReverseNo]);
                        AccF = Max(AccF, dTmp.sig[DATA::Acc].weight[Signal::Forward]-
                                   dTmp.sig[DATA::Acc].weight[Signal::ForwardNo]);
                    }
                    else
                    {
                        fprintf(stderr,"   WARNING: cDNA hits ignored."
                                " No splices sites predicted !\n");
                        exit(2);
                    }
                }

                //	printf("Extreme splices: %f %f - %f %f\n",DonF,AccF,DonR,AccR);
                WorstSpliceF = Min(WorstSpliceF,DonF);
                WorstSpliceF = Min(WorstSpliceF,AccF);
                WorstSpliceR = Min(WorstSpliceR,DonR);
                WorstSpliceR = Min(WorstSpliceR,AccR);
            }
            ThisBlock = ThisBlock->Next;
        }

        //    printf("Extreme splices: %e %e\n",WorstSpliceF,WorstSpliceR);

        // Tous les blocs ont ete traites
        if (WorstSpliceF == NINFINITY)
            TheStrand &= (~HitForward);
        if (WorstSpliceR == NINFINITY)
            TheStrand &= (~HitReverse);

        // next iteration on the same EST
        ThisBlock = ThisEST->Match;
        while (TheStrand && ThisBlock)
        {
            // Check for consistency with already read Hits
            // The Inc flag will keep the Inconsistency status
            // 1 - inconsistent with a previous EST
            // 2 - no splice site on the borders of a gap
            // 3 - an exon contains STRONG donor on both strands

            for (i = ThisBlock->Start+EstM; !Inc && i <= ThisBlock->End-EstM; i++)
            {

                if (((ESTMatch[i] & Hit) && !(ESTMatch[i] & TheStrand)) ||
                        (Inconsistent(ESTMatch[i] | TheStrand)))
                {
                    fprintf(stderr,"   [%s]: inconsistent hit [%d-%d]\n",
                            ThisEST->Name,ThisBlock->Start+1,ThisBlock->End+1);
                    Inc = 1;
                }
            }

            // si on a un gap
            if ( (ThisBlock->Prev != NULL) &&
               ( ( (ThisBlock->Prev->LStart <=  ThisBlock->LStart) &&
                   (abs(ThisBlock->LStart-ThisBlock->Prev->LEnd) <= estM)) || 
                 ( (ThisBlock->Prev->LStart > ThisBlock->LStart) &&
                   (abs(ThisBlock->Prev->LStart-ThisBlock->LEnd <= estM) ))
               ))
            {
                for (i=ThisBlock->Prev->End+1+EstM; !Inc && i<ThisBlock->Start-EstM; i++)
                    if (((ESTMatch[i] & Gap) && !(ESTMatch[i] & (TheStrand << HitToGap))) ||
                            (Inconsistent(ESTMatch[i] | (TheStrand << HitToGap))))
                    {
                        fprintf(stderr,"   [%s]: inconsistent gap [%d-%d]\n",
                                ThisEST->Name,ThisBlock->Prev->End+2,ThisBlock->Start);
                        Inc = 1;
                    }
            }

            DonF = NINFINITY;
            DonR = NINFINITY;

            // calcul des sites d'epissage internes a l'exon
            for (i = ThisBlock->Start+EstM+1; !Inc && i <= ThisBlock->End-EstM-1; i++)
            {
                MS->GetInfoSpAt(Type_Acc|Type_Don, X, i, &dTmp);
                DonF = Max(DonF, dTmp.sig[DATA::Don].weight[Signal::Forward]-
                           dTmp.sig[DATA::Don].weight[Signal::ForwardNo]);
                DonR = Max(DonR, dTmp.sig[DATA::Don].weight[Signal::Reverse]-
                           dTmp.sig[DATA::Don].weight[Signal::ReverseNo]);
            }
            if (DonF > DonorThreshold)
                ExonInc |= 1;
            if (DonR > DonorThreshold)
                ExonInc |= 2;
            if (ExonInc == 3 && !Inc && !ThisEST->NGaps)
            {
                fprintf(stderr,"   [%s]: Gapless EST with strong donor [%d-%d]\n",
                        ThisEST->Name,ThisBlock->Start+1,ThisBlock->End+1);
                Inc = 3;
            }

            ThisBlock = ThisBlock->Next;
        }

        if (!TheStrand)
        {
            fprintf(stderr, "   [%s]: no matching splice site\n",ThisEST->Name);
            Inc =2;
        }

        // Si une incoherence est detectee, on va jeter la sequence

        if (Inc)
        {
            Rejected++;
            ThisEST->Rejected = 1;

            ThisBlock = ThisEST->Match;
            while (ThisBlock)
            {

                if (TheStrand & HitForward)
                    PlotESTHit(ThisBlock->Start,ThisBlock->End,1,1);
                if (TheStrand & HitReverse)
                    PlotESTHit(ThisBlock->Start,ThisBlock->End,-1,1);

                if (PAR.getI("Output.graph") && (ThisBlock->Prev != NULL) &&
                   ( ( (ThisBlock->Prev->LStart <=  ThisBlock->LStart) &&
                       (abs(ThisBlock->LStart-ThisBlock->Prev->LEnd) <= estM)) || 
                     ( (ThisBlock->Prev->LStart > ThisBlock->LStart) &&
                       (abs(ThisBlock->Prev->LStart-ThisBlock->LEnd <= estM) ))
                   ))
                {
                    if (TheStrand & HitForward)
                        PlotESTGap(ThisBlock->Prev->End,ThisBlock->Start,1,1);
                    if (TheStrand & HitReverse)
                        PlotESTGap(ThisBlock->Prev->End,ThisBlock->Start,-1,1);
                }
                ThisBlock = ThisBlock->Next;
            }
        }
        // sinon on l'exploite
        else
        {
            int LBoundary;
            int RBoundary;
            unsigned char Info;

            ThisBlock = ThisEST->Match;

            while (ThisBlock)
            {

                // Aligners tend to extend beyond the true hit on
                // extremities: we remove EstM on frontiers

                LBoundary = ((ThisBlock->Prev == NULL) ?
                             Min(X->SeqLen, ThisBlock->Start+EstM) :
                             ThisBlock->Start);

                RBoundary = ((ThisBlock->Next == NULL) ?
                             Max(0,ThisBlock->End-EstM) :
                             ThisBlock->End);

                for (i = LBoundary; i <= RBoundary; i++)
                {
                    if ((Info = (ESTMatch[i] & Hit)))
                        ESTMatch[i] |= (Info & TheStrand);
                    else
                        ESTMatch[i] |= TheStrand;
                }

                // Do we want to mark extreme region or alignements as margins ?
#ifdef EXTREMEMARGIN
                if (ThisBlock->Prev == NULL)
                    for (i = ThisBlock->Start; i<LBoundary; i++)
                        ESTMatch[i] |= TheStrand << HitToMLeft;

                if (ThisBlock->Next == NULL)
                    for (i = RBoundary+1; i <= ThisBlock->End; i++)
                        ESTMatch[i] |= TheStrand << HitToMRight;
#endif

                if (PAR.getI("Output.graph"))
                {
                    if (TheStrand & HitForward)
                        PlotESTHit(ThisBlock->Start,ThisBlock->End,1,0);
                    if (TheStrand & HitReverse)
                        PlotESTHit(ThisBlock->Start,ThisBlock->End,-1,0);
                }

                if ( (ThisBlock->Prev != NULL) &&
                     ( ( (ThisBlock->Prev->LStart <=  ThisBlock->LStart) &&
                         (abs(ThisBlock->LStart-ThisBlock->Prev->LEnd) <= estM)) || 
                       ( (ThisBlock->Prev->LStart > ThisBlock->LStart) &&
                         (abs(ThisBlock->Prev->LStart-ThisBlock->LEnd <= estM) ))
                   ))
                {
                    for (i = Max(0,ThisBlock->Prev->End+1-EstM);
                            i < Min(X->SeqLen, ThisBlock->Prev->End+1+EstM); i++)
                        if ((Info = (ESTMatch[i] & MLeft)))
                            ESTMatch[i] |= (Info & (TheStrand << HitToMLeft));
                        else
                            ESTMatch[i] |= TheStrand << HitToMLeft;

                    if (abs(ThisBlock->Prev->End - ThisBlock->Start) < MaxIntIntron)
                    {

                        for (i = ThisBlock->Prev->End+1; i <= ThisBlock->Start-1; i++)
                            if ((Info = (ESTMatch[i] & Gap)))
                                ESTMatch[i] |= (Info & (TheStrand << HitToGap));
                            else
                                ESTMatch[i] |= TheStrand << HitToGap;
                    }
                    else
                        fprintf(stderr, "   [%s]: long intron ignored (%d bases)\n",
                                ThisEST->Name,abs(ThisBlock->Prev->End - ThisBlock->Start));

                    for (i =  Max(0,ThisBlock->Start-1-EstM);
                            i <  Min(X->SeqLen, ThisBlock->Start-1+EstM); i++)
                        if ((Info = (ESTMatch[i] & MRight)))
                            ESTMatch[i] |= (Info & (TheStrand << HitToMRight));
                        else
                            ESTMatch[i] |= TheStrand << HitToMRight;


                    if (PAR.getI("Output.graph"))
                    {
                        if (TheStrand & HitForward)
                            PlotESTGap(ThisBlock->Prev->End,ThisBlock->Start,1,0);
                        if (TheStrand & HitReverse)
                            PlotESTGap(ThisBlock->Prev->End,ThisBlock->Start,-1,0);
                    }
                }
                ThisBlock = ThisBlock->Next;
            }
        }
    }

    if (Rejected)
        fprintf(stderr,"A total of %d/%d sequences rejected\n",Rejected,*NumEST);
    return HitTable;
}

// ----------------------------
//  Plot Sensor information
// ----------------------------
void SensorEst :: Plot(DNASeq *TheSeq)
{}

// ------------------
//  Post analyse
// ------------------
void SensorEst :: PostAnalyse(Prediction *pred, FILE *MINFO)
{
    int trStart = 0, trEnd = 0, cdsStart = 0, cdsEnd = 0;
    int pprocess = PAR.getI("Est.PostProcess",GetNumber());

    if (!pprocess)
        return;
    if (NumEST == 0)
        return;

    fprintf(MINFO, "#=========================================#\n");
    fprintf(MINFO, "#=             Est evidences             =#\n");
    fprintf(MINFO, "#=========================================#\n");

    qsort((void *)HitTable, NumEST, sizeof(void *), HitsCompareLex);

    // Reset static in EST Support or FEA Support
    if (pprocess == 1)
        ESTSupport(NULL,NULL,100,0,100,0,NULL,0);
    else
        FEASupport(NULL,NULL,100,0,100,0,NULL,0,0);

    for(int i=0; i<pred->nbGene; i++)
    {
        trStart  = pred->vGene[i]->trStart;
        trEnd    = pred->vGene[i]->trEnd;
        cdsStart = pred->vGene[i]->cdsStart;
        cdsEnd   = pred->vGene[i]->cdsEnd;

        if (pprocess == 1)
	{
            // Analyse des ESTs par rapport � la pr�diction
            ESTSupport(pred, MINFO, trStart, trEnd, cdsStart, cdsEnd, HitTable, NumEST);
	}
        else
	{
	    //printf("Gene number %d\n",i);
	    //pred->SanityCheck();
            // Analyse de la pr�diction (des features) par rapport aux EST
            FEASupport(pred, MINFO, trStart, trEnd, cdsStart, cdsEnd, HitTable, NumEST, i+1);
	}
    }
}

// -------------------------------------------------------------------------
//  Post analyse (Analyse des ESTs par rapport � la pr�diction)
//    Tdebut/fin = debut/fin de transcript
//    debut/fin  = debut/fin de traduit
// -------------------------------------------------------------------------
void SensorEst :: ESTSupport(Prediction *pred, FILE *MINFO, int Tdebut, int Tfin,
                             int debut, int fin,  Hits **HitTable, int Size)
{
    static int EstIndex;

    int supported = 0;
    int CDSsupported = 0;
    unsigned char *Sup;
    Block *ThisBlock;
    int ConsistentEST,i;
    int from = 0, to = 0, ESTEnd = 0;

    if (pred == NULL)
    {
        EstIndex = 0;
        return;
    }

    Sup = new unsigned char[Tfin-Tdebut+1];

    for (i=0; i <= Tfin-Tdebut; i++)
        Sup[i]=0;

    // si la fin du codant n'a pas ete rencontree
    if (fin == -1)
        fin = Tfin;
    if ((debut == -1) || (debut > Tfin))
        debut = Tfin+1;

    //si l'iteration precedente a atteint l'extremite
    if (EstIndex >= Size)
        EstIndex = Max(0,Size-1);

    // on rembobine....
    while ((EstIndex > 0) && (HitTable[EstIndex]->End > Tdebut))
        EstIndex--;

    if (EstIndex >= 0  &&  HitTable[EstIndex]->End < Tdebut)
        EstIndex++;

    while (EstIndex >=0  &&  EstIndex < Size)
    {
        // le dernier transcrit exploitable est passe
        if (HitTable[EstIndex]->Start > Tfin)
            break;

        ConsistentEST = 1;
        ThisBlock = HitTable[EstIndex]->Match;

        while (ThisBlock && ConsistentEST)
        {
            // si on a un gap
            if (ThisBlock->Prev != NULL)
            {
                // intersection
                from = Max(Tdebut,ThisBlock->Prev->End+1);
                to = Min(Tfin,ThisBlock->Start-1);

                for (i = from; i <= to; i++)
                {
                    //if it's not an intron
		    if (!pred->GetStateAtPos(i+1)->IsIntron())
                        ConsistentEST = 0;
                }
            }

            from = Max(Tdebut,ThisBlock->Start);
            to = Min(Tfin,ThisBlock->End);
            ESTEnd = ThisBlock->End;

            for (i = from; i <= to; i++)
            {
                //if it's not (transcribed and not spliced)
		if (!pred->GetStateAtPos(i+1)->IsTranscribedAndUnspliced())
                    ConsistentEST = 0;
            }
            ThisBlock = ThisBlock->Next;
        }
        fprintf(MINFO, "cDNA  %-12s %7d %7d     %4d     %2d introns    ",
                HitTable[EstIndex]->Name,
                HitTable[EstIndex]->Start+1,HitTable[EstIndex]->End+1,
                HitTable[EstIndex]->Length,HitTable[EstIndex]->NGaps);

        if (HitTable[EstIndex]->Rejected)
            fprintf(MINFO, "Filtered ");
        else if (!ConsistentEST)
            fprintf(MINFO, "Inconsistent");
        else
        {
            if ((HitTable[EstIndex]->Start) <= Tdebut && (ESTEnd >= Tfin))
                fprintf(MINFO, "Full Transcript Support");
            else if ((HitTable[EstIndex]->Start) <= debut && (ESTEnd >= fin))
                fprintf(MINFO, "Full Coding Support");
            else
                fprintf(MINFO, "Support");

            ThisBlock = HitTable[EstIndex]->Match;
            while (ThisBlock)
            {
                if (ThisBlock->Prev != NULL)
                {
                    // intersection
                    from = Max(Tdebut,ThisBlock->Prev->End+1);
                    to = Min(Tfin,ThisBlock->Start-1);

                    for (i = from; i <= to; i++)
                        if (!Sup[i-Tdebut])
                        {
                            Sup[i-Tdebut] = 1;
                            supported++;   // no need to test intron/gap consistency (consistent EST)
                            if ((i >=debut) && (i <=fin))
                                CDSsupported++; // !! introns, should not be counted
                        }
                }
                from = Max(Tdebut,ThisBlock->Start);
                to = Min(Tfin,ThisBlock->End);

                for (i= from; i <= to; i++)
                    if (!Sup[i-Tdebut])
                    {
                        Sup[i-Tdebut] = 1;
                        supported++;    // no need to test exon/hit compatibility (consistent EST)
                        if ((i >=debut) && (i <=fin))
                            CDSsupported++;
                    }
                ThisBlock = ThisBlock->Next;
            }
        }
        fprintf(MINFO, "\n");
        EstIndex++;
    }
    if (fin >= debut)
        fprintf(MINFO, "      CDS          %7d %7d    %5d     supported on %d bases\n",
                debut+1,fin+1,fin-debut+1,CDSsupported);
    fprintf(MINFO, "      Gene         %7d %7d    %5d     supported on %d bases\n",
            Tdebut+1,Tfin+1,Tfin-Tdebut+1,supported);
    delete [] Sup;
    return;
}

// -------------------------------------------------------------------------
//  Post analyse (Analyse de la pr�diction (features) par rapport aux EST)
//    Tdebut/fin = debut/fin de transcript
//    debut/fin  = debut/fin de traduit
// -------------------------------------------------------------------------
void SensorEst :: FEASupport(Prediction *pred, FILE *MINFO,int Tdebut,int Tfin,
                             int debut,int fin,Hits **HitTable,int Size,int NumGene)
{
    static int EstIndex;

    Block *ThisBlock;
    int ConsistentEST, i, j;
    int from = 0, to = 0;
    int start, end, len;
    State *featState;
    std::vector <int> vSupEstI;  // index des transcrits supportant la pred

    if (pred == NULL)
    {
        EstIndex = 0;
        return;
    }

    /*************************************************************************/
    /** First, create the coordinates of the gene associated transcript     **/
    /** vTranscriptStarts will contain the start of each region of the mRNA **/
    /*************************************************************************/
   int currentStateTranscribed = -1;
   std::vector<int> vTranscriptStarts;

   //printf("Gene ");

   // Check first feature: transcribed or not (partial gene)
   i = 0;
   if (!pred->vGene[NumGene-1]->vFea[i]->IsTranscribedAndUnspliced()) i++;
	
   for (; i < pred->vGene[NumGene-1]->nbFea(); i++)
   {
	start = pred->vGene[NumGene-1]->vFea[i]->start; 
	assert(pred->vGene[NumGene-1]->vFea[i]->IsIntergenic() == false); // no intergenic region here
	
	// si on change d etat : si on passe de IG ou INTRON vers transcrit episse (UTR ou EXON), ou le contraire
	if (pred->vGene[NumGene-1]->vFea[i]->IsTranscribedAndUnspliced() != currentStateTranscribed)
	{
	    vTranscriptStarts.push_back(start-1);
	    currentStateTranscribed =  pred->vGene[NumGene-1]->vFea[i]->IsTranscribedAndUnspliced();
	}
   }

    if (pred->vGene[NumGene-1]->vFea[i-1]->IsTranscribedAndUnspliced() )
    {
	vTranscriptStarts.push_back(pred->vGene[NumGene-1]->vFea[i-1]->end+1); // last exon
    }


    /***********************************************************************/
    /** Objectif : obtenir un vecteur contenant les index des transcrits  **/
    /**            qui supportent la prédiction                           **/
    /***********************************************************************/
    // si la fin du codant n'a pas ete rencontree
    if (fin == -1)
        fin = Tfin;
    if ((debut == -1) || (debut > Tfin))
        debut = Tfin+1;

    // si l'iteration precedente a atteint l'extremite
    if (EstIndex >= Size)
        EstIndex = Max(0, Size-1);

    // rewinding
    while ((EstIndex > 0) && (HitTable[EstIndex]->End > Tdebut))
        EstIndex--;

    if (EstIndex >= 0  &&  HitTable[EstIndex]->End < Tdebut)
        EstIndex++;

    while (EstIndex >= 0  &&  EstIndex < Size)
    {
        // le dernier transcrit exploitable est passé
        if (HitTable[EstIndex]->Start > Tfin)
            break;

        ConsistentEST = 1;
        ThisBlock = HitTable[EstIndex]->Match;
        from   = Max(Tdebut, ThisBlock->Start-1);
	//printf("Est %s starts at %d\n",HitTable[EstIndex]->Name,from);
	int exonIndex = 0;

	// find the exon that overlaps the first match start.
	for (exonIndex=0; ((exonIndex+1 < vTranscriptStarts.size()) && (from > vTranscriptStarts[exonIndex+1]-1)); exonIndex += 2)
	{
	//printf("at exonindex %d %d %d\n",exonIndex,vTranscriptStarts.size(),vTranscriptStarts[exonIndex]);
	}

	if (exonIndex >= vTranscriptStarts.size())
	{
	    // no matching exon found
            //printf("No exon matching first match starting at %d found\n",start);
	    ConsistentEST = 0;
	    break;
	}

	// We can assume that the first match starts in the current ExonIndex.
        while (ThisBlock && ConsistentEST && exonIndex+1 < vTranscriptStarts.size())
        {
            // The GAP case
            if (ThisBlock->Prev != NULL)
            {
                from = Max(Tdebut, ThisBlock->Prev->End+1);
                to   = Min(Tfin,   ThisBlock->Start-1);
	    
                //printf("Est GAP (%d-%d) compared to (%d,%d)\n",from,to,vTranscriptStarts[exonIndex],vTranscriptStarts[exonIndex+1]-1);
	        ConsistentEST &= ((from >= vTranscriptStarts[exonIndex]) && to < vTranscriptStarts[exonIndex+1]);
	        //printf("G exonindex %d size %d\n",exonIndex,vTranscriptStarts.size());

	        exonIndex++;
            }

 	    // the match case
            from   = Max(Tdebut, ThisBlock->Start);
            to     = Min(Tfin,   ThisBlock->End);
	    
	    //printf("Est MATCH (%d-%d) compared to in (%d,%d)\n",from,to,vTranscriptStarts[exonIndex],vTranscriptStarts[exonIndex+1]-1);
	    //printf("M exonindex %d size %d\n",exonIndex,vTranscriptStarts.size());
	    ConsistentEST &= ((from >= vTranscriptStarts[exonIndex]) && (to < vTranscriptStarts[exonIndex+1]));

            ThisBlock = ThisBlock->Next;
	    exonIndex++;
        }

        // Si EST "coherente"
        if (ConsistentEST)
            vSupEstI.push_back( EstIndex );
        EstIndex++;
    }
    vTranscriptStarts.clear();

    // Sup will contain for each base of the transcript which are supported or not.
    unsigned char *Sup = new unsigned char[Tfin-Tdebut+1];
    for (i=0; i<=Tfin-Tdebut; i++) Sup[i] = 0;

    // CDSSupport and GeneSupport will store the number of base supported by each protein, overall 
    int* CDSSupport = new int[vSupEstI.size()];
    int* GeneSupport = new int[vSupEstI.size()];
    for (j=0; j<(int)vSupEstI.size(); j++) 
    {
	CDSSupport[j] = GeneSupport[j] = 0;
    }

    /***********************************************************************/
    /* now, vSupEstI is a set of EST consistent with the analyzed gene     */
    /* Objectif : Analyser chaque feature predite -> supportée ?           */
    /***********************************************************************/
    char fea[5];
    bool inCDS;
    int nlen, transcLen=0, translLen=0;
    char strand = pred->vGene[NumGene-1]->vFea[0]->strand;
    Hits **TMPHitTable = new Hits *[vSupEstI.size()]; // need for sorting by % support

    for (i = 0; i < pred->vGene[NumGene-1]->nbFea(); i++)
    {
	featState = pred->vGene[NumGene-1]->vFea[i]->featureState; 
        start     = pred->vGene[NumGene-1]->vFea[i]->start;
        end       = pred->vGene[NumGene-1]->vFea[i]->end;

	// Coord of tDebut/Tfin start at 0. start/end start at 1
	//printf("NumGene %i state %i start %d Tdebut+1 %d end %d Tfin+1 %d\n", NumGene, state,start,Tdebut+1,end,Tfin+1);

	assert((start >= Tdebut+1) && (end <= Tfin+1));
        len      = 0;
        int numF = -1;

	if (pred->vGene[NumGene-1]->vFea[i]->IsTranscribedAndUnspliced())
            numF = pred->vGene[NumGene-1]->vFea[i]->number;

	inCDS = ((start-1 >= debut) && (end-1 <= fin));

	if (featState->IsCodingExon())
        {
            strcpy(fea, "Exon");
	}
	else if (featState->IsUTR5())
            strcpy(fea, "UTR5");
	else if (featState->IsUTR3())
            strcpy(fea, "UTR3");

        len = 0;
	for (j=0; j<(int)vSupEstI.size(); j++)
        {
	    nlen = 0;
	    //printf("offset in Sup (start-debut-1) %d start %i debut %i\n",start-Tdebut-1,start,debut);
	    HitTable[vSupEstI[j]]->Support = LenSup(HitTable, Sup+start-Tdebut-1, vSupEstI, nlen , j, start, end);
            len += nlen;
	    GeneSupport[j] += HitTable[vSupEstI[j]]->Support;
            transcLen+= nlen; //transcript overall count
	    if (inCDS)
	    {
	        translLen += nlen;
		CDSSupport[j] += HitTable[vSupEstI[j]]->Support;
	    }
	}

        if ((numF != -1) && (len > 0))
        {
            fprintf(MINFO, "%s.%d.%d\tEuGene_cDNA\t%s\t%d\t%d\t%d\t%c\t.\t",
                    pred->seqName, (((NumGene-1)*stepid)+1), numF,
                    fea, start, end, len, strand);
            fprintf(MINFO, "%d\t", (int)((float)len/(end-start+1)*100));
             // On copie la hittable pour trier sur le % support
            for (int k=0; k<vSupEstI.size(); k++)
                TMPHitTable[k] = HitTable[vSupEstI[k]];
            qsort((void*)TMPHitTable, vSupEstI.size(), sizeof(void*), HitsCompareSup);
            // On affiche les ppNumber premiers hits supportant
            for (j=0; j<vSupEstI.size() && j<ppNumber && TMPHitTable[j]->Support!=0; j++)
                fprintf(MINFO, "%s(%d,%d) ", TMPHitTable[j]->Name, (int)((float)TMPHitTable[j]->Support/(end-start+1)*100), !TMPHitTable[j]->Rejected);
            fprintf(MINFO, "\n");
        }
    }
 

    /***********************************************************************/
    /* Now overall Gene/Transcript statistics                              */
    /***********************************************************************/
    int * support;

    for (int i=0; i<2; i++)
    {
        if (i==0)
        {
            start = debut;
            end = fin;
            strcpy(fea, "CDS");
	    len = translLen;
	    support = CDSSupport;
        }
        else
        {
            start = Tdebut;
            end = Tfin;
            strcpy(fea, "Gene");
	    len = transcLen;
	    support = GeneSupport;
        }
        assert(end >= start);

        if (len > 0)
        {
            fprintf(MINFO, "%s.%d  \tEuGene_cDNA\t%s\t%d\t%d\t%d\t%c\t.\t",
                    pred->seqName, (((NumGene-1)*stepid)+1), fea,
                    start+1, end+1, len, strand);
            fprintf(MINFO, "%d\t", (int)((float)len/(end-start+1)*100));

            for (j=0; j<(int)vSupEstI.size(); j++)
            {
		//printf("EST %s Support %i\n",HitTable[vSupEstI[j]]->Name,support[j]);
                HitTable[vSupEstI[j]]->Support = (int)((float)support[j]/(end-start+1)*100);
            }
            // On copie la hittable pour trier sur le % support�
            for (int k=0; k<vSupEstI.size(); k++)
                TMPHitTable[k] = HitTable[vSupEstI[k]];
            qsort((void*)TMPHitTable, vSupEstI.size(), sizeof(void*), HitsCompareSup);

            // On affiche les ppNumber premiers hits supportant
            for (j=0; j<vSupEstI.size() && j<ppNumber && TMPHitTable[j]->Support!=0; j++)
                fprintf(MINFO, "%s(%d,%d) ", TMPHitTable[j]->Name, TMPHitTable[j]->Support, !TMPHitTable[j]->Rejected);
            fprintf(MINFO, "\n");
        }
    }

    vSupEstI.clear();

    if (TMPHitTable != NULL)
        delete [] TMPHitTable;
    TMPHitTable = NULL;
    delete [] CDSSupport;
    CDSSupport = NULL;
    delete [] GeneSupport;
    GeneSupport = NULL;
    return;
}

// -------------------------------------------------------------------------
//  Length supported by EST. The EST are assumed to be 100% consistent
// returns the number of nuc matched by the EST on the feature region
// and increments the overall (all EST) number of nucleotides matching
// -------------------------------------------------------------------------
int SensorEst :: LenSup(Hits **HitTable, unsigned char *Sup, 
			std::vector<int> vSupEstI, int &additionalsup,
                        int index, int beg, int end)
{
    int supported = 0;
    Block *ThisBlock;
    int from = 0, to = 0;
    int i;
    int j;

    ThisBlock = HitTable[vSupEstI[index]]->Match;

    while (ThisBlock)
    {

        if (ThisBlock->Prev != NULL)
        {
            // intersection
            from = Max(beg-1, ThisBlock->Prev->End+1);
            to   = Min(end-1, ThisBlock->Start-1);

	    //if (from <= to) printf("GAP %d-%d => %d\n",from,to,to-from+1); 
            for (j=from; j<=to; j++) 
	    {
                if (!Sup[j-from])
                {
                    Sup[j-from] = 1;
                    additionalsup++;
                }
	        supported++;
	    }
        }
        //printf("G Supported %d Additional %i\n",supported,additionalsup);

        from = Max(beg-1, ThisBlock->Start);
        to   = Min(end-1, ThisBlock->End);

	//if (from <= to) printf("MATCH %d-%d => %d\n",from,to,to-from+1); 
        for (j=from; j<=to; j++)
	{
	    //printf("j-from %d\n",j-from);
            if (!Sup[j-from])
            {
                Sup[j-from] = 1;
                additionalsup++;
            } 
	    supported++;
	}
	//printf("M Supported %d Additional %i\n",supported,additionalsup);
        ThisBlock = ThisBlock->Next;
    }

    return supported;
}


void SensorEst :: Print (char name[FILENAME_MAX+1])
{
  FILE *fp;
  strcat (name, ".out");
  if (!(fp = fopen(name, "w"))) {
    fprintf(stderr, "cannot write in %s\n",  name);
    exit(2);
  }

  fprintf(fp, "vPos %d\n",  vPos.size());
  for (int i=0; i< vPos.size();i++ )
  {
    fprintf(fp, "vPos %d\t%d\n",i,  vPos[i]);
  }
  
  fprintf(fp, "vESTMatch %d\n",  vESTMatch.size());
  for (int i=0; i< vESTMatch.size();i++ )
  {
    fprintf(fp, "vESTMatch %d\t\n", vESTMatch[i]);
  }

  fclose(fp);
}
