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
// $Id: Const.h,v 1.25 2008-10-27 09:26:50 tschiex Exp $
// ------------------------------------------------------------------
// File:     Const.h
// Contents: global variables
// ------------------------------------------------------------------


#ifndef  CONST_H_INCLUDED
#define  CONST_H_INCLUDED

#include <math.h>
// #include <stdlib.h>
#include <climits>
#include <cstdlib>
#include <cstring>

#define EUGENE_DIR getenv("EUGENEDIR")
#define DEFAULT_PARA_FILE "cfg/eugene.par"
#define PLUGINS_DIR "plugins"
#define MODELS_DIR "models"
#define WEB_DIR "web"

#define POSSIBLE_ARGUMENTS "FGREBadrsghm:w:f:n:o:p:x:y:c:u:v:U:V:b::l:O:D:t::M:Z::A:"

const int FASTA_Len = 50;
const int MAX_LINE  = 300;
const int MAX_GFF_LINE  = 2048;
const double NINFINITY   = log(0.0);

// Les Hits EST
const unsigned char HitForward    = 0x1;
const unsigned char MLeftForward  = 0x2;
const unsigned char GapForward    = 0x4;
const unsigned char MRightForward = 0x8; 

// Left shift to go from Hits to ...
const  unsigned int HitToMLeft  = 1;
const  unsigned int HitToGap    = 2;
const  unsigned int HitToMRight = 3;

// Right shift to go from Reverse to Forward
const unsigned int  ReverseToForward = 4; 

const unsigned char HitReverse    = 0x10;
const unsigned char MLeftReverse  = 0x20;
const unsigned char GapReverse    = 0x40;
const unsigned char MRightReverse = 0x80; 

const unsigned char Hit           = HitForward    | HitReverse;
const unsigned char MLeft         = MLeftForward  | MLeftReverse;
const unsigned char MForward      = MLeftForward  | MRightForward;
const unsigned char MReverse      = MLeftReverse  | MRightReverse;
const unsigned char Gap           = GapForward    | GapReverse;
const unsigned char MRight        = MRightForward | MRightReverse;
const unsigned char Margin        = MRight        | MLeft;
const unsigned char AllForward    = HitForward | MLeftForward | GapForward | MRightForward;
const unsigned char AllReverse    = HitReverse | MLeftReverse | GapReverse | MRightReverse;

const unsigned char NotAHit       = Margin | Gap;

#endif
