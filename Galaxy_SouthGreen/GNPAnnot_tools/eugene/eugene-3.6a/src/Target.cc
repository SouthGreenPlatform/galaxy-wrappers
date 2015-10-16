/***************************************************************************
 *   Copyright (C) 2007 by Noirot   *
 *   cnoirot@milieu   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "Target.h"

// -----------------------
//  Default constructeur
// -----------------------
Target::Target ( ) 
{
  name_ = "";
  locus_ = NULL;
  sequenceData_ = "";
  targetLength_ = 0;
  isFullLength_ = UNKNOWN; 
  frameHit_=0;
  scoreHit_=-1;
}

Target::Target (std::string name , std::string sequenceData,  int targetLength, int isFullLength, Locus * locus, int frameHit, int scoreHit) 
{
  name_ = name;
  locus_ = new Locus (*locus);
  sequenceData_ = sequenceData;
  targetLength_ = targetLength;
  isFullLength_ = isFullLength; 
  frameHit_     = frameHit;
  scoreHit_     = scoreHit;
}

Target::Target ( std::string name, std::string sequenceData,  int targetLength, int isFullLength, int start, int end, char strand, int frameHit, int scoreHit)
{
        name_ = name;
	locus_ = new Locus (start,end,strand);
	sequenceData_ = sequenceData;
	targetLength_ = targetLength;
	isFullLength_ = isFullLength;
	frameHit_     = frameHit;
	scoreHit_     = scoreHit;
	//cout << "Constructeur TARGET :  "<<getString()<<endl;
}
		
Target::Target (Target * target) 
{
  name_ = target->getName();
  locus_ = new Locus (*target->locus_);
  sequenceData_ = target->sequenceData_;
  targetLength_ = target->targetLength_;
  isFullLength_ = target->isFullLength_;
  frameHit_     = target->frameHit_;
  scoreHit_     = target->scoreHit_;
}
// -----------------------
//  Destructeur
// -----------------------

Target::~Target ( ) 
{ 
  delete locus_;
}

bool Target::hasLocus ( ) const 
{
  bool res = false ;
  if (locus_ != NULL )
  {
    res = true;
  }
  return res;
}

// -----------------------
//Accessor
// -----------------------

// Name
//
void Target::setName ( std::string name ) 
{
 	name_ = name;
}
std::string Target::getName ( )
{
	return name_;
}

// Locus
//
void Target::setLocus ( Locus * locus ) 
{
 	locus_ = locus;
}
Locus * Target::getLocus ( )
{
	return locus_;
}

// SequenceData_
//
void Target::setSequenceData ( std::string sequenceData ) 
{
 	sequenceData_ = sequenceData;
}
std::string Target::getSequenceData ( )
{
	return sequenceData_;
}


//targetLength_
//
void Target::setGlobalLength ( int targetLength )
{
	targetLength_=targetLength;
}
int Target::getGlobalLength ( )
{
	return targetLength_;
}


//isFullLength_
//
void Target::setIsFullLength ( int isFullLength )
{
	isFullLength_=isFullLength;
}
int Target::getIsFullLength ( )
{
	return isFullLength_;
}

int Target::getFrameHit ( )
{
  return frameHit_;
}

int Target::getScoreHit ( )
{
  return scoreHit_;
}


//getString
std::string Target::getString ( )  const 
{
	std::string my_target = "";
	if ( name_ != ""  && locus_ != NULL)
	{	
		 my_target = "Target="+name_+" "+locus_->getString()+";";
	}
	if ( isFullLength_ != UNKNOWN )
	{
		my_target += "is_full_length=" + to_string(isFullLength_) + ";";
	}
	if ( targetLength_ != -1 && targetLength_ != 0)
	{
		my_target += "target_length=" + to_string(targetLength_)+ ";";
	}
	if ( sequenceData_ != ""  )
	{
		my_target += "target_sequence=" + sequenceData_+ ";";
	}
	if ( frameHit_ != 0 )
	{
	  my_target += "frame_hit=" + to_string(frameHit_)+ ";";
	}
	if ( scoreHit_ > -1 )
	{
	  my_target += "score_hit=" + to_string(scoreHit_)+ ";";
	}
	return my_target;
}
