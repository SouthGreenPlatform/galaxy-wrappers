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

#ifndef TARGET_H
#define TARGET_H

#include <string>
#include "Locus.h"

/**
  * class Target
  * define the hit targeted of an alignment
  */

class Target
{



private:

  // Private attributes
  //  
  std::string name_;
  Locus * locus_;
  std::string sequenceData_;
  int targetLength_;
  int isFullLength_; 
  // 0 -> Est or fragment protein
  // 1 -> Riken (have the 5' and 3' Est)
  // 2 -> cDNA or Protein full length
  int frameHit_;
  int scoreHit_;

public:

  enum LengthType { UNKNOWN = -1, NOT_FULL_LENGTH = 0, RIKEN = 1, FULL_LENGTH=2 };
  // Constructors/Destructors

   Target ( );
   Target ( std::string name, std::string sequenceData,  int targetLength, int isFullLength, Locus * locus,int frameHit,int scoreHit);
   Target ( std::string name, std::string sequenceData,  int targetLength, int isFullLength, int start, int end, char strand, int frameHit, int scoreHit);
   Target ( Target * target);
   virtual ~Target ( );


  // Public attribute accessor methods

  void setName ( std::string name );
  std::string getName ( );

  void setLocus ( Locus * range );
  Locus* getLocus ( );

  void setSequenceData ( std::string sequenceData );
  std::string getSequenceData ( );

  void setGlobalLength ( int targetLength );
  int getGlobalLength ( );

  void setIsFullLength ( int isFullLength );
  int getIsFullLength ( );
  int getFrameHit ( );
  int getScoreHit ( );
  bool hasLocus ( ) const;
  std::string getString ( ) const ;
  friend std::ostream& operator<<(std::ostream& o, const Target & target )
  {
	return o << target.getString();
  }

};

#endif // TARGET_H
