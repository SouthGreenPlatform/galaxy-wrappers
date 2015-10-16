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
#ifndef RANGE_H
#define RANGE_H
#include <iostream>
#include <string>
#include <iostream>
    using namespace std; 
#include "System.h"  //to_string()
/**
  * class Locus
  */


class Locus
{
private:

  // Private attributes
  //  
  int start_ ;
  int end_;
  char strand_;

public:

  // Constructors/Destructors
  //  

  Locus ( );
  Locus (int start, int end, char strand);
  Locus (const Locus & locus );
  virtual ~Locus ( );

  // Puplic attribute accessor methods
  //  

  //start
  void setStart ( int start );
  int getStart ( ) const;

  void setEnd ( int end );
  int getEnd ( ) const;

  void setStrand ( char strand );
  char getStrand ( ) const;
  std::string getString ( ) const;

  friend std::ostream& operator<<(std::ostream& o, const Locus & range )
  {
	return o << range.getString();
  }
};




#endif // RANGE_H
