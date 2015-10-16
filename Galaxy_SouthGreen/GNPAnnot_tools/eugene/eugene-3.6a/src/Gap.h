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

#ifndef GAP_H
#define GAP_H

#include <string>
#include "System.h"  //to_string()
/**
  * class Gap
  * describe a gap in protein and nucleotide alignment
  */

class Gap
{


private:

  // Private attributes
  //  
  char type_;
  int length_;


public:

   
  // Constructors/Destructors

   Gap ( );
   Gap (char type , int length );
   virtual ~Gap ( );


  // Public attribute accessor methods
  std::string getString ( ) const ;

  friend std::ostream& operator<<(std::ostream& o, const Gap & gap )
  {
	return o << gap.getString();
  }

};

#endif // GAP_H
