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
#include "Range.h"

// -----------------------
//  Default constructeur
// -----------------------
Range::Range ( ) 
{
  start_ =0;
  end_   =0;
  strand_=' ';
}

Range::Range (int start, int end, char strand) 
{
  start_ =start;
  end_   =end;
  strand_=strand;
}

Range::Range (Range *range )
{
  start_ =range->getStart();
  end_   =range->getEnd();
  strand_=range->getStrand();
}

// -----------------------
//  Destructeur
// -----------------------
Range::~Range ( ) { }


//getString
std::string Range::getString ( )
{
	std::string my_range = to_string(start_)+' '+to_string(end_);
	if ( strand_ != ' ' )
	{
		my_range +=' '+to_string(strand_);
	}
	return my_range;
}

// surcharge <<
// ostream& operator<<(ostream &o, const Range *range )
// {
// 	return o << range->getString();
// }

// -----------------------
//  Start
// -----------------------
void Range::setStart ( int start )
{
	start_=start;
}
int Range::getStart ( )
{
	return start_;
}

// -----------------------
//  End
// -----------------------
void Range::setEnd ( int end )
{
	end_=end;
}
int Range::getEnd ( )
{
	return end_;
}


// -----------------------
// Strand
// -----------------------
void Range::setStrand ( char strand )
{
	strand_=strand;
}
char Range::getStrand ( )
{
	return strand_;
}


