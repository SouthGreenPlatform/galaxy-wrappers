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
#include "Gap.h"

// -----------------------
//  Default constructeur
// -----------------------
Gap::Gap ( ) 
{
  type_ = ' ';
  length_ = 0;
}

Gap::Gap (char type, int length ) 
{
  type_ = type;
  length_ = length;
}

// -----------------------
//  Destructeur
// -----------------------

Gap::~Gap ( ) { }


//getString
std::string Gap::getString ( ) const 
{
	std::string my_gap = "";
	
	if ( type_ != ' ')
	{	
		 my_gap = to_string(type_)+to_string(length_);
	}
	return my_gap;
}
