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
// $Id: EndianConv.h,v 1.1 2004-09-16 13:13:05 cros Exp $
// ------------------------------------------------------------------
// File:     EndianConv.h
// Contents: utilitary functions
// ------------------------------------------------------------------

#ifndef  ENDIAN_CONV_H_INCLUDED
#define ENDIAN_CONV_H_INCLUDED

inline unsigned int LEndianReverse (unsigned int N)
{
  return ((N & 0x000000FF) << 24) |
         ((N & 0x0000FF00) << 8)  |
         ((N & 0x00FF0000) >> 8)  |
         ((N & 0xFF000000) >> 24);
}

inline unsigned short int SEndianReverse (unsigned short int N)
{
  return ((N & 0x00FF) << 8) |
         ((N & 0xFF00) >> 8);
}

#endif

