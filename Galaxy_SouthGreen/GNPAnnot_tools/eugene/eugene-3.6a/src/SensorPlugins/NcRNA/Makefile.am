# ------------------------------------------------------------------
# Copyright (C) 2004 INRA <eugene@ossau.toulouse.inra.fr>
#
# This program is open source; you can redistribute it and/or modify
# it under the terms of the Artistic License (see LICENSE file).
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#
# You should have received a copy of Artistic License along with
# this program; if not, please see http://www.opensource.org
#
# $Id: Makefile.am,v 1.2 2010-01-25 16:53:52 sallet Exp $
# ------------------------------------------------------------------
# File:     Makefile.am
# Contents: compilation of the NcRNA directory
# ------------------------------------------------------------------

AM_CXXFLAGS = $(eugene_cxxflags)
AM_CFLAGS = 


noinst_PROGRAMS = Sensor.NcRNA.so
Sensor_NcRNA_so_SOURCES = Sensor.NcRNA.cc Sensor.NcRNA.h 
Sensor_NcRNA_so_LDFLAGS = $(sensors_ldflags) 


CLEANFILES = *~
