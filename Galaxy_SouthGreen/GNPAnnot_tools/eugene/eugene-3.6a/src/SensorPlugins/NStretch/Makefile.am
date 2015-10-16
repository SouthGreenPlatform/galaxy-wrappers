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
# $Id: Makefile.am,v 1.1 2009-04-22 13:53:01 tschiex Exp $
# ------------------------------------------------------------------
# File:     Makefile.am
# Contents: compilation of the NStretch directory
# ------------------------------------------------------------------

AM_CXXFLAGS = $(eugene_cxxflags)
AM_CFLAGS = 


noinst_PROGRAMS = Sensor.NStretch.so
Sensor_NStretch_so_SOURCES = Sensor.NStretch.cc Sensor.NStretch.h 
Sensor_NStretch_so_LDFLAGS = $(sensors_ldflags) 


CLEANFILES = *~
