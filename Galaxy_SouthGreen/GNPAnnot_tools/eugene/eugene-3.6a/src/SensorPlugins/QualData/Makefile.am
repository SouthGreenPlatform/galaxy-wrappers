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
# $Id: Makefile.am,v 1.1 2008-09-22 12:54:20 dleroux Exp $
# ------------------------------------------------------------------
# File:     Makefile.am
# Contents: compilation of the QualData directory
# ------------------------------------------------------------------

AM_CXXFLAGS = $(eugene_cxxflags)
AM_CFLAGS = 


noinst_PROGRAMS = Sensor.QualData.so
Sensor_QualData_so_SOURCES = Sensor.QualData.cc Sensor.QualData.h 
Sensor_QualData_so_LDFLAGS = $(sensors_ldflags)


CLEANFILES = *~
