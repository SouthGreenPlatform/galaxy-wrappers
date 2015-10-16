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
# $Id: Makefile.am,v 1.3 2005-03-17 15:02:23 cros Exp $
# ------------------------------------------------------------------
# File:     Makefile.am
# Contents: compilation of the Riken directory
# ------------------------------------------------------------------

AM_CXXFLAGS = $(eugene_cxxflags)
AM_CFLAGS = 


noinst_PROGRAMS = Sensor.Riken.so
Sensor_Riken_so_SOURCES = Sensor.Riken.cc  Sensor.Riken.h
Sensor_Riken_so_LDFLAGS = $(sensors_ldflags) 


CLEANFILES = *~
