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
# $Id: Makefile.am,v 1.2 2005-03-17 14:55:53 cros Exp $
# ------------------------------------------------------------------
# File:     Makefile.am
# Contents: compilation of the SPred directory
# ------------------------------------------------------------------

AM_CXXFLAGS = $(eugene_cxxflags)
AM_CFLAGS = 


noinst_PROGRAMS = Sensor.SPred.so
Sensor_SPred_so_SOURCES = Sensor.SPred.cc Sensor.SPred.h
Sensor_SPred_so_LDFLAGS = $(sensors_ldflags) 


CLEANFILES = *~
