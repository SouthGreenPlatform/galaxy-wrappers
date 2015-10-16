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
# $Id: Makefile.am,v 1.2 2005-03-17 14:22:26 cros Exp $
# ------------------------------------------------------------------
# File:     Makefile.am
# Contents: compilation of the Parametrization directory
# ------------------------------------------------------------------

AM_CXXFLAGS = $(eugene_cxxflags)

noinst_LIBRARIES = libParametrization.a
libParametrization_a_SOURCES = ParaOptimization.cc ParaOptimization.h \
			OptiAlgorithm.cc OptiAlgorithm.h \
			Random.cc Random.h \
			LineSearch.cc LineSearch.h \
			Genetic.cc Genetic.h \
			CrossingOver.cc CrossingOver.h \
			Mutating.cc Mutating.h \
			Scaling.cc Scaling.h \
			Selecting.cc Selecting.h \
			Sharing.cc Sharing.h \
			SimulatedAnnealing.cc SimulatedAnnealing.h \
			Chromosome.cc Chromosome.h


CLEANFILES = *~
