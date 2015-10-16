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
# $Id: Makefile.am,v 1.16 2010-01-26 09:47:47 sallet Exp $
# ------------------------------------------------------------------
# File:     Makefile.am
# Contents: compilation of the src directory
# ------------------------------------------------------------------


SUBDIRS = Parametrization GDIF . SensorPlugins

AM_CXXFLAGS = $(eugene_cxxflags) -DDEFAULT_EUGENE_DIR=\"${pkgdatadir}\"
AM_CFLAGS = 

bin_PROGRAMS = eugene

eugene_SOURCES = AltEst.cc AltEst.h Attributes.cc Attributes.h BackP.cc BackP.h \
	Const.h DAG.cc DAG.h DNASeq.cc DNASeq.h Dll.cc Dll.h EuGene.cc Gap.cc Gap.h \
	GeneFeature.cc GeneFeature.h GeneFeatureSet.cc GeneFeatureSet.h Gff3Line.cc Gff3Line.h \
	Hits.cc Hits.h Locus.cc Locus.h MSensor.cc MSensor.h Param.cc Param.h \
	PenaltyDist.cc PenaltyDist.h Prediction.cc Prediction.h Prediction_cte.h Range.cc \
	Range.h Sensor.cc Sensor.h SensorIF.h State.h State.cc SoTerms.cc SoTerms.h System.cc \
	System.h Target.cc Target.h getopt.h

eugene_LDADD = Parametrization/libParametrization.a  GDIF/libGDIF.a -ldl

CLEANFILES = *~
