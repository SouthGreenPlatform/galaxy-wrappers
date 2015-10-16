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
# $Id: configure.ac,v 1.31 2010-01-26 08:30:25 sallet Exp $
# ------------------------------------------------------------------
# File:     configure.ac
# Contents: compilation configuration
# ------------------------------------------------------------------

# Autoconf requirements
AC_INIT(, , [EuGene Team eugene@ossau.toulouse.inra.fr])
AC_CONFIG_AUX_DIR(config)
AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE(eugene, 3.6a)

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_AWK
AC_PROG_LN_S 
AC_PROG_RANLIB


# Checks for libraries.
AC_CHECK_LIB([c], [main], [], [ \
	echo "Error! You need to have the c library to compile eugene."
        exit -1
        ])
AC_CHECK_LIB([dl], [main], [], [ \
	echo "Error! You need to have the dl library to compile eugene."
        exit -1
        ])
AC_CHECK_LIB([gcc], [main], [], [ \
	echo "Error! You need to have the gcc library to compile eugene."
        exit -1
        ])
AC_CHECK_LIB([gd], [main], [], [ \
	echo "Error! You need to have the gd library (>= 2.0) to compile eugene."
	echo "See http://www.boutell.com/gd/http/ for a download."
        exit -1
        ])
AC_CHECK_LIB([m], [main], [], [ \
	echo "Error! You need to have the m library to compile eugene."
        exit -1
        ])
AC_CHECK_LIB([png], [main], [], [ \
	echo "Error! You need to have the png library for the gd library."
        exit -1
        ])
AC_CHECK_LIB([z], [main], [], [ \
	echo "Error! You need to have the z library for the gd library."
        exit -1
        ])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([netdb.h stdlib.h string.h strings.h unistd.h values.h getopt.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_STRUCT_TM
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_STAT
AC_CHECK_FUNCS([floor pow sqrt strchr strspn uname getopt])

#check for OS
AM_CONDITIONAL(SOLARIS_OS, test "$OSTYPE" = solaris)
AM_CONDITIONAL(MACOS_OS, test "$OSTYPE" = darwin)


AC_CONFIG_FILES([ 
	Makefile
	README
	src/Makefile
	src/GDIF/Makefile
	src/Parametrization/Makefile
	src/SensorPlugins/Makefile
	src/SensorPlugins/0_SensorTk/Makefile
	src/SensorPlugins/0_SensorTk/GetData/Makefile
	src/SensorPlugins/AnnotaStruct/Makefile
	src/SensorPlugins/BlastX/Makefile
	src/SensorPlugins/Est/Makefile
	src/SensorPlugins/EuStop/Makefile
	src/SensorPlugins/FrameShift/Makefile
	src/SensorPlugins/GCPlot/Makefile
	src/SensorPlugins/GFF/Makefile
	src/SensorPlugins/GSplicer/Makefile
	src/SensorPlugins/Homology/Makefile
	src/SensorPlugins/IfElse/Makefile
	src/SensorPlugins/MarkovConst/Makefile
	src/SensorPlugins/MarkovIMM/Makefile
	src/SensorPlugins/MarkovIMM/GetData/Makefile
	src/SensorPlugins/MarkovProt/Makefile
	src/SensorPlugins/MarkovProt/GetData/Makefile
	src/SensorPlugins/NG2/Makefile
	src/SensorPlugins/NStart/Makefile
	src/SensorPlugins/PatConst/Makefile
	src/SensorPlugins/PepSignal/Makefile
	src/SensorPlugins/Plotter/Makefile
	src/SensorPlugins/Repeat/Makefile
	src/SensorPlugins/Riken/Makefile
	src/SensorPlugins/SMachine/Makefile
	src/SensorPlugins/SPred/Makefile
	src/SensorPlugins/SpliceWAM/Makefile
	src/SensorPlugins/StartWAM/Makefile
	src/SensorPlugins/Tester/Makefile
	src/SensorPlugins/Transcript/Makefile
	src/SensorPlugins/QualData/Makefile
	src/SensorPlugins/NStretch/Makefile
	src/SensorPlugins/NcRNA/Makefile
	plugins/Makefile
	doc/Makefile
])

AC_OUTPUT
