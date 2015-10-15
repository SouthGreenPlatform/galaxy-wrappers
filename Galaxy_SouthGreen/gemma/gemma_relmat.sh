#!/usr/bin/env bash

# file hosted online: https://github.com/SouthGreenPlatform/galaxy-wrappers

relmatfile=$1
logfile=$2
outprefix=`basename $3`
shift 3

SHDIR=`dirname $0`
TMPDIR=`mktemp -d ${SHDIR}/gemma.XXXXXXXXXX`
cd $TMPDIR

${SHDIR}/gemma $* -o ${outprefix}

mv output/${outprefix}.cXX.txt ${relmatfile}
mv output/${outprefix}.log.txt ${logfile}

rm -rf $TMPDIR
