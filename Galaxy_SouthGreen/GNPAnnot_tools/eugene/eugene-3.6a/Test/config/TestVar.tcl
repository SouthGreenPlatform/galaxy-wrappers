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
# $Id: TestVar.tcl,v 1.24 2010-01-25 16:34:16 sallet Exp $
# ------------------------------------------------------------------
# File:     TestVar.tcl
# Contents: Environment variables definition
# ------------------------------------------------------------------


############################# Environment variables ###########################
set EUGENE_TEST_PAR eugeneTest.par
set EUGENE src/eugene

if {$action=="Test"} {
# ask tests in the root directory
    set EUGENE_DIR [exec pwd]
    set OUTPUT_DIR Test/Outputs
    set SEQ_DIR Test/Sequences
    set TRACE_DIR Test/TestTrace
} elseif {$action=="Generate"} {
# ask generation of ouput files in Test/config
    set EUGENE_DIR [exec pwd]/../..
    set OUTPUT_DIR ../Outputs
    set SEQ_DIR ../Sequences
    set TRACE_DIR ../TestTrace
}

############################################################################
set AllSensorsList {AnnotaStruct BlastX Est EuStop FrameShift GCPlot GFF GSplicer \
		    Homology IfElse MarkovConst MarkovIMM MarkovProt NG2 NStart PatConst  \
		    PepSignal Plotter Repeat Riken SMachine SPred SpliceWAM StartWAM  \
		    Tester Transcript NcRNA}

############################################################################
# At the moment, NcRNA  is not in this list because the gff3 format for this 
# sensor is already tested in unit tests
set AllGff3Sensors {AnnotaStruct BlastX Est GFF GSplicer NG2 NStart \
		    PepSignal Repeat SMachine SPred Homology}
############################################################################
##################### Units tests variables ################################
set SEQ(Sensor) {seq14ac002535g4g5}
set OPTIONS(Sensor) "-pd"

############################################################################# 
set FunctionalTestList {SeqAra SeqDoc SeqHom SeqAlt SeqRest SeqNcRNA}

##################### SeqAra test variables #################################
set SensorsList(SeqAra) {MarkovIMM MarkovConst EuStop NStart IfElse GSplicer Est BlastX}
set SEQ(SeqAra) {seq25ab005234g10g11.tfa}
set IMG(SeqAra) {seq25ab005234g10g11.000.png}
set FILE_REF(SeqAra) Output_SeqAra
set OPTIONS(SeqAra) "-po -g"

##################### SeqDoc test variables ############################
set SensorsList(SeqDoc) {MarkovIMM MarkovConst EuStop NStart IfElse Transcript \
		     FrameShift }
set SEQ(SeqDoc) {SYNO_ARATH.fasta}
set IMG(SeqDoc) {SYNO_ARATH.000.png}
set FILE_REF(SeqDoc) Output_SeqDoc
set OPTIONS(SeqDoc) "-po -g"

##################### SeqHom test variables ################################
set SensorsList(SeqHom) {MarkovProt MarkovConst EuStop StartWAM SpliceWAM Homology \
			     Transcript GFF}
set SEQ(SeqHom) {exSeqHom.fasta}
set IMG(SeqHom) {exSeqHom.000.png}
set FILE_REF(SeqHom) Output_SeqHom
set OPTIONS(SeqHom) "-po -g"
set NewValueSeqHom(Transcript.Start*) 1e-9
set NewValueSeqHom(Transcript.Stop*)  1e-9
set NewValueSeqHom(EuStop.stopP*)     3.1439

##################### SeqRest test variables #################################
set SensorsList(SeqRest) {MarkovIMM EuStop PatConst IfElse Riken Repeat}
set SEQ(SeqRest) {seq25ab005234g10g11.tfa}
set IMG(SeqRest) {seq25ab005234g10g11.000.png}
set FILE_REF(SeqRest) Output_SeqRest
set OPTIONS(SeqRest) "-po -g"


############################################################################
set ArabidopsisTestList {Araset Big}

##################### Araset test variables ################################
set SEQBASE_DIR "/Annotation/Arabidopsis/araset/Genes"
#"/home/bardou/araset/Genes"
set SensorsList(Araset) {MarkovIMM MarkovConst EuStop NStart IfElse GSplicer Est BlastX Transcript}
set SEQ(Araset) {}
set seq [eval exec ls $SEQBASE_DIR]
foreach f $seq {
    if { [string match *.tfa $f] } {
	set ff ${SEQBASE_DIR}/$f
	lappend SEQ(Araset) $ff
    }
}
set FILE_REF(Araset) Output_Araset
set OPTIONS(Araset) "-po"
set FILE_REF(ArasetSpSn) Output_Araset_SpSn 
set FILE_COORD(ArasetSpSn) "/Annotation/Arabidopsis/araset/Truth/araset.coord"
set PRG_EVAL_PRED $EUGENE_DIR/Procedures/Eval/egn_evalpred.pl

##################### Big test variables ################################
set SensorsList(Big) {MarkovIMM  Transcript MarkovConst EuStop NStart  GSplicer Est BlastX}
set SEQ(Big) $SEQ_DIR/chr1_2002.tfa_22540000-23040000.tfa
set FILE_REF(Big) Output_Big
set OPTIONS(Big) "-po"

##################### SeqAlt test variables #################################
set SensorsList(SeqAlt) {MarkovIMM MarkovConst EuStop NStart IfElse Est NG2 SPred}
set SEQ(SeqAlt) {At5g18830.fasta.genomicAJ011613.fasta}
set IMG(SeqAlt) {At5g18830.fasta.genomicAJ011613.000.png}
set FILE_REF(SeqAlt) Output_SeqAlt
set OPTIONS(SeqAlt) "-po -g -a"

##################### SeqNcRNA test variables ###############################
set SensorsList(SeqNcRNA) {MarkovIMM MarkovConst NStart Transcript NcRNA EuStop IfElse}
set SEQ(SeqNcRNA) {seq14ac002535g4g5.tfa}
set IMG(SeqNcRNA) {seq14ac002535g4g5.000.png}
set FILE_REF(SeqNcRNA) Output_SeqNcRNA
set OPTIONS(SeqNcRNA) "-po -g -a"
######################################################################################### Parameters optimization ##############################
set FILE_REF(Optimization) Output_Optimization


