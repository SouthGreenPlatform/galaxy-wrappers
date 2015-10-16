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
# $Id: TestProc.tcl,v 1.32 2010-02-01 16:18:43 sallet Exp $
# ------------------------------------------------------------------
# File:     TestProc.tcl
# Contents: Procedures testing the eugene software
# ------------------------------------------------------------------




############################################################################
# Procedure ModifyParaValue
# Description : in the parameter file specified by the 1st argument, 
#               modify the value of the parameters given as index 
#               of the 2sd argument (array)
#               which contains the new values of parameters
# Example     : tcl> set V(EuGene.minEx) 4
#               tcl> set V(Output.graph) 1
#               tcl> ModifyParaValue eugeneTest.par V
#               After execution, in the file eugeneTest.par 
#               the line 'EuGene.minEx   3' will be ' EuGene.minEx   4' 
#               the line 'Output.graph   0' will be 'Output.graph  1'   
# Note that if there is no ambiguity, it is sufficient to give 
# the end of parameter name
# In the example, set V(minEx) 4 would produces the same result 
# if there is not an other 'minEx' in the file.
# BEWARE in case of parameter name included in an other: to distinguish
# the shortest it is necessary to add a space at the end of its name
############################################################################
proc ModifyParaValue {FileName NewValues} {
# To pass an array as 2nd argument
    upvar $NewValues V

    set f [open $FileName r]
    set new_content ""

    foreach line [split [read $f] \n] {
	if {[string length $line] != 0} {
	    set newline ""
	    foreach para [array names V] {
		set Pos [string first ${para} $line]
		if {$Pos != -1} {
		    set newline [string range $line \
				     0 [expr $Pos - 1 + [string length $para]]]
		    set newline "$newline $V($para)"
		}
	    }
	    if {[string length $newline] == 0} {
		set newline $line
	    }
	    set new_content "$new_content$newline\n"
	}
    }
    close $f

# Write the new parameter file
    set f [open $FileName w]
    puts -nonewline $f $new_content
    close $f
}


##############################################################################
# Procedure   : InitParameterFile
# Description : Update values of parameters in a parameter file
#               Arguments: FileName = name of the parameter file
#                          SensorsList = list of sensors that will be specified 
#                                        not to use
#                          EuGeneDir = directory of the executable
# BEWARE      : Priority of sensors are not updated
###############################################################################
proc InitParameterFile {FileName SensorsList EuGeneDir} {

# when name of parameter is include in an other
# (example SpliceConst.accP is included in SpliceConst.accPNo)
# it is necessary to add a space at the name of the sensor
# (example SpliceConst.accP$space)
set space " "

# Set the parameter to wanted values
#################################################################
###################### GENERAL PARAMETERS #######################
#################################################################
# Avoid to set set EuGene.version
set NewValue1(EuGene.organism) 		Arabidopsis
set NewValue1(EuGene.sloppy)            0
##### Lengths #####
set NewValue1(EuGene.InitExDist)	init.dist
set NewValue1(EuGene.IntrExDist)	intr.dist
set NewValue1(EuGene.TermExDist)	term.dist
set NewValue1(EuGene.SnglExDist)	sngl.dist
set NewValue1(EuGene.IntronDist)	intron.dist
set NewValue1(EuGene.InterGDist)	intergenic.dist
set NewValue1(EuGene.5PrimeDist)	utr.dist
set NewValue1(EuGene.3PrimeDist)	utr.dist
set NewValue1(EuGene.RnaDist)	    rna.dist
##### Priors #####
set NewValue1(EuGene.SplicedStopPen)    1e999.0
set NewValue1(EuGene.ExonPrior)	        0.18	
set NewValue1(EuGene.IntronPrior)	0.17	
set NewValue1(EuGene.InterPrior)	0.4
set NewValue1(EuGene.FivePrimePrior)	0.03
set NewValue1(EuGene.ThreePrimePrior)	0.07
set NewValue1(EuGene.RnaPrior)	0.15
##### Output control ######
set NewValue1(Output.MinCDSLen)		60
set NewValue1(Output.truncate)		5
set NewValue1(Output.stepid)		1
set NewValue1(Output.graph)		0	
set NewValue1(Output.resx)		900
set NewValue1(Output.resy)		400
set NewValue1(Output.glen)		-1
set NewValue1(Output.golap)		-1
set NewValue1(Output.gfrom)		-1
set NewValue1(Output.gto)		-1
set NewValue1(Output.window)		48
set NewValue1(Output.format)		l
set NewValue1(Output.offset)		0
set NewValue1(Output.normopt)		1
set NewValue1(Output.Prefix)		./
set NewValue1(Output.webdir)            LOCAL
set NewValue1(Output.intron)		0
#################################################################
################### SIGNAL SENSORS PARAMETERS####################
#################################################################
##### EuStop parameters #####
set NewValue1(EuStop.stopP*)            4.155
##### FrameShift parameters #####
set NewValue1(FrameShift.Ins*)	1e999.0
set NewValue1(FrameShift.Del*)	1e999.0
##### GeneSplicer parameters #####
set NewValue1(GSplicer.coefAcc*)	1
set NewValue1(GSplicer.penAcc*)	        0
set NewValue1(GSplicer.coefDon*)	1
set NewValue1(GSplicer.penDon*) 	0
##### PepSignal #####
set NewValue1(PepSignal.startP*)	1
set NewValue1(PepSignal.startB*)	0
##### SpliceMachine parameters #####
set NewValue1(SMachine.cmd)		"splicemachine.pl "
set NewValue1(SMachine.accP*)		0.102032725565
set NewValue1(SMachine.accB*)		5.585
set NewValue1(SMachine.donP*)		0.020202707318
set NewValue1(SMachine.donB*)		27.670
set NewValue1(SMachine.startP*)	        0.052
set NewValue1(SMachine.startB*)	        0.308
##### NetGene2 parameters #####
set NewValue1(NG2.accP*\[0\])     0.903
set NewValue1(NG2.accB*\[0\])     5.585
set NewValue1(NG2.donP*\[0\])     0.980
set NewValue1(NG2.donB*\[0\])     27.670
set NewValue1(NG2.accP*\[1\])	0.903
set NewValue1(NG2.accB*\[1\])	5.585
set NewValue1(NG2.donP*\[1\])	0.980
set NewValue1(NG2.donB*\[1\])	27.670
##### NetStart parameters #####
set NewValue1(NStart.startP*)	0.052
set NewValue1(NStart.startB*)	0.308
##### PatConst sensor (uniform penalties) #####
set NewValue1(PatConst.type\[0\])	start
set NewValue1(PatConst.pat\[0\])	ATG
set NewValue1(PatConst.newStatePos\[0\]) 1
set NewValue1(PatConst.patP*\[0\])	2.897949
set NewValue1(PatConst.patPNo*\[0\])	0
##### Sensor SpliceWAM #####
set NewValue1(SpliceWAM.MarkovianOrder)	1
set NewValue1(SpliceWAM.donmodelfilename)	WAM/WAM.ARA.DON.L9
set NewValue1(SpliceWAM.NbNtBeforeGT)	3
set NewValue1(SpliceWAM.NbNtAfterGT)	4
set NewValue1(SpliceWAM.DonScaleCoef*)	2.9004
set NewValue1(SpliceWAM.DonScalePenalty*)	-7.5877
set NewValue1(SpliceWAM.accmodelfilename)	WAM/WAM.ARA.ACC.L7
set NewValue1(SpliceWAM.NbNtBeforeAG)		2
set NewValue1(SpliceWAM.NbNtAfterAG)		1
set NewValue1(SpliceWAM.AccScaleCoef*)		2.9004
set NewValue1(SpliceWAM.AccScalePenalty*)		-7.5877
##### SplicePredictor parameters #####
set NewValue1(SPred.accP*\[0\])   0.987
set NewValue1(SPred.accB*\[0\])  3.850
set NewValue1(SPred.donP*\[0\])   0.929
set NewValue1(SPred.donB*\[0\])   10.800
set NewValue1(SPred.accP*\[1\])	0.987
set NewValue1(SPred.accB*\[1\])	3.850
set NewValue1(SPred.donP*\[1\])	0.929
set NewValue1(SPred.donB*\[1\])	10.800
##### Sensor StartWAM #####
set NewValue1(StartWAM.modelfilename)	WAM/WAM.ARA.START9
set NewValue1(StartWAM.NbNtBeforeATG)	3
set NewValue1(StartWAM.NbNtAfterATG)	3
set NewValue1(StartWAM.MarkovianOrder)		1
set NewValue1(StartWAM.ScaleCoef*)		0.1594
set NewValue1(StartWAM.ScalePenalty*)		-3.1439
##### Transcript parameters #####
set NewValue1(Transcript.Start*)	4.155
set NewValue1(Transcript.Stop*)		4.155
#################################################################
################# CONTENT SENSORS PARAMETERS ####################
#################################################################
##### Proteic similarity sensor parameters #####
set NewValue1(BlastX.PostProcess) 0
set NewValue1(BlastX.PPNumber) 5
set NewValue1(BlastX.levels)	0
set NewValue1(BlastX.level0*)	0.2
set NewValue1(BlastX.level1*)	0.05
set NewValue1(BlastX.level2*)	0.0
set NewValue1(BlastX.level3*)	0.0
set NewValue1(BlastX.level4*)	0.0
set NewValue1(BlastX.level5*)	0.0
set NewValue1(BlastX.level6*)	0.0
set NewValue1(BlastX.level7*)	0.0
set NewValue1(BlastX.level8*)	0.0
set NewValue1(BlastX.level9*)	0.0
set NewValue1(BlastX.blastxM*)	10	
set NewValue1(BlastX.minIn) 50
##### Est sensor parameters #####
set NewValue1(Est.PostProcess\[0\])	0
set NewValue1(Est.PPNumber\[0\])     2
set NewValue1(Est.estP*\[0\])	-0.4
set NewValue1(Est.estM\[0\])	        6
set NewValue1(Est.utrP*\[0\])	0.35
set NewValue1(Est.utrM\[0\])	        5
set NewValue1(Est.SpliceBoost*\[0\]) 0.0
set NewValue1(Est.StrongDonor\[0\])	0.95
set NewValue1(Est.MinDangling\[0\])  10
set NewValue1(Est.MaxIntron\[0\])    15000
set NewValue1(Est.FileExtension\[0\]) .est
set NewValue1(Est.PostProcess\[1\])	0
set NewValue1(Est.PPNumber\[1\])     2
set NewValue1(Est.estP*\[1\])	-0.4
set NewValue1(Est.estM\[1\])	        6
set NewValue1(Est.utrP*\[1\])	0.35
set NewValue1(Est.utrM\[\])	        5
set NewValue1(Est.SpliceBoost*\[1\]) 0.0
set NewValue1(Est.StrongDonor\[1\])	0.95
set NewValue1(Est.MinDangling\[1\])  10
set NewValue1(Est.MaxIntron\[1\])    15000
set NewValue1(Est.FileExtension\[1\]) .est2
##### Homology Sensor parameters #####
set NewValue1(Homology.TblastxP*\[0\]) 	0
set NewValue1(Homology.TblastxB*\[0\]) 	0.0595
set NewValue1(Homology.protmatname\[0\])	BLOSUM80
set NewValue1(Homology.MaxHitLen\[0\])	15000
set NewValue1(Homology.FileExtension\[0\]) .tblastx

set NewValue1(Homology.TblastxP*\[1\]) 	0
set NewValue1(Homology.TblastxB*\[1\]) 	0.0595
set NewValue1(Homology.protmatname\[1\])	BLOSUM80
set NewValue1(Homology.MaxHitLen\[1\])	15000
set NewValue1(Homology.FileExtension\[1\]) .tblastx2

##### State penalties (exp length distributions) #####
set NewValue1(MarkovConst.minGC\[0\])	0
set NewValue1(MarkovConst.maxGC\[0\])	100
set NewValue1(MarkovConst.Coding*) 	1.0
set NewValue1(MarkovConst.IntronUTR*)	0.98
set NewValue1(MarkovConst.Intron*) 	1.0
set NewValue1(MarkovConst.UTR5*)	0.999
set NewValue1(MarkovConst.UTR3*) 	0.999
set NewValue1(MarkovConst.Inter*) 	1.0
##### Interpolated Markov Models parameters #####
set NewValue1(MarkovIMM.matname\[0\])	Ara2UTR.mat
set NewValue1(MarkovIMM.minGC\[0\])	0
set NewValue1(MarkovIMM.maxGC\[0\])	100
set NewValue1(MarkovIMM.useM0asIG\[0\])	0
set NewValue1(MarkovIMM.maxOrder\[0\])	8
##### Markov proteic model parameters #####
set NewValue1(MarkovProt.matname\[0\])	swissprot.maxorder2.bin
set NewValue1(MarkovProt.minGC\[0\])	0
set NewValue1(MarkovProt.maxGC\[0\])	100
set NewValue1(MarkovProt.maxorder)    2
set NewValue1(MarkovProt.order)        2
##### Repeat sensor parameters #####
set NewValue1(Repeat.UTRPenalty*)	0.0
set NewValue1(Repeat.IntronPenalty*)	0.1
set NewValue1(Repeat.ExonPenalty*)	1.0
#### NcRNA sensor parameters ####
set NewValue1(NcRNA.FileExtension) ncrna
set NewValue1(NcRNA.NpcRna*)   1
set NewValue1(NcRNA.TStartNpc*)  1
set NewValue1(NcRNA.TStopNpc*)  1
set NewValue1(NcRNA.Format) GFF3
#################################################################
############## SIGNAL/CONTENT SENSORS PARAMETERS ################
#################################################################
##### Sensors AnnotaStruct #####
set NewValue1(AnnotaStruct.FileExtension)     gff
set NewValue1(AnnotaStruct.Start*)            0.1
set NewValue1(AnnotaStruct.StartType)         p 
set NewValue1(AnnotaStruct.Stop*)             0.2
set NewValue1(AnnotaStruct.StopType)          p
set NewValue1(AnnotaStruct.Acc*)              0.3
set NewValue1(AnnotaStruct.AccType)           p
set NewValue1(AnnotaStruct.Don*)              0.4
set NewValue1(AnnotaStruct.DonType)           p
set NewValue1(AnnotaStruct.TrStart*)          0.5
set NewValue1(AnnotaStruct.TrStartType)       p
set NewValue1(AnnotaStruct.TrStop*)           0.6
set NewValue1(AnnotaStruct.TrStopType)        p
set NewValue1(AnnotaStruct.Exon*)             1
set NewValue1(AnnotaStruct.Intron*)           2
set NewValue1(AnnotaStruct.CDS*)              3
##### IfElse #####
set NewValue1(IfElse.SensorIf)		NG2
set NewValue1(IfElse.SensorElse)	SPred
##### Riken sensor parameters #####
set NewValue1(Riken.StrandRespect)		0
set NewValue1(Riken.Min\_est\_diff)		100
set NewValue1(Riken.Max\_overlap		60
set NewValue1(Riken.Max\_riken\_length)		60000
set NewValue1(Riken.Max\_riken\_est\_length)	3000
set NewValue1(Riken.Min\_riken\_length)		120 
set NewValue1(Riken.Min\_riken\_est\_length)	10
set NewValue1(Riken.RAFLPenalty*)		-120.0
#################################################################
################## OTHERS SENSORS PARAMETERS ####################
#################################################################
##### Sensor GCPlot #####
set NewValue1(GCPlot.Color)	5
set NewValue1(GCPlot.Zoom)	2.0
set NewValue1(GCPlot.Zoom3)	1.0
set NewValue1(GCPlot.Up)	GC
set NewValue1(GCPlot.Over)	ATGC
set NewValue1(GCPlot.Smooth)	100
##### GFF sensor parameters #####
set NewValue1(GFF.PostProcess)		0
##### Sensor Plotter
set NewValue1(Plotter.GC\[0\])	 1 
set NewValue1(Plotter.GC3\[0\])	 1 
set NewValue1(Plotter.A|T/A+T\[0\])	 0 
##### Sensor Tester #####
set NewValue1(Tester.Make)		SPSN
set NewValue1(Tester.Sensor)		EuStop
set NewValue1(Tester.Sensor.Instance)	0
set NewValue1(Tester.SPSN.MinNumbers)	100
set NewValue1(Tester.SPSN.Eval)         STOP
#################################################################
################# SENSORS CONFIGURATION PARAMETERS ##############
#################################################################
##### Sensors desactivation #####
# SIGNAL SENSORS
set NewValue1(Sensor.EuStop.use)	0
set NewValue1(Sensor.FrameShift.use)	0
set NewValue1(Sensor.GSplicer.use)	0
set NewValue1(Sensor.SMachine.use)	0
set NewValue1(Sensor.NG2.use)		0
set NewValue1(Sensor.NStart.use)	0
set NewValue1(Sensor.PatConst.use)	0
set NewValue1(Sensor.PepSignal.use)	0
set NewValue1(Sensor.SpliceWAM.use) 	0
set NewValue1(Sensor.SPred.use)	        0
set NewValue1(Sensor.StartWAM.use)	0
set NewValue1(Sensor.Transcript.use)	0
# CONTENT SENSORS
set NewValue1(Sensor.BlastX.use)	0
set NewValue1(Sensor.Est.use)		0
set NewValue1(Sensor.Homology.use)	0
set NewValue1(Sensor.MarkovConst.use)	0
set NewValue1(Sensor.MarkovIMM.use)	0
set NewValue1(Sensor.MarkovProt.use)	0
set NewValue1(Sensor.Repeat.use)	0
set NewValue1(Sensor.NStretch.use) 0
set NewValue1(Sensor.NcRNA.use) 0
# SIGNAL/CONTENT SENSORS
set NewValue1(Sensor.AnnotaStruct.use)  0
set NewValue1(Sensor.IfElse.use)	0
set NewValue1(Sensor.Riken.use)	        0
# OTHERS SENSORS
set NewValue1(Sensor.GCPlot.use)	0
set NewValue1(Sensor.GFF.use)		0
set NewValue1(Sensor.Plotter.use)	0
set NewValue1(Sensor.Tester.use)	0
#
##### Sensor priorities	 #####
# SIGNAL SENSORS
set NewValue1(Sensor.EuStop$space)	1
set NewValue1(Sensor.FrameShift$space)	1
set NewValue1(Sensor.GSplicer$space)	1
set NewValue1(Sensor.NG2$space)		1
set NewValue1(Sensor.NStart$space)	1
set NewValue1(Sensor.PatConst$space)	1
set NewValue1(Sensor.PepSignal$space) 	1
set NewValue1(Sensor.SMachine$space)	1
set NewValue1(Sensor.SpliceWAM$space) 	1
set NewValue1(Sensor.SPred$space)	1
set NewValue1(Sensor.StartWAM$space)	1
set NewValue1(Sensor.Transcript$space)	1
# CONTENT SENSORS
set NewValue1(Sensor.BlastX$space)	1
set NewValue1(Sensor.Est$space)		20
set NewValue1(Sensor.Homology$space)  	1
set NewValue1(Sensor.MarkovConst$space) 1
set NewValue1(Sensor.MarkovIMM$space) 	1
set NewValue1(Sensor.MarkovProt$space)	1
set NewValue1(Sensor.Repeat$space)	1
set NewValue1(Sensor.NStretch$space)   1
set NewValue1(Sensor.NcRNA$space)  1
# SIGNAL/CONTENT SENSORS
set NewValue1(Sensor.AnnotaStruct$space) 1
set NewValue1(Sensor.IfElse$space)	1
set NewValue1(Sensor.Riken$space)	1
# OTHERS SENSORS
set NewValue1(Sensor.GCPlot$space)	1
set NewValue1(Sensor.GFF$space)		1
set NewValue1(Sensor.Plotter$space)	1
set NewValue1(Sensor.Tester$space)	1
#################################################################
################### PARAMETERS OPTIMIZATION #####################
#################################################################
set NewValue1(ParaOptimization.Use)	        0
set NewValue1(ParaOptimization.TrueCoordFile) 	---
set NewValue1(ParaOptimization.EvalPredDir)     ../Procedures/Eval
set NewValue1(ParaOptimization.Algorithm)	GENETIC+LINESEARCH
set NewValue1(ParaOptimization.Test)	        1
set NewValue1(ParaOptimization.Trace)		1
#
set NewValue1(ParaOptimization.NbCluster) 3
set NewValue1(ParaOptimization.Cluster\[0\]) LINKED
set NewValue1(ParaOptimization.Cluster\[1\]) IDENTICAL
set NewValue1(ParaOptimization.Cluster\[2\]) IDENTICAL
#
set NewValue1(ParaOptimization.NbParameter)   	5
#
set NewValue1(ParaOptimization.Para.Name\[0\])	para1*
set NewValue1(ParaOptimization.Para.Max\[0\])	1	
set NewValue1(ParaOptimization.Para.Min\[0\])	0
set NewValue1(ParaOptimization.Para.Cluster\[0\]) 0
#
set NewValue1(ParaOptimization.Para.Name\[1\])	para2*
set NewValue1(ParaOptimization.Para.Max\[1\])	1	
set NewValue1(ParaOptimization.Para.Min\[1\])	0
set NewValue1(ParaOptimization.Para.Cluster\[1\]) 0
#
set NewValue1(ParaOptimization.Para.Name\[2\])	para3*
set NewValue1(ParaOptimization.Para.Max\[2\])	1	
set NewValue1(ParaOptimization.Para.Min\[2\])	0
set NewValue1(ParaOptimization.Para.Cluster\[2\]) 1
#
set NewValue1(ParaOptimization.Para.Name\[3\])	para4*
set NewValue1(ParaOptimization.Para.Max\[3\])	1	
set NewValue1(ParaOptimization.Para.Min\[3\])	0
set NewValue1(ParaOptimization.Para.Cluster\[3\]) 1
#
set NewValue1(ParaOptimization.Para.Name\[4\])	para5*
set NewValue1(ParaOptimization.Para.Max\[4\])	1	
set NewValue1(ParaOptimization.Para.Min\[4\])	0
set NewValue1(ParaOptimization.Para.Cluster\[4\]) 2
#
################## Genetic ######################################
set NewValue1(Genetic.NbRun)		2
set NewValue1(Genetic.NbGeneration)	2
set NewValue1(Genetic.NbElement)	10
set NewValue1(Genetic.CrossOverProbability)	0.6
set NewValue1(Genetic.MutationProbability)	0.2
set NewValue1(Genetic.SelectionType)	1    
set NewValue1(Genetic.ScalingType)	1
set NewValue1(Genetic.Sharing)		0.9
set NewValue1(Genetic.Clustering)	1
set NewValue1(Genetic.Elitism)		0.9
set NewValue1(Genetic.SA.Mutation)	0
set NewValue1(Genetic.SA.CrossOver)	0
set NewValue1(Genetic.Seed)		4
#
#
######### LINESEARCH ###########################################
set NewValue1(LineSearch.NbMaxCycle)	1
set NewValue1(LineSearch.NbMinCycle)	1
set NewValue1(LineSearch.NbMaxStab)	2
set NewValue1(LineSearch.DivInter)	10
set NewValue1(LineSearch.Alpha)	0.6
set NewValue1(LineSearch.EvolutionMini) 0.001
set NewValue1(LineSearch.Seed)		1
#
set NewValue1(LineSearch.Para.Step\[0\])	0.01
set NewValue1(LineSearch.Para.Init\[0\])	0.5
set NewValue1(LineSearch.Para.MaxInit\[0\]) 	1
set NewValue1(LineSearch.Para.MinInit\[0\]) 	0
#
set NewValue1(LineSearch.Para.Step\[1\])	0.01
set NewValue1(LineSearch.Para.Init\[1\])	0.5
set NewValue1(LineSearch.Para.MaxInit\[1\])	1
set NewValue1(LineSearch.Para.MinInit\[1\])	0
#
set NewValue1(LineSearch.Para.Step\[2\])	0.01
set NewValue1(LineSearch.Para.Init\[2\])	0.5
set NewValue1(LineSearch.Para.MaxInit\[2\])	1
set NewValue1(LineSearch.Para.MinInit\[2\])	0
#
set NewValue1(LineSearch.Para.Step\[3\])	0.01
set NewValue1(LineSearch.Para.Init\[3\])	0.5
set NewValue1(LineSearch.Para.MaxInit\[3\])	1
set NewValue1(LineSearch.Para.MinInit\[3\])	0
#
set NewValue1(LineSearch.Para.Step\[4\])	0.01
set NewValue1(LineSearch.Para.Init\[4\])	0.5
set NewValue1(LineSearch.Para.MaxInit\[4\])	1
set NewValue1(LineSearch.Para.MinInit\[4\])	0

ModifyParaValue $FileName  NewValue1
}




#############################################################################
# Procedure   : BackQuoteLine
# Description : Add a '\\' before '+' '(' ')' '|' '[' ']'
#               Argument : OldLine = string to consider
# Evaluation  : modified string 
#############################################################################
proc BackQuoteLine {OldLine} {
    set l $OldLine
    foreach sign { + ( ) | [ ]} {
	set NewLine ""
	set sign_pos [string first $sign $l]
	while { $sign_pos != -1} {
	    set NewLine "$NewLine[string range $l 0 [expr $sign_pos - 1]]\\"
	    set NewLine "$NewLine$sign"
	    set l [string range $l [expr $sign_pos + 1] [string length $l]]
	    set sign_pos [string first $sign $l]
	}
	set NewLine "$NewLine$l"
	set l $NewLine
    }
    return $NewLine
}

#############################################################################
# Procedure   : RemoveFirstLines
# Description : Remove the 2 first lines in the file given in argument
#############################################################################
proc RemoveFirstLines {file_name} {
   exec cp $file_name RemoveFirstLines.tmp
   exec tail -n +5 RemoveFirstLines.tmp > $file_name
   exec rm RemoveFirstLines.tmp
}


#############################################################################
# Procedure   : GetSeqLength
# Description : returns the number of caracters of the file given in argument
#############################################################################
proc GetSeqLength {file_name} {
    set l [exec wc -c $file_name]
    set l [string trim $l]
    set l [lindex [split $l] 0]
    return $l
}

