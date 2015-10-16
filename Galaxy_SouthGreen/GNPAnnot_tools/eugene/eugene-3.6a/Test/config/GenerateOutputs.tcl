#!/usr/bin/tclsh

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
# $Id: GenerateOutputs.tcl,v 1.21 2010-01-25 16:34:16 sallet Exp $
# ------------------------------------------------------------------
# File:     GenerateOutput.tcl
# Contents: Generation of the reference files for the test suite
# Save the reference output (both stdout and stderr) for each test
# ------------------------------------------------------------------



############################# Environment variables ############################
set action "Generate"

source ./TestProc.tcl
source ./TestVar.tcl
################################################################################


################################################################################
proc AskReplace { file_name } {
	puts "There is a difference between the old reference file $file_name and the new reference file generated.
Replace the old one ? (enter 'Y' for yes else another key for no)" 
}
################################################################################


# Erase all old reference files
puts "Remove all reference files ? (enter 'Y' if you are OK else another key)"
set key [gets stdin]
if {$key=="Y"} {
    foreach f [exec ls ${OUTPUT_DIR}] {
	if {![string match CVS $f]} {
	    exec rm ${OUTPUT_DIR}/$f
	}
    }
    set erase 1
} else {
    set erase 0
}


########################## Init the parameter file ########################
# Convention: if a test requires different parameters values,
#             the initial values are restored after the test
#             BEWARE does not concern sensor use parameters
###########################################################################
# Copy locally the default parameter file
exec cp  ${EUGENE_DIR}/cfg/eugene.par $EUGENE_TEST_PAR
# Init parameters values
InitParameterFile $EUGENE_TEST_PAR $AllSensorsList $EUGENE_DIR

########################################################################
##################        Units tests       ############################
########################################################################
foreach sensor $AllSensorsList {

   # Get stderr and stdout
    if {$sensor == "Est"} {
	catch { eval exec $EUGENE_DIR/$EUGENE -A $EUGENE_TEST_PAR \
		    $OPTIONS(Sensor) -D Sensor.${sensor}.use=2 -D Sensor.NG2.use=1 \
		    $SEQ_DIR/$SEQ(Sensor).tfa 2> tmp%std } 
    } else {
	if {$sensor == "Homology"} {
	    catch { eval exec $EUGENE_DIR/$EUGENE -A $EUGENE_TEST_PAR \
			$OPTIONS(Sensor) -D Sensor.${sensor}.use=2 \
			$SEQ_DIR/$SEQ(Sensor).tfa 2> tmp%std }
	} else {
	    if {$sensor == "Tester"} {
		catch { eval exec $EUGENE_DIR/$EUGENE -A $EUGENE_TEST_PAR \
			    $OPTIONS(Sensor) -D Sensor.${sensor}.use=1 \
			    $SEQ_DIR/exSeqHom.fasta 2> tmp%stderr > tmp%stdout }
		catch {exec cat tmp%stderr tmp%stdout > tmp%std} 
		catch {exec rm  tmp%stderr tmp%stdout} 
	    } else {
		catch { eval exec $EUGENE_DIR/$EUGENE -A $EUGENE_TEST_PAR \
			    $OPTIONS(Sensor) -D Sensor.${sensor}.use=1 \
			    $SEQ_DIR/$SEQ(Sensor).tfa 2> tmp%std }
	    }
	}
    }
    
    # Open files
    if {$sensor == "Tester"} {
	set out  [open exSeqHom.egn.debug {RDONLY}]
    } else {
	set out [open $SEQ(Sensor).egn.debug {RDONLY}]
    }
    set err [open tmp%std {RDONLY}]
    set std [open tmp%GenerateOutputs w+]

    # Copy stderr and stdout in the reference file
    set f [read -nonewline $err]
    puts $std $f
    set flux [read -nonewline $out]
    puts $std $flux

    # Close files
    close $out
    close $err
    close $std

    # Remove the 2 first lines related to version number
    RemoveFirstLines tmp%GenerateOutputs
    
    if {$erase == 1 || ![file exists $OUTPUT_DIR/Output_${sensor}]} {
	exec cp tmp%GenerateOutputs $OUTPUT_DIR/Output_${sensor}
    } elseif {[catch {exec diff $OUTPUT_DIR/Output_${sensor} tmp%GenerateOutputs}]} {
	AskReplace $OUTPUT_DIR/Output_${sensor} 
	if {[gets stdin] == "Y"} {
	    exec cp tmp%GenerateOutputs $OUTPUT_DIR/Output_${sensor}
	} else {
	    exec cp tmp%GenerateOutputs $OUTPUT_DIR/Output_${sensor}.new
	}
    }

    # Remove all temporary files
    catch {exec rm tmp%GenerateOutputs tmp%stderr tmp%std}
    catch {exec rm $SEQ(Sensor).egn.debug $SEQ(Sensor).misc_info exSeqHom.egn.debug exSeqHom.misc_info}
    # Remove created files, note that eugeneTest.par.<date>.OPTI remains
    catch {exec rm Sensor.EuStop.SpSn}

    puts "Reference files for $sensor unit test created or checked."
}




########################################################################
########################      FUNCTIONAL TESTS      ####################
########################################################################
foreach TEST $FunctionalTestList {

    # Preparation of the parameter file with the correct sensors
    foreach sensor $SensorsList($TEST) \
	{set NewValue${TEST}(Sensor.${sensor}.use) 1}
    ModifyParaValue $EUGENE_TEST_PAR  NewValue${TEST}

    # Get the sequence length to have only one png file
    set l [GetSeqLength $SEQ_DIR/$SEQ($TEST)]
    
    # Save output of software and treat them
    catch {eval exec $EUGENE_DIR/$EUGENE -A $EUGENE_TEST_PAR $OPTIONS($TEST) -l $l $SEQ_DIR/$SEQ($TEST) > tmp%FunctionalTest}

    # 1/ image file
    if {$erase == 1 || ![file exists $OUTPUT_DIR/Output_${TEST}.png]} {
	exec cp $IMG($TEST) $OUTPUT_DIR/Output_${TEST}.png
    } elseif {[catch {exec diff $IMG($TEST) $OUTPUT_DIR/Output_${TEST}.png}]} {
	AskReplace $OUTPUT_DIR/Output_${TEST}.png
	if {[gets stdin] == "Y"} {
	    eval exec cp $IMG($TEST) $OUTPUT_DIR/Output_${TEST}.png
	} else {
	    eval exec cp $IMG($TEST) $OUTPUT_DIR/Output_${TEST}.png.new
	}
    }
    # Remove temporary file
    exec rm $IMG($TEST)

    # 2/ Remove the first lines related to version number
    #RemoveFirstLines tmp%FunctionalTest

    # 3/ reference file in sp
    if {$erase == 1 || ![file exists $OUTPUT_DIR/Output_${TEST}]} {
	exec cp tmp%FunctionalTest $OUTPUT_DIR/Output_${TEST}
    } elseif {[catch {exec diff $OUTPUT_DIR/Output_${TEST} tmp%FunctionalTest}]} {
	AskReplace $OUTPUT_DIR/Output_${TEST}
	if {[gets stdin] == "Y"} {
	    exec cp tmp%FunctionalTest $OUTPUT_DIR/Output_${TEST}
	} else {
	    exec cp tmp%FunctionalTest $OUTPUT_DIR/Output_${TEST}.new
	}
    }
    # Remove temporary files tmp%stderr tmp%stdout 
    exec rm tmp%FunctionalTest

    # Restore initial parameters values
    InitParameterFile $EUGENE_TEST_PAR $AllSensorsList $EUGENE_DIR

    puts "Reference files for $TEST functional test created or checked."
}

# Remove temporary files (misc_info)
set seq [eval exec ls .]
foreach f $seq {
    if { [string match *.misc_info $f]} {
	exec rm $f
    }
}


########################################################################
######################## Arabidopsis Sequences base Tests ##########################
########################################################################

foreach TEST $ArabidopsisTestList {

    # Preparation of the parameter file with the correct sensors
    foreach sensor $SensorsList($TEST) \
	{set NewValueBase(Sensor.${sensor}.use) 1}
    ModifyParaValue $EUGENE_TEST_PAR  NewValueBase
    

    catch {eval exec $EUGENE_DIR/$EUGENE -A $EUGENE_TEST_PAR $OPTIONS($TEST) $SEQ($TEST) > tmp%stdout}
   
    if {$erase == 1 || ![file exists $OUTPUT_DIR/$FILE_REF($TEST)]} {
	exec cp tmp%stdout $OUTPUT_DIR/$FILE_REF($TEST)
    } elseif {[catch {exec diff $OUTPUT_DIR/$FILE_REF($TEST) tmp%stdout}]} {
	AskReplace $OUTPUT_DIR/$FILE_REF($TEST)
	if {[gets stdin] == "Y"} {
	    exec cp tmp%stdout $OUTPUT_DIR/$FILE_REF($TEST)
	} else {
	    exec cp tmp%stdout $OUTPUT_DIR/$FILE_REF($TEST).new
	}
    }
    
    if {$TEST == "Araset"} {
		set TEST2 "ArasetSpSn"
		eval exec $PRG_EVAL_PRED $FILE_COORD($TEST2) tmp%stdout > tmp%stdout2
		RemoveFirstLines tmp%stdout2
		if {[catch {exec diff $OUTPUT_DIR/$FILE_REF($TEST2) tmp%stdout2}]} {
			AskReplace $OUTPUT_DIR/$FILE_REF($TEST2)
			if {[gets stdin] == "Y"} {
	    			exec cp tmp%stdout2 $OUTPUT_DIR/$FILE_REF($TEST2)
			} else {
	    			exec cp tmp%stdout2 $OUTPUT_DIR/$FILE_REF($TEST2).new
			}
		exec rm tmp%stdout2
	    }
    }
    # remove temporary file	
    exec rm tmp%stdout
    set seq [eval exec ls .]
    foreach f $seq {
	if { [string match *.misc_info $f] || [string match *.egn.debug $f]} {
	    exec rm $f
	}
    }
    
    # Restore initial parameters values
    InitParameterFile $EUGENE_TEST_PAR $AllSensorsList $EUGENE_DIR

    puts "Reference files for $TEST test created or checked."

}


########################################################################
##################### Parameters optimization ##########################
########################################################################

catch {eval exec $EUGENE_DIR/$EUGENE test -A $EUGENE_TEST_PAR -D ParaOptimization.Use=1 > tmp%stdout}

if {$erase == 1 || ![file exists $OUTPUT_DIR/$FILE_REF(Optimization)]} {
    exec cp tmp%stdout $OUTPUT_DIR/$FILE_REF(Optimization)
} elseif {[catch {exec diff $OUTPUT_DIR/$FILE_REF(Optimization) tmp%stdout}]} {
    AskReplace $OUTPUT_DIR/$FILE_REF(Optimization)
    if {[gets stdin] == "Y"} {
	exec cp tmp%stdout $OUTPUT_DIR/$FILE_REF(Optimization)
    } else {
	exec cp tmp%stdout $OUTPUT_DIR/$FILE_REF(Optimization).new
    }
}
# remove temporary file	
exec rm tmp%stdout

puts "Reference files for Optimization test created or checked."


# Indicate the end of the reference files generation
puts "Reference files generated in the $OUTPUT_DIR directory."

catch {eval exec rm ./$EUGENE_TEST_PAR}

