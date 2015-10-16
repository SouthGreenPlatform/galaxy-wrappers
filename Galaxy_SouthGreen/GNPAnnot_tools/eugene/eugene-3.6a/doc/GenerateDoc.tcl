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
#  $Id: GenerateDoc.tcl,v 1.12 2009-01-20 15:35:38 tschiex Exp $
# ------------------------------------------------------------------
# File:     GenerateDoc.tcl
# Contents: Generation of the eugene documentation
# ------------------------------------------------------------------


#===========================================================================
# definition of variables related to I/O

set SEQ Sequences/SYNO_ARATH.fasta
set SEQALT Sequences/At5g18830.fasta.genomicAJ011613.fasta
set EUGENE ../src/eugene
set env(EUGENEDIR) [pwd]/..

set FIC_TEX "Doc.tex"
set FIC_TEX_TMP "EuGeneDoc"
set FIC_TMP "tmp"
set CMDFLAGS_INDEX "CmdFlags"

set Cmd_end " >& $FIC_TMP"
set Flag(1) EXECUTION_TRACE1; set Cmd_begin(1) ""; set Cmd(1) "$EUGENE -s -po $SEQ"
set Flag(2) EXECUTION_TRACE2; set Cmd_begin(2) ""; set Cmd(2) "$EUGENE -s -po -d $SEQ"
set Flag(3) EXECUTION_TRACE3; set Cmd_begin(3) ""; set Cmd(3) "$EUGENE -s -po -d -E $SEQ"
set Flag(4) EXECUTION_TRACE4; set Cmd_begin(4) ""; set Cmd(4) "$EUGENE -s -po -d -b012 -B $SEQ"
set Flag(5) EXECUTION_TRACE5; set Cmd_begin(5) ""; set Cmd(5) "$EUGENE -s -a -po $SEQALT"

set nbflags 5
#===========================================================================


# read tex doc file
set f [open $FIC_TEX r]
set content [read $f]
close $f


# substitution of flags with execution trace
set new_content ""
for {set i 1} {$i<= $nbflags} {incr i} {
    set FlagPos [string first $Flag($i) $content]
    if { $FlagPos == -1 } {
	puts "ERROR no $Flag($i) set."
	set i $nbflags
    } else {
	set begin_new_content [string range $content 0 [expr $FlagPos - 1]]
	set new_content "$new_content $begin_new_content" 
	set begin_content [expr $FlagPos + [string length $Flag($i)]]
	set content [string range $content $begin_content [string length $content]]

	if {[string length $Cmd_begin($i)] != 0} {
	    eval exec $Cmd_begin($i)
	}

	eval exec $Cmd($i) $Cmd_end

	# write the executed command in the documentation before the result
	set new_content "$new_content>$Cmd($i)\n"
	set f [open $FIC_TMP r]
	set new_content "$new_content[read $f]"
	close $f

	# for (3) & (4) print the misc_info file
	if { $i==3 || $i==4 } {
	    set new_content "$new_content> cat SYNO_ARATH.misc_info\n"
	    set f [open SYNO_ARATH.misc_info r]
	    set new_content "$new_content[read $f]"
	    close $f
	}
    }	
}
set new_content "$new_content $content"

# write new tex doc file
set f [open $FIC_TEX_TMP.tex w]
puts $f $new_content
close $f

# clean directory (beware except image .png)
exec rm  $FIC_TMP
exec rm  SYNO_ARATH.misc_info
exec rm At5g18830.fasta.genomicAJ011613.misc_info

# ask for compilation
exec pdflatex -interaction=nonstopmode $FIC_TEX_TMP.tex
catch {exec makeindex $CMDFLAGS_INDEX}
exec pdflatex -interaction=nonstopmode $FIC_TEX_TMP.tex
puts "eugene documentation generated."
