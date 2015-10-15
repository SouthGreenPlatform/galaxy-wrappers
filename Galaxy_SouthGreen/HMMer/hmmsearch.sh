#!/bin/bash
echo $* > $HOME/galaxy/tools/SouthGreen/HMMer/trace.txt
$HOME/galaxy/tools/SouthGreen/HMMer/hmmsearch $*;
