#!/bin/bash
directory=`dirname $0`

module load system/java/jdk6
java -mx4096m -jar $directory/RapGreen.jar $*;
