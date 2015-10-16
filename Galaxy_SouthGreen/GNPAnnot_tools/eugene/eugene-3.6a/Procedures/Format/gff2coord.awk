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
# $Id: gff2coord.awk,v 1.2 2004-09-16 11:56:08 cros Exp $
# ------------------------------------------------------------------
# File:     gff2coord.awk
# Contents: convert an annotation file from gff to "coord" format (.gff to .coord)
# see http://www.sanger.ac.uk/Software/formats/GFF/  for gff infos
# and evalpred.pl for coord format infos
# WARNING: 
#    gff file must contain only complete genes
#               can contain multiple genes
#    only exons marked as "First","Internal","Terminal","Single","Exon" or "exon" 
#    are reported in coord file.
# launch with:
# awk -f gff2coord.awk annotationfile.gff
# ------------------------------------------------------------------


{if(($3=="First")||($3=="Internal")||($3=="Terminal")||($3=="Single")||($3=="exon")||($3=="Exon")){beg=$4;end=$5; if($7=="-"){beg=0-beg;end=0-end}; if(ID!=$1){if(NR!=1){print ""};printf $1};printf " "beg" "end;ID=$1}}END{print"";print""}
