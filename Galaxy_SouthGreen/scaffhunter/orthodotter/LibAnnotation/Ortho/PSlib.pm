#!/usr/local/bin/perl -w

################################################################################
# * 
# * Copyright Jean-Marc Aury / Institut de Genomique / DSV / CEA  
# *                            <jmaury@genoscope.cns.fr>
# *
# * 
# * This software is governed by the CeCILL license under French law and
# * abiding by the rules of distribution of free software.  You can  use, 
# * modify and/ or redistribute the software under the terms of the CeCILL
# * license as circulated by CEA, CNRS and INRIA at the following URL
# * "http://www.cecill.info". 
# * 
# * As a counterpart to the access to the source code and  rights to copy,
# * modify and redistribute granted by the license, users are provided only
# * with a limited warranty  and the software's author,  the holder of the
# * economic rights,  and the successive licensors  have only  limited
# * liability. 
# * 
# * In this respect, the user's attention is drawn to the risks associated
# * with loading,  using,  modifying and/or developing or reproducing the
# * software by the user in light of its specific status of free software,
# * that may mean  that it is complicated to manipulate,  and  that  also
# * therefore means  that it is reserved for developers  and  experienced
# * professionals having in-depth computer knowledge. Users are therefore
# * encouraged to load and test the software's suitability as regards their
# * requirements in conditions enabling the security of their systems and/or 
# * data to be ensured and,  more generally, to use and operate it in the 
# * same conditions as regards security. 
# * 
# * The fact that you are presently reading this means that you have had
# * knowledge of the CeCILL license and that you accept its terms.
################################################################################

use strict;

sub printHeader {
    my $rgbFile = shift;
    my $res = "%!PS-Adobe-3.0\n";
    $res .= "%%DocumentData: Clean8bit\n";
    $res .= "%%PageOrder: Ascend\n";
    $res .= "%%Pages: 1\n";
    $res .= "%%DocumentFonts: Helvetica\n";
    $res .= "%%EndComments\n";
    
    $res .=  "/cm {28.3464567 mul} def\n";
    $res .=  "%%Bounding-Box: 0 cm 0 cm 21 cm 29.7 cm\n";
    
    $res .=  "\n";
    $res .=  "/Arial findfont\n";
    $res .=  "10 scalefont\n";
    $res .=  "setfont\n";

    if(defined $rgbFile) { $res .= loadRGBColor($rgbFile); }
    else { $res .= getColorDef(); }
	
    return $res;
}

sub loadRGBColor {
    my $file = shift;
    my $fh = new FileHandle($file, "r"); 
    my $colorMap = "";
    my $cpt = 0;
    if (defined $fh) {
	while(my $l = $fh->getline()) {
	    if($l =~ /^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d\.\d+)\s+(\d\.\d+)\s+(\d\.\d+)\s*$/) {
		$colorMap .= "/$1 [ $5 $6 $7 ] def\n";
	    }
	}
    }
    if($colorMap eq "") { warn "colorMap is empty\n"; }
    $colorMap .= "/color{\n";
    $colorMap .= "aload pop setrgbcolor\n";
    $colorMap .= "}def\n";
    $colorMap .= "\n";
    return $colorMap;
}

sub printFormer {
    my $res = "showpage\n";
    return $res;
}

sub getColorDef {
    my $res = "/grey62 [ 0.62 0.62 0.62 ] def\n";
    $res .= "/black [ 0.00 0.00 0.00 ] def\n";
    $res .= "\n";
    $res .= "/color{\n";
    $res .= "aload pop setrgbcolor\n";
    $res .= "}def\n";
    $res .= "\n";
}

sub text {
    my($X, $Y, $text, $angle) = @_;
    $angle = $angle || 0;
    my $res .= "newpath\n";
    $res .= "$X $Y moveto\n";
    if($angle != 0) {
	$res .= "$angle rotate\n";
    }
    $res .= "($text) show\n";
    if($angle != 0) {
	$res .= "-$angle rotate\n";
    }
    $res .= "closepath\n";
    return $res;
}

sub rline {
    my($X1, $Y1, $X2, $Y2) = @_;
    my $res = "newpath\n";
    $res .= "$X1 $Y1 moveto\n";
    $res .= "$X2 $Y2 rlineto\n";
    $res .= "closepath\n";
    $res .= "stroke\n";
    return $res;
}

sub cercle {
    my($X, $Y, $R, $A, $B, $color, $fill) = @_;
    my $res = "newpath\n";
    if($fill) {
	$res .= "0 setlinewidth\n";
	$res .= "$X $Y $R $A $B arc closepath\n";
	$res .= "gsave\n";
	$res .= "$color color fill\n";
	$res .= "grestore\n";
    } else {
	$res .= "$X $Y $R $A $B arc closepath\n";
    }
    $res .= "stroke\n";
    $res .= "grestore\n";
    return $res;
}

sub box {
    my($X, $Y, $L, $H, $color, $fill) = @_;
    my $res = "newpath\n";
    $res .= "$X $Y moveto\n";
    $res .= "$L 0 rlineto\n";
    $res .= "0 $H rlineto\n";
    $res .= "-$L 0 rlineto\n";
    $res .= "closepath\n";
    if($fill) { 
	$res .= "gsave\n";
	$res .= "$color color fill\n";
	$res .= "grestore\n";
    }
    $res .= "stroke\n";
    return $res;
}

1;
