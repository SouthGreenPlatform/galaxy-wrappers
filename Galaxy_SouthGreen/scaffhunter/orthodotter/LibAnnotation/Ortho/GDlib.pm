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

use GD;

sub getNewImage {
    my($width, $height) = @_;
    return new GD::Image($width, $height);
}

sub font2GD {
    my $font = shift;
    if($font == 0) { return GD::Font->Tiny; }
    if($font == 1) { return GD::Font->Small; }
    if($font == 2) { return GD::Font->MediumBold; }
    if($font == 3) { return GD::Font->Large; }
    if($font == 4) { return GD::Font->Giant; }
    return GD::Font->Small;
}

sub getColor {
    my($colorMap, $c, $img, $default) = @_;
    if(!defined $c) {
	if(defined $default) { return $default; }
	return $img->colorAllocate(0, 0, 0);
    }
    if(defined $colorMap->{$c}) {
	return $colorMap->{$c}; 
    }
    my($f, $s, $t) = (-1, -1, -1);
    if(defined $colorMap->{$c."firstComponent"}) {
	$f = $colorMap->{$c."firstComponent"}
    }
    if(defined $colorMap->{$c."secondComponent"}) {
	$s = $colorMap->{$c."secondComponent"}
    }
    if(defined $colorMap->{$c."thirdComponent"}) {
	$t = $colorMap->{$c."thirdComponent"}
    }
    my $correctCode = 1;
    if($f < 0 || $f > 255) { $correctCode = 0; }
    if($s < 0 || $s > 255) { $correctCode = 0; }
    if($t < 0 || $t > 255) { $correctCode = 0; }
    if($correctCode) { 
	my $ret = $img->colorAllocate($f, $s, $t);
	if($ret != -1) { 
	    $colorMap->{$c} = $ret;
	    return $ret; 
	}
    }
    if(defined $default) { return $default; }
    return $img->colorAllocate(0, 0, 0);
}

1;
