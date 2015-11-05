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

use LibAnnotation::Ortho::PSlib;
use LibAnnotation::Ortho::GDlib;
use FileHandle;
use GD;

use Cwd 'abs_path';
# my $dir = getcwd();

use File::Basename;
my $path = dirname(abs_path($0))."/LibAnnotation/Ortho/rgb.txt";


use constant RGB_FILE_PATH =>  dirname(abs_path($0)).'/LibAnnotation/Ortho/rgb.txt';

sub drawCmd {
    my($cmd, $xOffset, $yOffset, $height, $width, $resX, $resY, $format, $bgcolor, $fgcolor, $fhout) = @_;

    if($format eq "ps" || $format eq "pdf") {
	print $fhout printHeader(RGB_FILE_PATH);
	foreach(@$cmd) {
	    if($_->{'type'} eq "cercle") {
		my($x, $y, $r, $c) = @{$_->{'arg'}};
		$x += $xOffset;
		$y += $yOffset;
		print $fhout cercle($x, $y, $r, 0, 360, $c, 1);

	    } elsif($_->{'type'} eq "box") {
		my($x1, $y1, $x2, $y2, $color, $f) = @{$_->{'arg'}};
		if(!defined $f) { $f = 0; }
		$x1 += $xOffset;
		$y1 += $yOffset;
		print $fhout box($x1, $y1, $x2, $y2, $color, $f);

	    } elsif($_->{'type'} eq "line") {
		my($x1, $y1, $x2, $y2) = @{$_->{'arg'}};
		$x1 += $xOffset;
		$y1 += $yOffset;
		$x2 += $xOffset;
		$y2 += $yOffset;
		print $fhout rline($x1, $y1, $x2-$x1, $y2-$y1);

	    } elsif($_->{'type'} eq "text") {
		my($f, $x, $y, $text) = @{$_->{'arg'}};
		$x += $xOffset - length($text);
		$y += $yOffset;
		print $fhout text($x, $y, $text);

	    } elsif($_->{'type'} eq "textUp") {
		my($f, $x, $y, $text) = @{$_->{'arg'}};
		$x += $xOffset;
		$y += $yOffset - length($text);
		print $fhout text($x, $y, $text, 90);

	    } else {
		die "Unknown command : ", $_->{'type'}, "\n";
	    }
	}
	print $fhout printFormer();

    } elsif($format eq "jpg" || $format eq "gif" || $format eq "png") {
	my $img = new GD::Image($resX, $resY);
	my $colorMap = loadRGBFile($img);
	my $bg_color = getColor($colorMap, $bgcolor, $img);
	my $fg_color = getColor($colorMap, $fgcolor, $img);
	my $black = $img->colorAllocate(0, 0, 0);
	my $color = $black;
	if($bgcolor eq "transparent") { $img->transparent($bg_color); }
	foreach(@$cmd) {
	    if($_->{'type'} eq "cercle") {
		my($x, $y, $r, $c, $f) = @{$_->{'arg'}};
		if(!defined $f) { $f = 1; }
		$x += $xOffset;
		$y = $height - $y + $yOffset;
		$color = getColor($colorMap, $c, $img, $fg_color);
		if($f) { $img->filledArc($x, $y, $r, $r, 0, 360, $color); }
		else { $img->arc($x, $y, $r, $r, 0, 360, $color); }

	    } elsif($_->{'type'} eq "box") {
		my($x1, $y1, $x2, $y2, $c, $f) = @{$_->{'arg'}};
		if(!defined $f) { $f = 0; }
		$x1 += $xOffset;
		$y1 = $height - $y1 + $yOffset;
		$x2 += $xOffset;
		$y2 = $height - $y2 + $yOffset;
		$color = getColor($colorMap, $c, $img, $fg_color);
		if($f == 1) { $img->filledRectangle($x1, $y1, $x2, $y2, $color); }
		elsif($f == 2) { $img->rectangle($x1, $y1, $x2, $y2, $color); $img->fill(($x1+$x2)/2, ($y1+$y2)/2, $color); }
		else { $img->rectangle($x1, $y1, $x2, $y2, $color); }

	    } elsif($_->{'type'} eq "line") {
		my($x1, $y1, $x2, $y2, $c) = @{$_->{'arg'}};
		$x1 += $xOffset;
		$y1 = $height - $y1 + $yOffset;
		$x2 += $xOffset;
		$y2 = $height - $y2 + $yOffset;
		$color = getColor($colorMap, $c, $img, $fg_color);
		$img->line($x1, $y1, $x2, $y2, $color);

	    } elsif($_->{'type'} eq "text") {
		my($f, $x, $y, $text) = @{$_->{'arg'}};
		my $font = font2GD($f);
		$x += $xOffset - length($text)*($font->width/2);
		$y = $height - $y + $yOffset - ($font->height/2);
		$img->string($font, $x, $y, $text, $fg_color);

	    } elsif($_->{'type'} eq "textUp") {
		my($f, $x, $y, $text) = @{$_->{'arg'}};
		my $font = font2GD($f);
		$x += $xOffset - ($font->height/2);
		$y = $height - $y + $yOffset + length($text)*($font->width/2);
		$img->stringUp($font, $x, $y, $text, $fg_color);

	    } else {
		die "Unknown command : ", $_->{'type'}, "\n";
	    }
	}
	binmode $fhout;
	if($format eq "jpg") { print $fhout $img->jpeg(100); }
	if($format eq "gif") { print $fhout $img->gif(); }
	if($format eq "png") { print $fhout $img->png(); }
    }
}

sub loadRGBFile {
    my $img = shift;
    my $fh = new FileHandle(RGB_FILE_PATH, "r"); 
    my %colorMap = ();
    my $cpt = 0;
    if (defined $fh) {
	while(my $l = $fh->getline()) {
	    if($l =~ /^(\S+)\s+(\d+)\s+(\d+)\s+(\d+).*$/) {
		$colorMap{$1."firstComponent"} = $2;
		$colorMap{$1."secondComponent"} = $3;
		$colorMap{$1."thirdComponent"} = $4;
	    }
	}
    }
    return \%colorMap;
}

1;
