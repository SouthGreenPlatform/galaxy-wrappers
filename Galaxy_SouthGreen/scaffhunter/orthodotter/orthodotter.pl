#!/usr/bin/perl -w

################################################################################
# * 
# * Copyright Jean-Marc Aury / Institut de Genomique / DSV / CEA  
# *                            <jmaury@genoscope.cns.fr>
# *
# * 
# * This software, called orthodotter is a computer program whose purpose is to 
# * draw dotplot and create dots clusters.
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

use FileHandle;
use Getopt::Long;
use File::Temp "tempfile";
use File::Basename;
use FindBin;                 
use lib "$FindBin::Bin"; 

use LibAnnotation::Ortho::OrthoDotter;
use LibAnnotation::Ortho::GeneCollection;

my $A_FILE     = "-";
my $A_OUT      = ">-";
my $A_GEOM     = 600x600;
my $A_RESX     = 600;
my $A_RESY     = 600;
my $A_RADIUS   = 2;
my $A_FORMAT   = "png";
my $A_FGCOLOR  = "black";
my $A_BGCOLOR  = "transparent";
my $A_FONTSIZE = 1;
my $A_KEEPONLYGOOD = 0;
my @A_TOPLOT = ();
my $A_VERBOSE = 0;
my $A_HELP = 0;

my $result = &GetOptions("f=s"       => \$A_FILE,
			 "o=s"       => \$A_OUT,
                         "x=i"       => \$A_RESX,
			 "y=i"       => \$A_RESY,
                         "r=f"       => \$A_RADIUS,
                         "format=s"  => \$A_FORMAT,
                         "fg=s"      => \$A_FGCOLOR,
                         "bg=s"      => \$A_BGCOLOR,
                         "fsize=i"   => \$A_FONTSIZE,
			 "filter"    => \$A_KEEPONLYGOOD,
			 "toPlot=s@" => \@A_TOPLOT,
			 "v"         => \$A_VERBOSE,
                         "h"         => \$A_HELP,
                        );

if (!$result || $A_HELP) {
    usage();
}

my %arg = ( Kfilter => 0, resX => $A_RESX, resY => $A_RESY, width => $A_RESX-100, 
	    height => $A_RESY-100, flip => 0, format => $A_FORMAT,
	    xOffset => 75, yOffset => 25, radius => $A_RADIUS, 
	    keepGoodChrom => $A_KEEPONLYGOOD, fgcolor => $A_FGCOLOR, 
	    bgcolor => $A_BGCOLOR, fontSize => $A_FONTSIZE);

if($A_FORMAT eq "ps" || $A_FORMAT eq "pdf") {
    $A_RESX = 600;
    $A_RESY = 600;
    $arg{'yOffset'} = 50;
    $arg{'xOffset'} = 50;
    $arg{'resX'}    = 600;
    $arg{'resY'}    = 600;
    $arg{'height'}  = 550;
    $arg{'width'}   = 550;
}

my $fh_in = ($A_FILE eq "-") ? new FileHandle($A_FILE) : new FileHandle($A_FILE, "r");
if (!defined $fh_in) { usage("[Error] Cannot open $A_FILE in [r] mode."); }

my $fh_out = ($A_OUT eq ">-") ? new FileHandle($A_OUT) : new FileHandle($A_OUT, "w");
if (!defined $fh_out) { usage("[Error] Cannot open $A_OUT in [w] mode."); }

if (scalar(@A_TOPLOT) < 1) { usage("[Error] No dataset to plot. Use -toPlot argument"); }

warn "Load collections ...\n" if $A_VERBOSE;
my ($collections, $nb) = load_collections_from_file($fh_in);
warn "$nb collections loaded\n" if $A_VERBOSE;

my @dotters = ();
foreach my $opt (@A_TOPLOT) {
    if($opt =~ /^(\w+)=?([\w|\,]*):(\w+)=?([\w|\,]*):?(\w*):?(\w*)=?(\d*),?(\d*)$/) {
	if(!defined $collections->{$1}) { usage("[Error] Could not find dataset : $1"); }
	if(!defined $collections->{$3}) { usage("[Error] Could not find dataset : $3"); }
	my $color = (defined $5) ? $5 : "black";
	my $slc = (defined $6) ? $6 : 0;
	my $max_dist = (defined $7) ? $7 : 60;
	my $nb_genes = (defined $8) ? $8 : 10;
	my $dotter = new LibAnnotation::Ortho::OrthoDotter($A_RESX, $A_RESY, $A_RADIUS, $A_FORMAT, $A_BGCOLOR, 
							   $A_FGCOLOR, $color, $A_FONTSIZE, $A_KEEPONLYGOOD);
	my (%filter_x, %filter_y);
	if(defined $2 && $2 ne "") { 
	    print "filter : $2\n";
	    my @filterX = split(",", $2);
	    foreach(@filterX) { $filter_x{$_}=1; }
	    $dotter->filterX(\%filter_x);
	}
	if(defined $4 && $4 ne "") {
	    print "filter : $4\n";
	    my @filterY = split(",", $4);
	    foreach(@filterY) { $filter_y{$_}=1; }
	    $dotter->filterY(\%filter_y);
	}
	warn "plot $1 against $3\n" if $A_VERBOSE;
	$dotter->collectionX($collections->{$1});
	$dotter->collectionY($collections->{$3});
	if($slc) { $dotter->_doCluster($max_dist, $nb_genes); }
	push(@dotters, $dotter);
    } else {
	usage("[Error] bad syntax in argument -toPlot : $opt");
    }
}

my $cpt = 0;
my @instructions = ();
foreach(@dotters) {
    $cpt++;
    if($cpt != 1) { $_->draw_limits(0); }
    if($cpt != scalar(@dotters)) { $_->draw(\@instructions, 0, $fh_out); }
    else { $_->draw(\@instructions, 1, $fh_out); }
}
	
sub usage {
    my $msg = shift;
    if(defined $msg) { warn "\n$msg\n"; }
    my $usage = "
--------------------------------------------------------------------------------
orthodotter - Plot orthologous genes on an oxford grid.
       -f <file>     : input file, containing orthologous genes, default is stdin
                       species chr-name start end species chr-name start end
       -toPlot <arg> : give the x and y sets and the color separated by double-dots,
                       for example set1:set2:red will plot set1 on x, set2 on y with
                       red points. Could give several -toPlot arguments.
                       To launch the clustering of dots, use extra-option 1=dist,min_nb_genes
                       where dist is the minimal distance (euclidian) between two points and min_nb_genes the minimal
                       number of genes in a cluster to be valid.
       -o <file>     : output file, default is stdout
       -x <int>      : resolution of x axis, default is 600
       -y <int>      : resolution on y axis, default is 600
       -r <int>      : radius of circle representing orthologous genes
       -format       : could be png, gif, jpg, pdf or ps. Default is png.
       -fg           : foreground color, default is black
       -bg           : background color, default is transparent
       -fSize <int>  : fontSize, default is 1
       -filter       : check chromosome names
       -h            : help
--------------------------------------------------------------------------------
orthodotter -f Vigne_Banane.ortho -toPlot Vigne:Banane:black:1=10,5 -x 1200 -y 1200 -bg white -o Vigne_vs_Banane.png > Vigne_vs_Banane.clusters
--------------------------------------------------------------------------------\n";
    die $usage;
}

sub load_collections_from_file {
    my $fh_in = shift;
    my %collection = ();
    my $nb = 0;
    while(my $line = $fh_in->getline()) {
	if($line =~ /^(\S+)\s+(\S+)\s+(\d+|\.)\s+(\d+|\.)\s+(\S+)\s+(\S+)\s+(\d+|\.)\s+(\d+|\.)\s*(\d*|\.)\s*(\S*)$/) {
	    my ($g1, $g2);
	    if($2 ne ".") {
		if(!defined $collection{$1}) { 
		    my $c = new LibAnnotation::Ortho::GeneCollection(); 
		    $c->name($1);
		    $collection{$1} = $c;
		    $nb++;
		}
		$g1 = ($collection{$1})->addGene($2, $3, $4);
	    }
	    if($6 ne ".") {
		if(!defined $collection{$5}) { 
		    my $c = new LibAnnotation::Ortho::GeneCollection(); 
		    $c->name($5);
		    $collection{$5} = $c;
		    $nb++;
		}
		$g2 = ($collection{$5})->addGene($6, $7, $8);
		if($2 ne ".") {
		    my $score = (defined $9 && $9 ne ".") ? $9 : 0;
		    my $color = (defined $10) ? $10 : "none";
		    my $hit = new LibAnnotation::Ortho::GeneHit($g1, $g2, $score, $color);
		    $g1->addmatch($hit);
		    $g2->addmatch($hit);
		}
	    }
	}
    }
    return (\%collection, $nb);
}
