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

#!/usr/bin/perl -w

=head1 NAME

OrthoDotter - Class to produce oxford grid of complete genome

=head1 AUTHORS

2005, Genoscope - CNS, Jean-Marc Aury, jmaury@genoscope.cns.fr

=cut

package LibAnnotation::Ortho::OrthoDotter;

use FindBin;                 
use strict;
use FileHandle;
use LibAnnotation::Ortho::GeneCollectionIterator;
use LibAnnotation::Ortho::drawLib;

use Graph;
use Graph::Undirected;

use Math::Cephes qw(sqrt pow);

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {};
    bless $self, $class;
    $self->resX(shift || 1000);
    $self->resY(shift || 1000);
    $self->gene_radius(shift || 2);
    $self->format(shift || "png");
    $self->bgcolor(shift || "white");
    $self->fgcolor(shift || "black");
    $self->gene_color(shift || "black");
    $self->fontsize(shift || 1);
    $self->filterseq(shift || 0);
    $self->collectionX();
    $self->collectionY();
    $self->legendX("");
    $self->legendY("");
    $self->filterX();
    $self->filterY();
    $self->filterOnX(0);
    $self->filterOnY(0);
    $self->draw_genes(1);
    $self->draw_limits(1);
    return $self;
}

sub resX { 
    my($self, $arg) = @_; 
    if(defined $arg) { $self->{'resX'} = $arg; $self->{'width'} = $arg - 100; } 
    else { return $self->{'resX'}; }
}
sub resY { 
    my($self, $arg) = @_; 
    if(defined $arg) { $self->{'resY'} = $arg; $self->{'height'} = $arg - 100; } 
    else { return $self->{'resY'}; }
}
sub collectionX { 
    my($self, $arg) = @_; 
    if(defined $arg) { 
	my $c = $arg;
	if($self->filterOnX()) {
	    $c = $arg->dup($self->filterX());
	}
	$self->{'collectionX'} = $c; 
    } 
    else { return $self->{'collectionX'}; }
}
sub collectionY { 
    my($self, $arg) = @_; 
    if(defined $arg) { 
	my $c = $arg;
	if($self->filterOnY()) {
	    $c = $arg->dup($self->filterY());
	}
	$self->{'collectionY'} = $arg; 
    } 
    else { return $self->{'collectionY'}; }
}

sub width { my($self, $arg) = @_; if(defined $arg) { $self->{'width'} = $arg; } else { return $self->{'width'}; }}
sub height { my($self, $arg) = @_; if(defined $arg) { $self->{'height'} = $arg; } else { return $self->{'height'}; }}
sub gene_radius { my($self, $arg) = @_; if(defined $arg) { $self->{'gene_radius'} = $arg; } else { return $self->{'gene_radius'}; }}
sub format { my($self, $arg) = @_; if(defined $arg) { $self->{'format'} = $arg; } else { return $self->{'format'}; }}
sub bgcolor { my($self, $arg) = @_; if(defined $arg) { $self->{'bgcolor'} = $arg; } else { return $self->{'bgcolor'}; }}
sub fgcolor { my($self, $arg) = @_; if(defined $arg) { $self->{'fgcolor'} = $arg; } else { return $self->{'fgcolor'}; }}
sub gene_color { my($self, $arg) = @_; if(defined $arg) { $self->{'gene_color'} = $arg; } else { return $self->{'gene_color'}; }}
sub fontsize { my($self, $arg) = @_; if(defined $arg) { $self->{'fontsize'} = $arg; } else { return $self->{'fontsize'}; }}
sub filterseq { my($self, $arg) = @_; if(defined $arg) { $self->{'filterseq'} = $arg; } else { return $self->{'filterseq'}; }}
sub legendX { my($self, $arg) = @_; if(defined $arg) { $self->{'legendX'} = $arg; } else { return $self->{'legendX'}; }}
sub legendY { my($self, $arg) = @_; if(defined $arg) { $self->{'legendY'} = $arg; } else { return $self->{'legendY'}; }}
sub filterX { my($self, $arg) = @_; if(defined $arg) { $self->{'filterX'} = $arg; $self->filterOnX(1); } else { return $self->{'filterX'}; }}
sub filterY { my($self, $arg) = @_; if(defined $arg) { $self->{'filterY'} = $arg; $self->filterOnY(1); } else { return $self->{'filterY'}; }}
sub filterOnX { my($self, $arg) = @_; if(defined $arg) { $self->{'filterOnX'} = $arg; } else { return $self->{'filterOnX'}; }}
sub filterOnY { my($self, $arg) = @_; if(defined $arg) { $self->{'filterOnY'} = $arg; } else { return $self->{'filterOnY'}; }}
sub draw_genes { my($self, $arg) = @_; if(defined $arg) { $self->{'draw_genes'} = $arg; } else { return $self->{'draw_genes'}; }}
sub draw_limits { my($self, $arg) = @_; if(defined $arg) { $self->{'draw_limits'} = $arg; } else { return $self->{'draw_limits'}; }}

# Return the grid
sub _drawGrid {
    my $self = shift;
    return {'type' => "box", 
	    'arg' => [ 0, 0, $self->width(), $self->height() ]};
}

sub _drawLegendX {
    my $self = shift;
    my $offset = -30;
    my $yOffset = 25;
    if(defined $self->format() && ($self->format() eq "ps" || $self->format() eq "pdf")) { $offset = 10; }
    return {'type' => "text", 
	    'arg' => [2, $self->width()/2, -$yOffset+$offset, $self->legendX()]};
}

sub _drawLegendY {
    my $self = shift;
    my $xOffset = 75;
    return {'type' => "textUp", 
	    'arg' => [2, -$xOffset+10, $self->height()/2, $self->legendY()]};
}

# Draw the orthologous genes
sub _drawOrtho {
    my($self, $XPOS, $YPOS) = @_;
    my @res = ();
    my ($filterX, $filterY) = ($self->filterX(), $self->filterY());
    my $it = new LibAnnotation::Ortho::GeneCollectionIterator($self->collectionX());
    while(my $g = $it->next()) {
	my $matches = $g->matches();
	foreach my $m (@$matches) {
	    my $ortho_gene = ($m->gene1() == $g) ? $m->gene2() : $m->gene1();
	    # On ne dessine que les genes qui correspondent a la collection a plotter en Y
	    if($ortho_gene->collection() == $self->collectionY()) {
		my $color = ($m->color() ne "none") ? $m->color() : $self->gene_color();
		my $x = $XPOS->{$g->seq()."@".$g->start()."@".$g->stop()};
		my $y = $YPOS->{$ortho_gene->seq()."@".$ortho_gene->start()."@".$ortho_gene->stop()};
		push(@res, {'type' => "cercle", 
			    'arg' => [ $x, $y, $self->gene_radius(), $color ]});
	    }
	}
    }
    return \@res;
}

sub _drawLimits {
    my ($self, $chromLimit, $axe) = @_;

    my @res;
    my ($prec, $posText) = (0, -10);
    foreach my $key (sort { $chromLimit->{$a} <=> $chromLimit->{$b} } (keys(%$chromLimit))) {
		my ($seq, $val) = ($key, $chromLimit->{$key});
		my ($st, $en) = ($val->{'start'}, $val->{'end'});
		my($xl1, $yl1) = ($st, 0);
		my($xl2, $yl2) = ($st, $self->height());
		if (!$axe) { $yl2 = $self->width(); ($xl1, $yl1) = ($yl1, $xl1); ($xl2, $yl2) = ($yl2, $xl2);}
		push(@res, {'type' => "line", 'arg' => [$xl1, $yl1, $xl2, $yl2]});

		my $pos = $st + (($en - $st + 1) / 2);
		my($xt, $yt) = ($pos, $posText);
		if (!$axe) { ($xt, $yt) = ($yt, $xt); }
		# Remove 4 characters of scaffold name
		my $label = substr($seq,4);  
		push(@res, {'type' => "text", 'arg' => [$self->fontsize(), $xt, $yt, $label]});

		if($posText == -24) { $posText = -10; } else { $posText = -24; }
		$prec = $en;
    }
    return \@res;
}

sub draw {
    my ($self, $instructions, $draw, $fh_out) = @_;
    
    push(@$instructions, $self->_drawGrid);
    if($self->legendX() ne "") { push(@$instructions, $self->_drawLegendX()); }
    if($self->legendY() ne "") { push(@$instructions, $self->_drawLegendY()); }

    my($xlimit, $ylimit, $xpos, $ypos);
    if($self->draw_genes() || $self->draw_limits()) {
	($self->collectionX())->axeSize($self->width());
	($self->collectionY())->axeSize($self->height());
	($xlimit, $xpos) = ($self->collectionX())->projectOnAxis();
	($ylimit, $ypos) = ($self->collectionY())->projectOnAxis();
    }

    if($self->draw_genes()) { 
	my $tmp = $self->_drawOrtho($xpos, $ypos);
	push(@$instructions, @$tmp); 
    }
    if($self->draw_limits()) { 
	my $tmp = $self->_drawLimits($xlimit, 1);
	push(@$instructions, @$tmp);
	$tmp = $self->_drawLimits($ylimit, 0);
	push(@$instructions, @$tmp);
    }

    if($draw) {
	drawCmd($instructions, 50, 50, $self->height(), $self->width(), $self->resX(), 
		$self->resY(), $self->format(), $self->bgcolor(), $self->fgcolor(), , $fh_out);
    } else {
	return $instructions;
    }
}

# SL Clustering
sub _doCluster {
    my($self, $max_dist, $nb_genes) = @_;
    my $seqx = ($self->collectionX())->getSeqNames();
    my $seqy = ($self->collectionY())->getSeqNames();
    my @clusters;
    foreach my $sx (@$seqx) {
	foreach my $sy (@$seqy) {
	    my($genesx, $matchesx, $uniquex) = ($self->collectionX())->getGenesFromSeq($sx, $sy);
	    if($uniquex < $nb_genes) { next; }
	    my($genesy, $matchesy, $uniquey) = ($self->collectionY())->getGenesFromSeq($sy, $sx);
	    if($uniquey < $nb_genes) { next; }
	    if(scalar(@$matchesx) != scalar(@$matchesy)) {
		warn "Different number of matches between sequences $sx and $sy\n";
	    }
	    my $c = $self->_clusterRect($genesx, $genesy, $matchesx, $max_dist, $nb_genes);
	    if(scalar(@$c) != 0) {
		warn scalar(@$matchesx), " matches between sequences $sx and $sy : ", scalar(@$c), " clusters created\n";
		push(@clusters, @$c);
	    }
	}
    }
    my @sortC = sort { scalar(@$b) <=> scalar(@$a) } @clusters;
    my $cpt=0;
    foreach my $cluster (@sortC) {
	foreach my $p (@$cluster) {
	    print "$cpt ", ($p->{'gx'})->collection()->name(), " ", ($p->{'gx'})->seq(), " ", 
	    ($p->{'gx'})->start(), " ", ($p->{'gx'})->stop(), " ", 
	    ($p->{'gy'})->collection()->name(), " ", ($p->{'gy'})->seq(), " ", 
	    ($p->{'gy'})->start(), " ", ($p->{'gy'})->stop(), 
	    " size= ", scalar(@$cluster), "\n";
	}
	$cpt++;
    }
    
}

sub _clusterRect {
    my($self, $genesx, $genesy, $matchesx, $max_dist, $nb_genes) = @_;
    my @sortx = sort { $a->start() <=> $b->start() } @$genesx;
    my @sorty = sort { $a->start() <=> $b->start() } @$genesy;
    my $idx = 0;
    
    # On cree les coordonnees de chaque gene
    my (%xidx, %yidx);
    foreach(@sortx) { $xidx{$_} = $idx; $idx++; }
    $idx = 0;
    foreach(@sorty) { $yidx{$_} = $idx; $idx++; }

    # On cree les points contenus dans le rectangle
    my @dots;
    foreach my $m (@$matchesx) {
	my $x = (defined $xidx{$m->gene1()}) ? $m->gene1() : $m->gene2();
	my $y = (defined $yidx{$m->gene1()}) ? $m->gene1() : $m->gene2();
	push(@dots, { 'gx' => $x , 'gy' => $y , 'm' => $m , 'x' => $xidx{$x} , 'y' => $yidx{$y}});
    }
    
    my $clusters = $self->_sl_cluster(\@dots, $max_dist, $nb_genes);
    return $clusters;
}

sub _sl_cluster {
    my($self, $dots, $max_dist, $nb_genes) = @_;
    my @color = ( "red", "blue", "green", "turquoise", "violet", "brown", "olivedrab", "purple", "skyblue" );

    # On calcule la matrice de distance
    my @dist;
    foreach my $p1 (@$dots) {
	my @_p;
	foreach my $p2 (@$dots) {
	    push(@_p, { 'p1' => $p1 , 'p2' => $p2 , 'dist' => sqrt(pow($p1->{'x'} - $p2->{'x'}, 2) + pow($p1->{'y'} - $p2->{'y'}, 2)) } );
	}
	push(@dist, \@_p);
    }

    # On fait une slc
    my %vertex = ();
    my $g = Graph::Undirected->new();
    my %edge;
    foreach my $lig (@dist) {
	foreach my $d (@$lig) {
	    if($d->{'dist'} <= $max_dist) {
		my($p1, $p2) = ($d->{'p1'}, $d->{'p2'});
		if(!defined $vertex{$p1}) {
		    $vertex{$p1} = 1;
		    $g->add_vertex($p1);
		    $edge{$p1} = $p1;
		}
		if(!defined $vertex{$p2}) {
		    $vertex{$p2} = 1;
		    $g->add_vertex($p2);
		    $edge{$p2} = $p2;
		}
		$g->add_edge($p1, $p2);
	    }
	}
    }
    my @cc = $g->connected_components();
    my @clusters;
    my $idx = 0;
    foreach my $cl (@cc) { 
	if( scalar(@$cl) > $nb_genes ) {
	    my @dot_list;
	    foreach my $dot (@$cl) {
		my $p = $edge{$dot};
		($p->{'m'})->color($color[$idx]);
		push(@dot_list, $p);
	    }
	    push(@clusters, \@dot_list);
	    $idx++;
	    if($idx == 9) { $idx = 0; }
	}
    }
    return \@clusters;
}


1;


