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

=head1 NAME

Gene - Gene structure

=head1 AUTHORS

2005, Genoscope - CNS, Jean-Marc Aury, jmaury@genoscope.cns.fr

=cut


package LibAnnotation::Ortho::Gene;

use strict;

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {};
    bless $self, $class;
    $self->id(shift || 0);
    $self->seq(shift || "");
    $self->start(shift || 0);
    $self->stop(shift || 0);
    $self->strand(shift || 1);
    $self->axe_pos(shift || 0);
    $self->collection(shift || undef);
    $self->matches([]);
    return $self;
}

sub id {
    my($self, $arg) = @_;
    if(!defined $arg) { return $self->{'id'}; }
    else { $self->{'id'} = $arg; }
}

sub seq {
    my($self, $arg) = @_;
    if(!defined $arg) { return $self->{'seq'}; }
    else { $self->{'seq'} = $arg; }
}

sub start {
    my($self, $arg) = @_;
    if(!defined $arg) { return $self->{'start'}; }
    else { $self->{'start'} = $arg; }
}

sub stop {
    my($self, $arg) = @_;
    if(!defined $arg) { return $self->{'stop'}; }
    else { $self->{'stop'} = $arg; }
}

sub axe_pos {
    my($self, $arg) = @_;
    if(!defined $arg) { return $self->{'axe_pos'}; }
    else { $self->{'axe_pos'} = $arg; }
}

sub strand {
    my($self, $arg) = @_;
    if(!defined $arg) { return $self->{'strand'}; }
    else { $self->{'strand'} = $arg; }
}

sub collection {
    my($self, $arg) = @_;
    if(!defined $arg) { return $self->{'collection'}; }
    else { $self->{'collection'} = $arg; }
}

sub matches {
    my($self, $arg) = @_;
    if(!defined $arg) { return $self->{'matches'}; }
    else { 
        $self->{'matches'} = $arg; 
    }
}

sub addmatch {
    my($self, $match) = @_;
    push(@{$self->matches()}, $match);
}

sub bestMatch {
    my $self = shift;
    if($self->hasMatches()) {
        my @sort = sort { $a->score() <=> $b->score() } @{$self->matches()};
	my $best = pop(@sort);
	return ($best->gene1() == $self) ? $best->gene2() : $best->gene1();
    }
    return undef;
}

sub hasMatches {
    my $self = shift;
    return (scalar(@{$self->{'matches'}}) > 0);
}

sub BRH_match {
    my $self = shift;
    my $best = $self->bestMatch();
    if(defined $best) {
        my $bestReciproq = $best->bestMatch();
        if( defined $bestReciproq && $bestReciproq == $self) {
            return $best;
        }
    }
    return undef;
}

package LibAnnotation::Ortho::GeneHit;

use strict;

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {};
    bless $self, $class;
    $self->gene1(shift || undef);
    $self->gene2(shift || undef);
    $self->score(shift || undef);
    $self->color(shift || "none");
    return $self;
}

sub gene1 {
    my($self, $arg) = @_;
    if(defined $arg) { $self->{'gene1'} = $arg; }
    else { return $self->{'gene1'}; }
}

sub gene2 {
    my($self, $arg) = @_;
    if(defined $arg) { $self->{'gene2'} = $arg; }
    else { return $self->{'gene2'}; }
}

sub score {
    my($self, $arg) = @_;
    if(defined $arg) { $self->{'score'} = $arg; }
    else { return $self->{'score'}; }
}

sub color {
    my($self, $arg) = @_;
    if(defined $arg) { $self->{'color'} = $arg; }
    else { return $self->{'color'}; }
}

1;
