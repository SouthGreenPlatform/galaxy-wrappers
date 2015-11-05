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

GeneCollection - Structure 

=head1 AUTHORS

2005, Genoscope - CNS, Jean-Marc Aury, jmaury@genoscope.cns.fr

=cut

package LibAnnotation::Ortho::GeneCollection;

use strict;
use FileHandle;
use LibAnnotation::Ortho::Gene;
use LibAnnotation::Ortho::GeneCollectionIterator;

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {};
    bless $self, $class;
    $self->file(shift);
    $self->{'genes'} = [];
    $self->{'indexes'} = {};
    $self->{'axeSize'} = 0;
    $self->{'lastIdx'} = 0;
    $self->{'existGene'} = {};
    $self->{'name'} = "";
    $self->{'cache'} = undef;
    $self->{'seqInCache'} = "";
    return $self;
}

sub dup {
    my($self, $seqFilter) = @_;
    my $class = ref($self);
    my $collec = {};
    bless $collec, $class;
    my $it = new LibAnnotation::Ortho::GeneCollectionIterator($self);
    my @genelist = ();
    while(my $el = $it->next()) {
	if(defined $seqFilter->{$el->seq()}) { push(@genelist, $el); }
    }
    $collec->genes(\@genelist);
    $collec->name($self->name());
    return $collec;
}

sub file {
    my($self, $arg) = @_;
    if(defined $arg) { $self->{'file'} = $arg; }
    else { return $self->{'file'}; }
}

sub genes {
    my ($self, $arg) = @_;
    if(defined $arg) { $self->{'genes'} = $arg; }
    else { return $self->{'genes'}; }
}

sub indexes {
    my($self, $arg) = @_;
    if(defined $arg) { $self->{'indexes'} = $arg; }
    else { return $self->{'indexes'}; }
}

sub axeSize {
    my($self, $arg) = @_;
    if(defined $arg) { $self->{'axeSize'} = $arg; }
    else { return $self->{'axeSize'}; }
}

sub lastIdx {
    my($self, $arg) = @_;
    if(defined $arg) { $self->{'lastIdx'} = $arg; }
    else { return $self->{'lastIdx'}; }
}

sub name {
    my($self, $arg) = @_;
    if(defined $arg) { $self->{'name'} = $arg; }
    else { return $self->{'name'}; }
}

sub existGene {
    my($self, $arg) = @_;
    if(defined $self->{'existGene'}->{$arg}) { return $self->{'existGene'}->{$arg}; }
    else { return 0; }
}

sub gene {
    my($self, $id) = @_;
    my $indexes = $self->indexes();
    if(defined $indexes->{$id}) { return $indexes->{$id}; }
    else { return undef; }
}

sub getGene {
    my ($self, $id) = @_;
    return ($self->indexes())->{$id};
}

sub addGene {
    my ($self, $seq, $start, $end) = @_;
    my $exist = $self->existGene($seq."@".$start."@".$end);
    if($exist) { return $exist; }
    my $gene = new LibAnnotation::Ortho::Gene($self->lastIdx(), $seq, $start, $end);
    $self->lastIdx($self->lastIdx() + 1);
    $gene->collection($self);
    push(@{$self->genes()}, $gene);
    $self->{'existGene'}->{$seq."@".$start."@".$end} = $gene;
    return $gene;
}

sub loadFromFile {
    my $self = shift;

    my $fhIndice = $self->_openFile();
    my $cpt = 0;
    my @geneSet;
    while(my $line = $fhIndice->getline()) {
        if($line =~ /^(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+([\+|\-|\.]?)/) {
            my $gene = new LibAnnotation::Ortho::Gene($1, $2, $3, $4, $5);
            $gene->collection($self);
            push(@geneSet, $gene);
            $cpt++;
        } else {
            warn "[Error] Bad format in geneSet file : $line\n";
            exit 1;
        }
    }
    $self->genes(\@geneSet);
    $self->indexGenes();
    $self->lastIdx($cpt);
    return $cpt;
}

sub loadFromHits {
    my ($self, $format, $collection) = @_;

    my $fhIndice = $self->_openFile();
    if($format == 1) {
	$self->_loadFromHits_simpleEntry($fhIndice, $collection);
    } elsif($format == 2) {
	$self->_loadFromHits_fullEntry($fhIndice, $collection);
    }
}

sub loadFromHits_simpleEntry {
    my($self, $fh, $collection) = @_;
    my $id = 1;
    my ($cpt1, $cpt2, $nbLink) = (0, 0, 0);
    my (@geneSet1, @geneSet2);
    while(my $line = $fh->getline()) {
        if($line =~ /^(\S+)\s+(\d+)\s+(\d+)\s+([\+|\-|\.]?)\s*(\S+)\s+(\d+)\s+(\d+)\s+([\+|\-|\.]?)$/) {
	    my($gene1, $gene2);
	    if($1 ne ".") {
		$gene1 = new LibAnnotation::Ortho::Gene($id++, $1, $2, $3, $4);
		$gene1->collection($self);
		push(@geneSet1, $gene1);
		$cpt1++;
	    }
	    if($5 ne ".") {
		$gene2 = new LibAnnotation::Ortho::Gene($id++, $5, $6, $7, $8);
		push(@geneSet2, $gene2);
		$gene2->collection($collection);
		$cpt2++;
		if($1 ne ".") {
		    $nbLink++;
		    my $hit = new LibAnnotation::Ortho::GeneHit($gene1, $gene2, 0);
		    $gene1->addmatch($hit);
		    $gene2->addmatch($hit);
		}
	    }
	    
        } else {
            warn "[Error] Bad format in geneHit file : $line\n";
            exit 1;
        }
    }
    $self->genes(\@geneSet1);
    $self->indexGenes();
    $self->lastIdx($cpt1);
    $collection->genes(\@geneSet2);
    $self->indexGenes();
    return ($cpt1, $cpt2, $nbLink);
}

sub loadFromHits_fullEntry {
    my($self, $fh, $collection) = @_;
    my ($cpt1, $cpt2, $nbLink) = (0, 0, 0);
    my (@geneSet1, @geneSet2);
    while(my $line = $fh->getline()) {
        if($line =~ /^(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+([\+|\-|\.])\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+([\+|\-|\.])\s+(\d+)$/) {
	    my($gene1, $gene2);
	    if($1 ne ".") {
		$gene1 = new LibAnnotation::Ortho::Gene($1, $2, $3, $4, $5);
		$gene1->collection($self);
		push(@geneSet1, $gene1);
		$cpt1++;
	    }
	    if($5 ne ".") {
		$gene2 = new LibAnnotation::Ortho::Gene($6, $7, $8, $9, $10);
		push(@geneSet2, $gene2);
		$gene2->collection($collection);
		$cpt2++;
		if($1 ne ".") {
		    $nbLink++;
		    my $hit = new LibAnnotation::Ortho::GeneHit($gene1, $gene2, $11);
		    $gene1->addmatch($hit);
		    $gene2->addmatch($hit);
		}
	    }
	    
        } else {
            warn "[Error] Bad format in geneHit file : $line\n";
            exit 1;
        }
    }
    $self->genes(\@geneSet1);
    $self->indexGenes();
    $self->lastIdx($cpt1);
    $collection->genes(\@geneSet2);
    $self->indexGenes();
    return ($cpt1, $cpt2, $nbLink);
}

sub _openFile {
    my $self = shift;
    my $fhIndice = new FileHandle($self->file(), "r");
    if(! defined $fhIndice) { die "[Error] Unable to open ", $self->file(), "\n"; }
    return $fhIndice;
}

sub indexGenes {
    my $self = shift;
    my %geneIndices;
    $self->sortGenes();
    foreach my $g (@{$self->genes()}) { $geneIndices{$g->id()} = $g; }
    $self->indexes(\%geneIndices);
}

sub sortGenes {
    my $self = shift;
    my @sortList = sort { ($a->seq() =~ /^\d+$/ && $b->seq() =~ /^\d+$/ && $a->seq() <=> $b->seq())
			      || $a->seq() cmp $b->seq()
#			      || (($a->seq() !~ /\d+/ || $b->seq() !~ /\d+/) && $a->seq() cmp $b->seq())
			      || $a->start() <=> $b->start() } @{$self->genes()};
    $self->genes(\@sortList);
}

sub projectOnAxis {
    my($self) = @_;
    my $nbGene = scalar(@{$self->genes()});
    my $offset = $self->axeSize() / $nbGene;
    $self->sortGenes();
    my $it = new LibAnnotation::Ortho::GeneCollectionIterator($self);
    my ($index, $currentS) = (0, "");
    my (%st, %pos);
    while(my $el = $it->next()) {
	if($el->seq() ne $currentS) {
	    if($currentS ne "") {
		$st{$currentS}->{'end'} = $index * $offset;
	    }
	    $currentS = $el->seq();
	    $st{$currentS}->{'start'} = $index * $offset;
	}
	$pos{$el->seq()."@".$el->start()."@".$el->stop()} = $index * $offset;
	$el->axe_pos($index * $offset);
	$index++;
    }
    $st{$currentS}->{'end'} = $index * $offset;
    return (\%st, \%pos);;
}    

sub getSeqNames {
    my($self) = @_;   
    my $nbGene = scalar(@{$self->genes()});
    $self->sortGenes();
    my $it = new LibAnnotation::Ortho::GeneCollectionIterator($self);
    my $currentS = "";
    my @seq = ();
    while(my $el = $it->next()) {
	if($el->seq() ne $currentS) {
	    if($currentS ne "") {
		push(@seq, $currentS);
	    }
	    $currentS = $el->seq();
	}
    }
    if($currentS ne "") { push(@seq, $currentS); }
    return \@seq;
}

sub getGenesFromSeq {
    my($self, $seqx, $seqy) = @_;
    if(!defined $seqy) { $seqy = ""; }

    my (@genes, @matches);
    if(defined $self->{'seqInCache'} && $self->{'seqInCache'} eq $seqx) {
	@genes = @{$self->{'cache'}};
    } else {
	my $it = new LibAnnotation::Ortho::GeneCollectionIterator($self);
	while(my $el = $it->next()) {
	    if($el->seq() eq $seqx) {
		push(@genes, $el);
	    }
	}
	$self->{'seqInCache'} = $seqx;
	$self->{'cache'} = \@genes;
    }
    
    my $unique_matchs = 0;
    my %vu;
    if($seqy ne "") {
	foreach my $el (@genes) {
	    my $matches = $el->matches();
	    foreach my $m (@$matches) {
		my $ortho_gene = ($m->gene1() == $el) ? $m->gene2() : $m->gene1();
		# On ne recupere que les matchs qui correspondent a la sequence $seqy
		if($ortho_gene->seq() eq $seqy) {
		    if(!defined $vu{$el}) { $vu{$el}=1; $unique_matchs++; }
		    push(@matches, $m); 
		}
	    }
	}
    }   
    return (\@genes, \@matches, $unique_matchs);
}

sub get_indexBased_coord {
    my($self) = @_; 
    my ($index, $currentS) = (1, "");
    my (%st, %pos);
    my $it = new LibAnnotation::Ortho::GeneCollectionIterator($self);
    while(my $el = $it->next()) {
	if($el->seq() ne $currentS) {
	    if($currentS ne "") {
		$st{$currentS}->{'nb'} = $index;
		$index = 1;
	    }
	    $currentS = $el->seq();
	    #$st{$currentS}->{'start'} = $index;
	}
	push(@{$st{$currentS}->{'genes'}}, $el);
	$pos{$el} = $index;
	$index++;
    }
    #$st{$currentS}->{'end'} = $index;
    return (\%st, \%pos);;
}

1;
