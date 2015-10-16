#!/usr/bin/perl -I../lib

#=============================================================================#
#=             Copyright (c) 2004 by INRA. All rights reserved.              =#
#=                 Redistribution is not permitted without                   =#
#=                 the express written permission of INRA.                   =#
#=                   Mail : eugene@ossau.toulouse.inra.fr                    =#
#=---------------------------------------------------------------------------=#
#= File         : joiner.pl                                                  =#
#= Description  : Build a multigenic sequence                                =#
#= Authors      : P.Bardou, M.J.Cros, S.Foissac, J.Gouzy, A.Moisan, T.Schiex =#
#=============================================================================#

BEGIN
{
    my ($dirprg,$nameprg) = $0 =~ /(.+)\/(.+)/;
    push(@INC,$dirprg);
}

use strict;
use warnings;
use Getopt::Long;
use ParamParser;
use IO::Handle;

=pod

=head1

=head1 NAME

  egn_joiner - Build a multigenic sequence and the corresponding annotation files

=head1 MANDATORY

  -f, --fof  [FILE NAME]
	File of files. Each line must contain a fasta file, an annotation file
	(gff), and optional files (gff).

  -i, --intergenic  [FILE NAME]
	Intergenic region (fasta).

=head1 OPTIONAL

  --utr5  [FILE NAME]
	UTR5 forward region (fasta)

  --utr3  [FILE NAME]
	UTR3 forward region (fasta).

  --minIntergLen  [INT]
	Minimal length of the intergenic region.
	Optional, default is 100.

  --maxIntergLen  [INT]
	Maxiamal length of the intergenic region.
	Optional, default is 5000.

  --utr5Len  [INT]
	Length of the UTR5 region if is missing in the annotation file (gff).
	Optional, default is 130.

  --utr3Len  [INT]
	Length of the UTR3 region if is missing in the annotation file (gff).
	Optional, default is 250.

  -p, --prefix  [FILE NAME]
        output prefix name.
	Optional, default is "seqs".

=head1 Copyright (c) 2004 by INRA. All rights reserved.

=head1 Mail : eugene@ossau.toulouse.inra.fr

=head1

=cut


# Default values
my $MININTERGLEN = 100;
my $MAXINTERGLEN = 5000;
my $UTR5LEN      = 130;
my $UTR3LEN      = 250;
my $PREFIX       = "seq";

# Global
my @a_catSeq;

MAIN:
{
    ## Parse and verify paramaters ##
    my $rh_param = New ParamParser();
    $rh_param->SetUsage(my $usage=sub { &Usage(); } );
    $rh_param->SetBehaviour('exit_on_getopt_error');
    $rh_param->SetBehaviour('assert_empty_file_allowed');
    $rh_param->Update("GETOPTLONG", "I",
		      ("fof=s", "intergenic=s", "utr5=s", "utr3=s", "minIntergLen=i",
		       "maxIntergLen=i", "utr5Len=i", "utr3Len=i", "prefix=s"));

    $rh_param->AssertFileExists('fof');
    $rh_param->AssertFileExists('intergenic');
    my $fof        = $rh_param->Get('fof');
    my $intergfile = $rh_param->Get('intergenic');

    if($rh_param->IsDefined('prefix'))       { $PREFIX       = $rh_param->Get('prefix');       }
    if($rh_param->IsDefined('minIntergLen')) { $MININTERGLEN = $rh_param->Get('minIntergLen'); }
    if($rh_param->IsDefined('maxIntergLen')) { $MAXINTERGLEN = $rh_param->Get('maxIntergLen'); }
    if($rh_param->IsDefined('utr5Len'))      { $UTR5LEN      = $rh_param->Get('utr5Len');      }
    if($rh_param->IsDefined('utr3Len'))      { $UTR3LEN      = $rh_param->Get('utr3Len');      }

    my ($utr5seq, $utr5seqR, $utr3seq, $utr3seqR) = ("", "", "", "");
    if($rh_param->IsDefined('utr5'))
    {
	$utr5seq  = substr(&GetSeq($rh_param->Get('utr5')), 0, $UTR5LEN);
	$utr5seqR = $utr5seq;
	$utr5seqR =~ tr/ATGC/TACG/;
	$utr5seqR = reverse($utr5seqR);
    }
    if($rh_param->IsDefined('utr3'))
    {
	$utr3seq = substr(&GetSeq($rh_param->Get('utr3')), 0, $UTR3LEN);
   	$utr3seqR = $utr3seq;
	$utr3seqR =~ tr/ATGC/TACG/;
	$utr3seqR = reverse($utr3seqR);
    }

    open(FOF, $fof) || die "ERROR, unable to open >$fof<\n";
    my $prevend = $MAXINTERGLEN - 1;
    my %h_gff   = ("gff" => "");
    my $coord   = "";
    my $seq     = "";
    while (my $fofline = <FOF>)
    {
	chomp $fofline;
	my @a_fofline = split (/\s+/, $fofline);
	if(@a_fofline < 2) { &Usage(); exit; }

	my ($startgene, $startcds, $endcds, $endgene, $strand) = &GetLimits($a_fofline[1]);
	my ($utr1, $utr2) = ("", "");

	# Patch because of bug in eugeneAdapt.pl.....
	# The end coord of the (right) utr can be > at the length of the seq.
	if (length(&GetSeq($a_fofline[0])) < $endgene) { $endgene = length(&GetSeq($a_fofline[0])) - 1; }

	my $lengene = $endgene  - $startgene + 1;
	my $lenutr1 = $startcds - $startgene + 1;

	if (($utr5seq eq "" || $utr3seq eq "") && ($startgene == $startcds || $endgene == $endcds))
	{
	    &Usage();
	    print "\nNo UTR coordinates found in >".$a_fofline[1]."<\n";
	    print "You must provide --utr5 and --utr3.\n\n";
	    exit;
	}

	if ($startgene == $startcds)
	{
	    if ($strand eq "+") { $utr1 = $utr5seq;  }
	    else                { $utr1 = $utr3seqR; }
	    $lenutr1 = length($utr1) + 1;
	}
	if ($endgene == $endcds)
	{
	    if ($strand eq "+") { $utr2 = $utr3seq;  }
	    else  	        { $utr2 = $utr5seqR; }
	}

	my $interglen = &ComputeIntergLen($prevend, $lenutr1);
	my $intergseq = &BuildIntergSeq($intergfile, $interglen);
	$seq .= $intergseq.$utr1.(substr(&GetSeq($a_fofline[0]), $startgene, $lengene)).$utr2;
	my $modif = ($prevend - $MAXINTERGLEN) + ($interglen + $lenutr1) - $startcds;
	$coord        .= &ComputeCOORD($a_fofline[1], $modif);
 	$h_gff{"gff"} .= &ComputeGFF  ($a_fofline[1], $modif);
	
	for (my $i=2; $i<=$#a_fofline; $i++)
	{
	    my @a_ext = split (/\./, $a_fofline[$i]);
	    my $ext   = "";
	    if ($a_ext[-1] ne "gff") { $ext = $a_ext[-1].".gff"; }
	    else                     { $ext = $a_ext[-2].".gff"; }

	    $h_gff{$ext} .= &ComputeGFF($a_fofline[$i], $modif);
	}
	$prevend += $interglen + $lengene + length($utr1) + length($utr2);
    }
    $seq .= &BuildIntergSeq($intergfile, $MININTERGLEN);
    close FOF;

    ## OUTPUT ##
    # Multigenic fasta sequence
    open(SEQ, ">$PREFIX.fasta") || die "ERROR, unable to create >$PREFIX.fasta<\n";
    print SEQ ">$fof\n$seq\n";
    close SEQ;

    # GFF files
    foreach my $k (keys(%h_gff))
    {
	open(GFF, ">$PREFIX.fasta.$k") || die "ERROR, unable to create >$PREFIX.$k<\n";
	print GFF $h_gff{$k};
	close GFF;
    }
    close GFF;

    # Coord.txt file for evalpred
    open(COORD, ">$PREFIX.coord.txt") || die "ERROR, unable to create >$PREFIX.coord.txt<\n";
    print COORD "$coord\n";
    close COORD;
}


=head1 SUBROUTINES

=cut

#-----------------------------------------------------------------------------#

=head2 Procedure Usage

 Title		: Print usage
 Usage		: &Usage()
 Function	: Print usage on stdout

=cut

sub Usage
{
    my @usage = `pod2text $0`;
    foreach my $u ( @usage )
    {
	if ($u =~ /Mail :/) { last; }
	else                { print "$u"; }
    }
}

#-----------------------------------------------------------------------------#

=head2 Function GetLimits

 Title		: Get limits : gene start, cds start, cds end, gene end.
 Usage		: &GetLimits()
 Function	: Extract limit coords for the coord file (gff).
 Returns	: gene start, cds start, cds end, gene end and strand.
		  -> WITH gene start <= cds start < cds end <= gene end
		  -> If no UTR geneStart = cdsStart and/or cdsEnd = geneEnd.
 Args		: Gff file name
 Globals	: none

=cut

sub GetLimits
{
    my ($gfffile) = @_;
    my ($startgene, $startcds, $endcds, $endgene) = (-1, -1, -1, -1);
    my $strand = "";

    open(GFF, $gfffile) || die "ERROR, unable to open >".$gfffile."<\n";
    while (my $l = <GFF>)
    {
	chomp $l;
	if ( $l =~ /^#/ ) { next; }
	my @a_splitL = split (/\s+/, $l);
	if($l =~ /UTR5/)
	{
	    if($l =~ /\+/) { $startgene = $a_splitL[3]; }
	    else           { $endgene   = $a_splitL[4]; }
	}
	if($l =~ /UTR3/)
	{
	    if($l =~ /\+/) { $endgene   = $a_splitL[4]; }
	    else           { $startgene = $a_splitL[3]; }
	}
	if($l =~ /E.Init/  ||  $l =~ /E.Sngl/)
	{
	    if($l =~ /\+/) { $startcds  = $a_splitL[3]; }
	    else           { $endcds    = $a_splitL[4]; }
	}
	if($l =~ /E.Term/  ||  $l =~ /E.Sngl/)
	{
	    if($l =~ /\+/) { $endcds    = $a_splitL[4]; }
	    else           { $startcds  = $a_splitL[3]; }
	}
	$strand = $a_splitL[6];
    }
    close GFF;
    if ($startgene == -1) { $startgene = $startcds; }
    if ($endgene   == -1) { $endgene   = $endcds;   }

    return ($startgene, $startcds, $endcds, $endgene, $strand);
}

#-----------------------------------------------------------------------------#

=head2 Function ComputeIntergLen

 Title		: Compute the intergenic lenght.
 Usage		: &ComputeIntergLen()
 Function	: Compute the intergenic lenght to add between two genes.
 Returns	: intergenic lenght.
 Args		: $prevend (coord of the end of the previous gene),
		  $lenutr
 Globals	: $MININTERGLEN, $MAXINTERGLEN

=cut

sub ComputeIntergLen
{
    my ($prevend, $lenutr) = @_;

    my $interglen = $MAXINTERGLEN - ($prevend % $MAXINTERGLEN) - $lenutr;
    while ($interglen < $MININTERGLEN) { $interglen += $MAXINTERGLEN; }
    return $interglen;
}

#-----------------------------------------------------------------------------#

=head2 Function BuildIntergSeq

 Title		: Build intergenic sequence.
 Usage		: &BuildIntergSeq()
 Function	: Build intergenic sequence of LEN nucleotides from a file.
 Returns	: sequence.
 Args		: Fasta filename, len.
 Globals	: none

=cut

sub BuildIntergSeq
{
    my ($seqfile, $len) = @_;

    my $seq = "";
    while (length($seq) < $len) { $seq .= &GetSeq($seqfile); }

    return (substr($seq, 0, $len));
}

#-----------------------------------------------------------------------------#

=head2 Function GetSeq

 Title          : Get sequence
 Usage          : &GetSeq()
 Prerequisite   : none
 Returns        : sequence
 Args           : fasta filename
 Globals        : none

=cut

sub GetSeq
{
    my($seq) = @_;
    my($lign);
    my($seq_seq) = "";
    my($multifasta) = -1;

    open(SEQ,"$seq") or die "ERROR, unable to open >$seq<\n";
    while($lign=<SEQ>)
    {
        chomp($lign);
        if ( $lign =~ /^>/ ) { $multifasta++;     }
        else                 { $seq_seq .= $lign; }
        if ( $multifasta )
        {
            print STDERR "FATAL ERROR: this script cannot be used".
                " on a multifasta file ($seq)\n";
            exit;
        }
    }
    close(SEQ);
    $seq_seq =~ s/[^ATGCNX]//gi;
    return "\U$seq_seq\E";
}

#-----------------------------------------------------------------------------#

=head2 Procedure ComputeGFF

 Title		: Compute new coordinates.
 Usage		: &ComputeGFF()
 Function	: Compute new coordinates.
 Returns	: $   ...lines of the new gff file
 Args		: coord filename (gff), lenght before the gene.
 Globals	: none

=cut

sub ComputeGFF
{
    my ($gfffile, $len) = @_;
    my $c = "";

    open(GFF,"$gfffile") or die "ERROR, unable to open >$gfffile<\n";
    while(my $lign=<GFF>)
    {
	chomp $lign;
	if ( $lign =~ /^#/ ) { next; }

	my @a_splitL = split (/\s+/, $lign);
	$a_splitL[3] += $len;
	$a_splitL[4] += $len;
	
	$c .= join("\t", @a_splitL);
	$c .= "\n";
    }
    close GFF;

    return $c;
}

#-----------------------------------------------------------------------------#

=head2 Procedure ComputeCOORD

 Title		: Compute new coordinates.
 Usage		: &ComputeCOORD()
 Function	: Compute new coordinates.
 Returns	: $   ...lines of coord.txt for evalpred
 Args		: coord filename (gff), lenght before the gene.
 Globals	: none

=cut

sub ComputeCOORD
{
    my ($gfffile, $len) = @_;
    my $c = "";

    open(GFF,"$gfffile") or die "ERROR, unable to open >$gfffile<\n";
    while(my $lign=<GFF>)
    {
	chomp $lign;
	if ( $lign !~ /E\./ ) { next; }

	my @a_splitL = split (/\s+/, $lign);
	if ($c eq "") { $c = $a_splitL[0]." "; }
	$a_splitL[3] += $len;
	$a_splitL[4] += $len;
	if ($a_splitL[6] eq '+') {$a_splitL[6] = ""; }

	$c .= $a_splitL[6].$a_splitL[3]." ".$a_splitL[6].$a_splitL[4]." ";
    }
    $c .= "\n";
    close GFF;

    return $c;
}

#-----------------------------------------------------------------------------#
