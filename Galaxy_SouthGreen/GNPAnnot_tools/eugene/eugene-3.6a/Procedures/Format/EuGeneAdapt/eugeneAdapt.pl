#!/usr/bin/perl

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
# $Id: eugeneAdapt.pl,v 1.7 2004-09-16 11:56:35 cros Exp $
# ------------------------------------------------------------------
# File:     eugeneAdapt.pl
# Contents: creates the different files need for an adaptation of
# eugene to a new specie
# Phil: Attention aux lignes UTRs dans les fichiers gff.
#         Ne sont notés que les UTRs ne contenant pas d'intron. (méthode =
#         alignement de chaque UTR contre le génomique (sim4) et parsing de
#         la sortie Cf. Sim4ParsingUTR.)
#         Du coup concernant les longueurs "#LenGene" (Exons+Introns+UTRs) et
#         "#LencDNA" (Exons+UTRs) elles sont mentionnées "aux UTRs près"...!
#
# First step for the adaptation of eugene to a new specie.
# History      : version 1.0 (Dec. 1, 2002)  -> sim2eugene.pl
#                version 2.0 (Aou. 5, 2003)  -> eugeneAdapt.pl
# Warning :
#  - This script doesn't accept redondances in each list. Thus, if one genomic
#    sequence "contains" more than one cDNA (2), or if the genomic and the gff
#    file contains more than one gene, it has to be splited.
#    The alphabetical order allow to lauch eugene using *.fasta and obtain the
#    predictions must be in the same order than the list files.
#  - For the second approch, each cDNA has to be the real cDNA corresponding
#    to the gene present in his genomic sequence, a good sequence quality is
#    request (> ~95%), in order to obtain with the sim4 alignment the true
#    exons coordinates.
#  - eugeneAdapt need :
#       -> EUGENEDIR        : root directory of the eugene distribution
#       -> two EMBOSS tools : extractseq and revseq.
#       -> Eu-imm2UTR       : Matrices builder
#       -> WAMbuilder       : WAM builder
#       -> seqlogo          : Plot signal consensus (WAM)
#       -> For the second approch sim4 and GeneSeqer.
#  - eugeneAdapt creates temporary files deleted at the end of procedure.
# ------------------------------------------------------------------

use strict;
use warnings;
use IO::Handle;
use Getopt::Long;
use gffUtils;


=pod

=head1

=head1 DESCRIPTION

  eugeneAdapt creates the different files need for an adaptation of
  eugene to a new specie.
  Two approaches could be used. The first one (1) is based  on  the
  gene  coordinates  knowledge about  genomic  sequences.  And  the
  second (2)  one  is  based  on  sim4  program (cDNA  and  Genomic
  DNA Alignment).

=head1 OPTIONS

  -g, --genomic [FILE NAME] (1) & (2)
       File containing a list of the fasta genomic sequence file(s)
       sorted by alphabetical order.

  -c, --coord [FILE NAME] (1)

  -m, --mRna [FILE NAME] (2)
       File containing a list of the  fasta  cDNA  file(s), a  cDNA
       from a line corresponding to the  DNA  from  the  same  line
       number in the file genomic list.

  -s, --ssCoord [FILE NAME] (2)
       File containing  ATG  and  STOP for each cDNA,  in  the same
       order than the others lists (separated by a tabulation).

  -i, --input [FILE NAME] (2)
       Each line must contain:
             genomic file name / cDNA file name / Start / Stop
          or genomic file name / cDNA file name
       In the second case the cDNA must start by ATG and finish  by
       a stop codon.

  -a, --auto (2)
       Optional: For  an  automatic  process   (each  sequence   is
       automatically excluded).

  -o, --offset [INT] (1) & (2)
       Optional: offset (default = 1000) max number  of  nucleotids
       on each side of the START and STOP position.

  -w, --wam [INT] (1) & (2)
       Optional: build  the  WAM  files  and  test  the  new  fasta
       sequences.

  -d, --dirOutput [DIR NAME] (1) & (2)
       Optional: output directory (default = EuAdaptOutput).

=head1 EXAMPLES

  (1) eugeneAdapt -g genoListFile -c coordFile -w
  (2) eugeneAdapt -g genoListFile -m cDnaListFile -s ATGStopFile -w
  (2) eugeneAdapt -i inputsFile   -w

=head1 Copyright (c) 2003 by INRA. All rights reserved.

=head1 Mail : tschiex@toulouse.inra.fr

=head1

=cut

# $EUGENEDIR test
my $EUGENEDIR = $ENV{EUGENEDIR};
if (!defined $EUGENEDIR)
  {
    print "\n\$EUGENEDIR (root directory of the eugene distribution) must".
      " be defined.\n\n";
    exit();
  }

# Program used
my $cmd_sim4       = "sim4";          # cDNA and Genomic DNA Alignment
my $cmd_gseqer     = "GeneSeqer";     # cDNA and Genomic DNA Alignment
my $cmd_extractseq = "extractseq";    # EMBOSS
my $cmd_revseq     = "revseq";        # EMBOSS

# Matrices builder
my $cmd_IMM =
  "$EUGENEDIR/SensorPlugins/MarkovIMM/GetData/TrainIMM -h";

# WAM builder
my $cmd_WAMBuilder =
  "$EUGENEDIR/SensorPlugins/0_SensorTk/GetData/WAMbuilder";

# seqlogo
my $cmd_seqlogo =
  "perl -I$EUGENEDIR/../Procedures/Format/EuGeneAdapt $EUGENEDIR/../Procedures/Format/EuGeneAdapt/seqlogo";

# For the ouput files (extension or file name)
my $sim4_out_ext   = ".sim4.A3";
my $sim4_out_ext2  = ".sim4.UTR5.A3";
my $sim4_out_ext3  = ".sim4.UTR3.A3";
my $gseqer_out_ext = ".gseqer";
my $newdna_ext    = ".fasta";
my $gff_ext       = ".fasta.gff";
my $info_ext      = ".fasta.info";
my $coordFile     = "coord.txt";
my $exonFile      = "exon.txt";
my $intFile       = "intron.txt";
my $utr5File      = "utr5.txt";
my $utr3File      = "utr3.txt";
my $matFile       = "matrice.mat";
my $matLogFile    = "matrice.log";
my $wam_sta_fp    = "WAM.START.FP";
my $wam_don_fp    = "WAM.DON.FP";
my $wam_acc_fp    = "WAM.ACC.FP";
my $wam_sto_fp    = "WAM.STOP.FP";
my $wam_sta_tp    = "WAM.START.TP";
my $wam_don_tp    = "WAM.DON.TP";
my $wam_acc_tp    = "WAM.ACC.TP";
my $wam_sto_tp    = "WAM.STOP.TP";
my $verif         = "verif.txt";
my $eugAdaptLog   = "EUGADAPT.LOG";

# For print info
my $seqInEALog   = 9;         # first method (seq in ealog when exon < N nt)
my $nbExclu      = 0;         # number of excluded sequence
my $nbNoUTR5     = 0;         # number of sequence without utr5
my $nbNoUTR3     = 0;         # number of sequence without utr3
my $nbNoIntron   = 0;         # number of sequence without intron
my $nbSeq        = 0;         # number of sequence
my $nbSeqRev     = 0;         # number of reverse sequence
my $exonMinLen   = 100000;    # min length for exon
my $exonMaxLen   = 0;         # max length for exon
my $intronMinLen = 100000;    # min length for intron
my $intronMaxLen = 0;         # max length for intron

# For WAM
my $worder  = 3;              # WAM order
my $wamstaB = 10;              # Context lenght, B before signal
my $wamstaA = 10;              # Context lenght, A after  signal
my $wamdonB = 10;              # Context lenght, B before signal
my $wamdonA = 10;              # Context lenght, A after  signal
my $wamaccB = 10;              # Context lenght, B before signal
my $wamaccA = 10;              # Context lenght, A after  signal
my $wamstoB = 10;              # Context lenght, B before signal
my $wamstoA = 10;              # Context lenght, A after  signal

# Options
my $genomic;
my $cDna;
my $ssCoord;
my $coord;
my $input;
my $auto;
my $offset    = 1000;
my $outputDir = "EuAdaptOutput";
my $wam;

# For the -i option (2)
my $i_geno = "tmp_geno";
my $i_cdna = "tmp_cdna";
my $i_ss   = "tmp_ss";

MAIN:
{
    if (
        !GetOptions(
                    "genomic=s"   => \$genomic,
                    "mRna=s"      => \$cDna,
                    "ssCoord=s"   => \$ssCoord,
                    "coord=s"     => \$coord,
                    "input=s"     => \$input,
		    "auto"        => \$auto,
		    "offset=i"    => \$offset,
                    "wam"         => \$wam,
                    "dirOutput=s" => \$outputDir
        )
      )
    {
        exit(1);
    }

    ## Verify the options ##
    OptVerif();

    ## Open output files ##
    my $error = "no";
    open(COORD_OUT, ">$coordFile")  || die $error = $coordFile;
    open(EXON,      ">$exonFile")   || die $error = $exonFile;
    open(INTRON,    ">$intFile")    || die $error = $intFile;
    open(UTR5,      ">$utr5File")   || die $error = $utr5File;
    open(UTR3,      ">$utr3File")   || die $error = $utr3File;
    open(EALOG,     ">$eugAdaptLog")|| die $error = $eugAdaptLog;
    if (defined($wam))
    {
        open(FPSTART, ">$wam_sta_fp") || die $error = $wam_sta_fp;
        open(FPDON,   ">$wam_don_fp") || die $error = $wam_don_fp;
        open(FPACC,   ">$wam_acc_fp") || die $error = $wam_acc_fp;
        open(FPSTOP,  ">$wam_sto_fp") || die $error = $wam_sto_fp;
        open(TPSTART, ">$wam_sta_tp") || die $error = $wam_sta_tp;
        open(TPDON,   ">$wam_don_tp") || die $error = $wam_don_tp;
        open(TPACC,   ">$wam_acc_tp") || die $error = $wam_acc_tp;
        open(TPSTOP,  ">$wam_sto_tp") || die $error = $wam_sto_tp;
        open(VERIF,   ">$verif")      || die $error = $verif;
    }
    if ($error ne "no") {die "ERROR: Could not open file $error\n";}

    ## Choose the approach
    if (defined($coord)) {GenoCoord();}    # first  method
    else                 {GenoSim4();}     # second method

    ## Build matrices ##
    print "------------------------------------";
    print "-------------------------------\n";
    print "Launch TrainImm (Build matrices)..";
    system "$cmd_IMM $exonFile $intFile -5 $utr5File -3 $utr3File > $matLogFile";
    system "mv $outputDir/exon.txt.bin $matFile";
    print "$matFile....done\n";

    ## Build WAM ##
    if (defined($wam))
    {
        print "------------------------------------";
	print "-------------------------------\n";
        print "Launch WAMBuilder (Build WAM)..";
        print "StartTF..";
        my $arg = "$worder " . ($wamstaB + $wamstaA + 3) . " " . ($worder + 1);
        system "$cmd_WAMBuilder $arg $wam_sta_tp $wam_sta_tp. >& $wam_sta_tp.log";
        system "$cmd_WAMBuilder $arg $wam_sta_fp $wam_sta_fp. >& $wam_sta_fp.log";
        print "DonTF..";
        $arg = "$worder " . ($wamdonB + $wamdonA + 2) . " " . ($worder + 1);
        system "$cmd_WAMBuilder $arg $wam_don_tp $wam_don_tp. >& $wam_don_tp.log";
        system "$cmd_WAMBuilder $arg $wam_don_fp $wam_don_fp. >& $wam_don_fp.log";
        print "AccTF..";
        $arg = "$worder " . ($wamaccB + $wamaccA + 2) . " " . ($worder + 1);
        system "$cmd_WAMBuilder $arg $wam_acc_tp $wam_acc_tp. >& $wam_acc_tp.log";
        system "$cmd_WAMBuilder $arg $wam_acc_fp $wam_acc_fp. >& $wam_acc_fp.log";
        print "StopTF..";
        $arg = "$worder " . ($wamstoB + $wamstoA + 3) . " " . ($worder + 1);
        system "$cmd_WAMBuilder $arg $wam_sto_tp $wam_sto_tp. >& $wam_sto_tp.log";
        system "$cmd_WAMBuilder $arg $wam_sto_fp $wam_sto_fp. >& $wam_sto_fp.log";
        print ".done\n";

        ## Build logos ##
        print "Launch $cmd_seqlogo (Build logo)....";
        $arg = "-b -c -n -Y -F PNG";
        print "StartTF..";
        system "$cmd_seqlogo $arg -x StartTP -f $wam_sta_tp -o $wam_sta_tp";
        system "$cmd_seqlogo $arg -x StartFP -f $wam_sta_fp -o $wam_sta_fp";
        print "DonTF..";
        system "$cmd_seqlogo $arg -x DonorTP -f $wam_don_tp -o $wam_don_tp";
        system "$cmd_seqlogo $arg -x DonorFP -f $wam_don_fp -o $wam_don_fp";
        print "AccTF..";
        system "$cmd_seqlogo $arg -x AcceptorTP -f $wam_acc_tp -o $wam_acc_tp";
        system "$cmd_seqlogo $arg -x AcceptorFP -f $wam_acc_fp -o $wam_acc_fp";
        print "StopTF..";
        system "$cmd_seqlogo $arg -x StopTP -f $wam_sto_tp -o $wam_sto_tp";
        system "$cmd_seqlogo $arg -x StopFP -f $wam_sto_fp -o $wam_sto_fp";
        print ".done\n";
    }

    PrintInfo();
    FileClose();
}

##################
#     SUB        #
##################
#-----------------------------------------------------------------------------#
sub Usage
{
    system("pod2text $0");
    exit();
}

#-----------------------------------------------------------------------------#
#-  Verify the options and open the files corresponding to the approach.     -#
#-----------------------------------------------------------------------------#
sub OptVerif
{
    # "input" option
    if (defined($input))
    {
      if (!open(INPUTS, $input))
      {
	print STDERR "ERROR: Could not open inputs file $input\n";
	exit(2);
      }
      open(TMP_GENO, ">$i_geno") || die "ERROR: Could not open file $i_geno\n";
      open(TMP_CDNA, ">$i_cdna") || die "ERROR: Could not open file $i_cdna\n";
      open(TMP_SS,   ">$i_ss")   || die "ERROR: Could not open file $i_ss\n";
      while (my $inputsLine = <INPUTS>)
      {
	my @tab = split (/\s+/,$inputsLine);
	if (@tab == 2  ||  @tab == 4)
	  {
	    print TMP_GENO "$tab[0]\n";
	    print TMP_CDNA "$tab[1]\n";
	    if (@tab == 4) { print TMP_SS "$tab[2]\t$tab[3]\n"; }
	    else           { print TMP_SS "\n"; }
	  }
	else { &Usage(); }
      }
      close INPUTS;
      close TMP_GENO;
      close TMP_CDNA;
      close TMP_SS;
      $genomic = $i_geno;
      $cDna    = $i_cdna;
      $ssCoord = $i_ss;
    }

    # "genomic" option
    if (defined($genomic))
    {
        if (!open(GENOMIC, $genomic))
        {
            print STDERR "ERROR: Could not open genomic file list $genomic\n";
            exit(2);
        }
    }
    else {&Usage();}

    # "cDna" "ssCoord" and "coord" options
    if (defined($cDna) && defined($ssCoord) && !defined($coord))
    {
        if (!open(CDNA, $cDna))
        {
            print STDERR "ERROR: Could not open cDNA file list $cDna\n";
            exit(2);
        }
        if (!open(SSCOORD, $ssCoord))
        {
            print STDERR "ERROR: Could not open start-stop list $ssCoord\n";
            exit(2);
        }
    }
    else
    {
        if (defined($coord) && !defined($cDna) && !defined($ssCoord))
        {
            if (!open(COORD_IN, $coord))
            {
                print STDERR "ERROR: Could not open coordinate list $coord\n";
                exit(2);
            }
        }
        else {&Usage();}
    }

    if ($offset < 0) {&Usage();}

    if (defined($outputDir))
    {
        mkdir $outputDir;
        mkdir $outputDir . "/Genes";
        $coordFile  = $outputDir . "/" . $coordFile;
        $exonFile   = $outputDir . "/" . $exonFile;
        $intFile    = $outputDir . "/" . $intFile;
        $utr5File   = $outputDir . "/" . $utr5File;
        $utr3File   = $outputDir . "/" . $utr3File;
        $matFile    = $outputDir . "/" . $matFile;
        $matLogFile = $outputDir . "/" . $matLogFile;
	$eugAdaptLog= $outputDir . "/" . $eugAdaptLog;
        if (defined($wam))
        {
            mkdir $outputDir . "/WAM";
            $wam_sta_fp = $outputDir . "/WAM/" . $wam_sta_fp;
            $wam_don_fp = $outputDir . "/WAM/" . $wam_don_fp;
            $wam_acc_fp = $outputDir . "/WAM/" . $wam_acc_fp;
            $wam_sto_fp = $outputDir . "/WAM/" . $wam_sto_fp;
            $wam_sta_tp = $outputDir . "/WAM/" . $wam_sta_tp;
            $wam_don_tp = $outputDir . "/WAM/" . $wam_don_tp;
            $wam_acc_tp = $outputDir . "/WAM/" . $wam_acc_tp;
            $wam_sto_tp = $outputDir . "/WAM/" . $wam_sto_tp;
            $verif      = $outputDir . "/" . $verif;
        }
    }
}

#-----------------------------------------------------------------------------#
#-  Close files and remove the temporary files                               -#
#-----------------------------------------------------------------------------#
sub FileClose
{

    # Close input files
    close GENOMIC;
    if (defined($cDna)) { close CDNA; close SSCOORD; }
    else                { close COORD_IN; }

    # Close output files
    close EXON;
    close INTRON;
    close UTR5;
    close UTR3;
    close EALOG;
    if (defined($wam))
    {
        close FPSTART;
        close FPDON;
        close FPACC;
        close FPSTOP;
        close TPSTART;
        close TPDON;
        close TPACC;
        close TPSTOP;
        close VERIF;
    }

    # Remove tmp files
    unlink("s2e_t.txt");
    unlink("s2e_t2.txt");
    unlink("fastaREV");
    if (defined($input))
    {
      unlink($i_geno);
      unlink($i_cdna);
      unlink($i_ss);
    }
}

#-----------------------------------------------------------------------------#
#-  First approach.                                                          -#
#-----------------------------------------------------------------------------#
sub GenoCoord
{
    while (my $genofile = <GENOMIC>)
    {
        print "------------------------------------";
	print "-------------------------------\n";
        print($nbSeq+ 1);
        flush STDOUT;
        print STDERR ".";
        $nbSeq++;
        chomp $genofile;

        ## ID extraction ##
        my $id = IDExtract($genofile);
        print " $id :\n";

        ## Read Coord (.gff) file ##
        print"  - Read coordinates file (.gff)...............................";
        my $gfffile = <COORD_IN>;
        chomp $gfffile;
        my @geneCoord = gffUtils::gffGetCoord("$gfffile");
        print "done\n";

        ## Launch extractseq (EMBOSS) : + offset (DNA) ##
        print "  - Launch extractseq (Create new genomic file,";
	print " offset:$offset)...";
        my $begin = ((abs($geneCoord[1]) - $offset > 0)
		     ? abs($geneCoord[1]) - $offset : 1);
        my $end = abs($geneCoord[$#geneCoord - 1]) + $offset;
        ExtractSeq($genofile, $begin, $end,
		   " $outputDir" . "/Genes/" . $id . $newdna_ext);
        print "done\n";

        ## Recomputing coord for the new genomic sequence ##
        my $newFastaFile = "$outputDir/Genes/$id$newdna_ext";
        my $sens = (($geneCoord[1] > 0) ? "" : "-");
        if ($sens eq "-") {$nbSeqRev++;}
        my @newCoord = ();

        print "  - Recompute coordinates :\n";
        print "\t-> Write in coord file...........................";
        my $newcoord = 0;
        my $modif = (($begin < $offset)
		     ? $begin - 1 : (abs($geneCoord[1]) - $offset - 1));
        print COORD_OUT "$id";
        if ($geneCoord[0] == 0) {push(@newCoord, 0);}
        else {push(@newCoord, abs($geneCoord[0]) - $modif);}
        for (my $i = 1; $i < $#geneCoord; $i++)
        {
            my $c = abs($geneCoord[$i]) - $modif;
            push(@newCoord, "$sens$c");
            print COORD_OUT " $sens$c";
        }
        if ($geneCoord[-1] == 0) {push(@newCoord, 0);}
        else {push(@newCoord, abs($geneCoord[-1]) - $modif);}
        print COORD_OUT "\n\n";
        print "done\n";
        print "\t-> Write the gff and info files..................";
        gffUtils::gffWrite($id, "$outputDir/Genes/$id$gff_ext",
			   "$outputDir/Genes/$id$info_ext",
			   "$newFastaFile", @newCoord);
        print "done\n";

        ## Extract Exons - UTR - Introns ##
        print "  - For build IMM matrices :\n";
        ExtractExon2($id, $newFastaFile, @newCoord);
        ExtractUTRs($id, $newFastaFile, $newCoord[0],
		    $newCoord[1], $newCoord[-2], $newCoord[-1]);
        ExtractIntr($id, $newFastaFile, $sens, @newCoord);

        ## Build WAM files ##
        if (defined($wam))
        {
            print "  - Build the WAM files....";
            WAM_FP_TP($id, $newFastaFile, @newCoord);
            print "..done\n";
        }
    }
}

#-----------------------------------------------------------------------------#
#-  Second approach.                                                         -#
#-----------------------------------------------------------------------------#
# About each step :
#   For each cDNA seq :
#     - Read Start & Stop
#     - Launch the EMBOSS extractseq on the cDNA file (to extract the CDS)
#     - Launch Sim4 against the extracted CDS (with options R=2 A=3 P=1 N=1
#       C=14) and generate an output file.
#     - Read in the Sim4 output file the gene coordinates (exons-introns).
#     - Test the Sim4 output file : exclusion if no exon found or if the CDS
#       is not completely aligned.
#     - Launch GeneSeqer against the extracted CDS and generate an output file.
#     - Read in the GeneSeqer output file the gene coordinates (exons-introns).
#     - Test the GeneSeqer output file : exclusion if no exon found or if the
#       CDS is not completely aligned.
#     - Choice the alignment Sim4 or GeneSeqer.
#     - Sim4 is use to find the UTR coordinates.
#     - Launch the EMBOSS extractseq on the genomic file (to extract the CDS +
#       a maximum offset on each side) and generate a new genomic file
#     - Add in the coord.txt file the gene coordinates (recomputed for the new
#       genomic sequence) and generate a GFF file (ID.fasta.gff) and a INFO
#       file  (ID.fasta.info).
#     - Extract Intron - Exon - UTR5/3 and add the sequence in each file :
#          -> Launch the EMBOSS extractseq on the cDNA file for exon and
#	      utr5/3
#          -> Launch the EMBOSS extractseq on the entire genomic file for
#             introns, and if the cDNA align by Sim4 is the reverse complement
#             launch the EMBOSS revseq on the extracted seq
#
# Informations about sequences (exclusion...) are in the $eugAdaptLog.
#-----------------------------------------------------------------------------#
sub GenoSim4
{
    while (my $cdnafile = <CDNA>)
    {
        print "------------------------------------";
	print "-------------------------------\n".($nbSeq + 1);
        flush STDOUT;
        print STDERR ".";

        my $genofile = <GENOMIC>;
	chomp $cdnafile;
        chomp $genofile;
	$nbSeq++;

	# cDNA Length ?
	my $lencDNA = "";
	if (!open(CDNAFILE, $cdnafile))
	  {
            print STDERR "ERROR: Could not open cDNA file $cdnafile\n";
            exit(2);
	  }
	while (my $line = <CDNAFILE>)
	  {
	    chomp $line;
	    if ($line =~ />/) { }
	    else              { $lencDNA .= $line; }
	  }
	$lencDNA = length($lencDNA);
	close CDNAFILE;

        ## ID extraction ##
        my $id = IDExtract($cdnafile);
        print " $id :\n";

        ## Read Start Stop file ##
        print"  - Read Start and Stop positions..............................";
        my $ssline = <SSCOORD>;
        my $start;
        my $stop;
	my $lenSS;       # Longueur Start->Stop (CDS)
        if ($ssline ne "\n")
	  {
	    chomp $ssline;
	    if ($ssline =~ /([0-9]+)\s+([0-9]+)/)
	      {
		$start = $1;
		$stop  = $2;
		$lenSS = $stop-$start+1;
	      }
	  }
	else
	  {
	    $start = 1;
	    $stop  = $lencDNA;
	    $lenSS = $stop-$start+1;
	  }
        print "done\n";

        ## Launch extractseq (EMBOSS) ##
	# start -> stop (cDNA)
        print"  - Launch extractseq (CDS from cDNA file).....................";
        ExtractSeq($cdnafile, $start, $stop, "s2e_mRNAStartStop");
        print "done\n";
	# utr5
        print"  - Launch extractseq (UTR5 from cDNA file)....................";
        ExtractSeq($cdnafile, 1, $start, "s2e_utr5");
        print "done\n";
	# utr3
        print"  - Launch extractseq (UTR3 from cDNA file)....................";
        ExtractSeq($cdnafile, $stop, 1000000, "s2e_utr3");
        print "done\n";

        my $simOut      = "$outputDir" . "/Genes/" . "$id" . "$sim4_out_ext";
	my @coordSim4   = ();
	my $sensSim4    = "";
        my $exclSim4    = 0;
	my $gseqerOut   = "$outputDir" . "/Genes/" . "$id" . "$gseqer_out_ext";
	my @coordGSeqer = ();
	my $sensGSeqer  = "";
        my $exclGSeqer  = 0;
        my @newCoord    = ();
	my $sens        = "";
        my $exclusion   = 0;
        my $simOutUTR5  = "$outputDir" . "/Genes/" . "$id" . "$sim4_out_ext2";
        my $simOutUTR3  = "$outputDir" . "/Genes/" . "$id" . "$sim4_out_ext3";
 	my $posUTR5     = 0;
	my $posUTR3     = 0;

        ## Launch sim4 ##
	# mRNA start stop contre genomique
	print"  - Launch sim4 (CDS against genomic file).....................";
	Sim4($simOut, $genofile, "s2e_mRNAStartStop");
	# utr5 contre genomique
	print"  - Launch sim4 (UTR5 against genomic file)....................";
	Sim4($simOutUTR5, $genofile, "s2e_utr5");
	# utr3 contre genomique
	print"  - Launch sim4 (UTR3 against genomic file)....................";
        Sim4($simOutUTR3, $genofile, "s2e_utr3");

        ## Sim4 output parsing ##
        Sim4Parsing($simOut, \$sensSim4, \@coordSim4);
	Sim4ParsingUTR($simOutUTR5, \$posUTR5, $id, "UTR5");
        Sim4ParsingUTR($simOutUTR3, \$posUTR3, $id, "UTR3");

	## Sim4 alignment test ##
	alignTest($lenSS, \$exclSim4, @coordSim4);
	##### End Sim4 #####

	## Launch GeneSeqer ##
	# mRNA start stop contre genomique
	print"  - Launch GeneSeqer (CDS against genomic file)................";
	GeneSeqer($gseqerOut, $genofile, "s2e_mRNAStartStop");
	
	## GeneSeqer output parsing ##
        GeneSeqerParsing($gseqerOut, \$sensGSeqer, \@coordGSeqer);
	
	## GeneSeqer alignment test ##
	alignTest($lenSS, \$exclGSeqer, @coordGSeqer);
	##### End GeneSeqer #####
	
	## GeneSeqer / Sim4 ?
	# Aucun alignement valide on exclue la seq
	if ($exclSim4 + $exclGSeqer == 2)
	  {
	    print EALOG $id . " :\t Sim4->";
	    if (@coordSim4 == 0) { print EALOG "NO EXON"; }
	    else { print EALOG "CDS isn't completely aligned";}
	    print EALOG "  GeneSeqer->";
	    if (@coordGSeqer == 0) { print EALOG "NO EXON\n"; }
	    else { print EALOG "CDS isn't completely aligned\n";}
	    $exclusion = 1;
	    $nbExclu++;
	  }
	# Seul GeneSeqer valide on utilise GeneSeqer
	elsif ($exclSim4   == 1)
	  {
	    @coordSim4 = @coordGSeqer;
	    $sens = $sensGSeqer;
	  }
	# Seul Sim4 valide on utilise Sim4
	elsif ($exclGSeqer == 1)
	  {
	    $sens = $sensSim4;
	  }
	# GeneSeqer et Sim4 valide mais differents
	elsif (join('', @coordGSeqer) ne join('', @coordSim4))
	  {
	    if (defined($auto))
	      {
		print EALOG $id . " :\t GeneSeqer / Sim4 ";
		print EALOG "(excluded by automatic process)\n";
	      }
	    else {
	      alignExpert($id, $simOut, $gseqerOut,
			  \@coordSim4,  \@coordGSeqer,
			  \$sensSim4,   \$sensGSeqer,  \$sens, \$exclusion);
	    }
	  }
	else { $sens = $sensSim4; }
	
	## Alignement sélectionné on continu...(si seq non exclue)
	if ($exclusion == 0)
	  {
	    for (my $i=0; $i<$#coordSim4; $i+=2)
	      {
		my $lenEx = $coordSim4[$i+1] - $coordSim4[$i] + 1;
		if ($lenEx > $exonMaxLen) { $exonMaxLen = $lenEx; }
		if ($lenEx < $exonMinLen) { $exonMinLen = $lenEx; }
		if ($lenEx < $seqInEALog)
		  {
		    print EALOG "$id :\t Contains exon < $seqInEALog";
		    print EALOG "        (seq not excluded !)\n";
		  }
	      }
	    if ($sens eq "-") {  $nbSeqRev++; }
	
	    ## Launch extractseq (EMBOSS) : + offset (DNA) ##
	    print "  - Launch extractseq (Create new genomic file,";
	    print " offset:$offset)...";
	    my $newFastaFile = "$outputDir/Genes/$id$newdna_ext";
	    my $begin = (($coordSim4[0] - $offset > 0)
			 ? $coordSim4[0] - $offset : 1);
	    my $end   = $coordSim4[$#coordSim4] + $offset;
	    ExtractSeq($genofile, $begin, $end, " $newFastaFile");
	    print "done\n";
	
	    ## Recomputing coord for the new genomic sequence ##
	    print "  - Recompute coordinates :\n";
	    print "\t-> Write in coord file...........................";
	    my $modif    = (($begin < $offset)
			    ? $begin - 1 : ($coordSim4[0] - $offset - 1));
	    my $newcoord = 0;
	
	    print COORD_OUT "$id";
	    if ($posUTR5 != 0) { $posUTR5 = $posUTR5 - $modif; }
	    if ($posUTR3 != 0) { $posUTR3 = $posUTR3 - $modif; }
	    if ($sens eq "")
	      {
		push(@newCoord, ($posUTR5 < 0) ? $sens."1" : "$sens$posUTR5");
	      }
	    else
	      {
		push(@newCoord, ($posUTR3 < 0) ? $sens."1" : "$sens$posUTR3");
	      }
	    foreach $newcoord (@coordSim4)
	      {
		my $c = $newcoord - $modif;
		push(@newCoord, "$sens$c");
		print COORD_OUT " $sens$c";
	      }
	    if ($sens eq "") { push(@newCoord, "$sens$posUTR3"); }
	    else             { push(@newCoord, "$sens$posUTR5"); }
	    print COORD_OUT "\n\n";
	    print "done\n";
	    print "\t-> Write the gff and info files..................";
	    gffUtils::gffWrite($id, "$outputDir/Genes/$id$gff_ext",
			       "$outputDir/Genes/$id$info_ext",
			       "$newFastaFile", @newCoord);
	    print "done\n";
	
	    ## Extract Exons - UTR - Introns ##
	    print "  - For build IMM matrices :\n";
	    ExtractExon($id, $cdnafile, $start, $stop);
	    ExtractUTRs($id, $cdnafile, $start, $start, $stop, 10000000);
	    unshift(@coordSim4, 0);
	    push(@coordSim4, 0);
	    ExtractIntr($id, $genofile, $sens, @coordSim4);
	
	    ## Build WAM files ##
	    if (defined($wam))
	      {
		print "  - Build the WAM files....";
		WAM_FP_TP($id, $newFastaFile, @newCoord);
		print "..done\n";
	      }
	  }
        else
	  {
            print "\n  ----> EXCLUSION :  !! WARNING !!\n";
	  }
      }

    unlink("s2e_mRNAStartStop");
    unlink("s2e_utr5");
    unlink("s2e_utr3");
}

#-----------------------------------------------------------------------------#
#-  ID extraction from genomic (1) or cdna (2) file.                         -#
#-----------------------------------------------------------------------------#
sub IDExtract
{
  my $lfile = shift;
  # JER, 19 avril 2004  je fais ma bidouille pour les ID MENS sans
  # toucher a l'existant
  if ( $lfile =~ /Mt[CDS]/ )
    {	
      my($id) = $lfile =~ /(Mt[CDS].+_GC)/;
      return ($id) if ( $id ne "" );
    }
  # Phi, 19 mai 2004 bidouille pour le jeu des belges
    if ( $lfile =~ /AC.+_AC/ )
    {	
      my($id) = $lfile =~ /(AC.+_AC.+)\.tfa/;
      return ($id) if ( $id ne "" );
    }
  if ($lfile =~ /\/?(.+\/)*([^\.]+)/) { return $2; }
  else { die "ERROR: Incorrect file list (line:$lfile)\n" }
}

#-----------------------------------------------------------------------------#
#-  Launch extractseq. (1) & (2)                                             -#
#-----------------------------------------------------------------------------#
sub ExtractSeq
{
    my ($seq, $begin, $end, $out) = @_;
    my $regions = "-regions $begin-$end";
    my $trash   = " >& /dev/null";
    system "$cmd_extractseq $seq $regions -outseq $out $trash";
}

#-----------------------------------------------------------------------------#
#-  Launch sim4. (2)                                                         -#
#-----------------------------------------------------------------------------#
sub Sim4
{
    my ($lout, $lgeno, $lcdna) = @_;
    system "$cmd_sim4 $lcdna $lgeno R=2 A=3 P=1 N=1 C=14 > $lout";
    print "done\n";
}

#-----------------------------------------------------------------------------#
#-  Launch GeneSeqer. (2)                                                    -#
#-----------------------------------------------------------------------------#
sub GeneSeqer
{
    my ($lout, $lgeno, $lcdna) = @_;
    system"$cmd_gseqer -s Arabidopsis -l $lgeno -e $lcdna -o $lout>&/dev/null";
    print "done\n";
}

#-----------------------------------------------------------------------------#
#-  Parse the sim4 output for extract the coordinates. (2)                   -#
#-----------------------------------------------------------------------------#
sub Sim4Parsing
{
  my ($lsimOut, $refsens, $refcoord) = @_;
  open(SIM, $lsimOut) || die "ERROR: Could not open file $lsimOut\n";

  while (my $line = <SIM>)
    {
      if ($line =~ /\(complement\)/) { $$refsens = "-"; }

      if ($line =~ /\(([0-9]+)\-([0-9]+)\)/)
        {
	  push(@{$refcoord}, $1);
	  push(@{$refcoord}, $2);
	}
    }
  close SIM;
}

#-----------------------------------------------------------------------------#
#-  Parse the sim4 output for extract the coordinates. UTRs for .gff !! (2)  -#
#-----------------------------------------------------------------------------#
#-  Prise en compte de tous les UTR ne contenant pas d'intron qqsoit le %    -#
#-  d'alignement !                                                           -#
#-----------------------------------------------------------------------------#
sub Sim4ParsingUTR
{
  my ($lsimOut, $refposUTR, $lid, $lutr) = @_;
  open(SIM, $lsimOut) || die "ERROR: Could not open file $lsimOut\n";
  my $nb_match = 0;
  my $sensUTR = "+";
  my $pos;

  while (my $line = <SIM>)
    {
      if ($line =~ /\(complement\)/) { $sensUTR = "-"; }

      if ($line =~ /\(([0-9]+)\-([0-9]+)\)/)
        {
	  if ( $line =~ /==/ ) # JER
	    {
	      print EALOG $lid . " :\t $lutr result not reliable".
		" (seq not excluded !)\n";
	      $nb_match = -1;
	      last;
	    }
	  if ($sensUTR eq "+")
	    {
	      if ($lutr eq "UTR5") { $pos = $1;	}
	      else                 { $pos = $2;	}
	    }
	  else
	    {
	      if ($lutr eq "UTR5") { $pos = $2; }
	      else    	           { $pos = $1;	}
	    }	
	  $nb_match++;
	}
    }
  if ( $nb_match != -1)
    {
      if ( $nb_match == 0 ) # JER
	{
	  print EALOG $lid . " :\t $lutr NO MATCH            ";
	  print EALOG "(seq not excluded !)\n";
	}
      elsif ( $nb_match != 1 )
	{
	  print EALOG $lid . " :\t $lutr with intron ?       ";
	  print EALOG "(seq not excluded !)\n";
	}
      else
	{
	  $$refposUTR = $pos;
	}
    }
  close SIM;
}

#-----------------------------------------------------------------------------#
#-  Parse the GeneSeqer output for extract the coordinates OF THE BEST       -#
#-  (%couverture*score moyen) ALIGNMENT. (2)                                 -#
#-----------------------------------------------------------------------------#
sub GeneSeqerParsing
  {
  my ($lgseqerOut, $refsens, $refcoord) = @_;
  open(GSEQER, $lgseqerOut) || die "ERROR: Could not open file $lgseqerOut\n";

  my $i = -1;
  my $j;
  my @coord;
  my @score;      # score
  my @lenAl;      # couverture
  my @nbVal;
  my $lencDNA = "";

  while (my $line = <GSEQER>)
    {
      if ($lencDNA eq ""  &&  $line =~ /EST sequence/)
	{
	  $line = <GSEQER>; $line = <GSEQER>;
	  while (1)
	    {
	      if ($line =~ /\s+[0-9]+\s+([ATCG ]+)/)
		{
		  $lencDNA .= $1;
		}
	      else
		{
		  $lencDNA =~ s/\s//g;
		  $lencDNA = length($lencDNA);
		  last;
		}
	      $line = <GSEQER>;
	    }
	}
      if ($line =~ /Exon/  &&  $line =~ /cDNA/)
	{
	  my @tab = split (/\s+/, $line);
	  if ($tab[2] == 1)
	    {
	      $i++;
	      $j = 0;
	      $lenAl[$i] = 0;
	      $score[$i] = 0;
	    }
	  $coord[$i][$j++] =  $tab[3];
	  $coord[$i][$j++] =  $tab[4];
	
	  if ($line =~ /([0-9]+) n\); score:/)   { $lenAl[$i] += $1; }
	  if ($line =~ /score: ([0-9]+.[0-9]+)/) { $score[$i] += $1; }
	  $nbVal[$i] =  $j;
	}
    }

  # On recherche l'indice correspondant au meilleur
  # alignement : $c_s = %couverture*score moyen.
  my $index = 0;
  my $c_s   = 0;

  for (my $k=0; $k<=$#score; $k++)
    {
      my $c = $lenAl[$k] / ($lencDNA * 100);
      my $s = $score[$k] / ($nbVal[$k]/2);
      if (($c * $s) > $c_s)
	{
	  $c_s   = $c * $s;
	  $index = $k;
	}
    }
  if (@score != 0)      # Si il y a au moins 1 match
    {
      if ($coord[$index][0] > $coord[$index][1])
	{
	  $$refsens = "-";
	  for (my $k=$nbVal[$index]-1; $k>=0; $k--)
	    {
	      push(@{$refcoord}, $coord[$index][$k]);
	    }
	}
      else
	{
	  for (my $k=0; $k<$nbVal[$index]; $k++)
	    {
	      push(@{$refcoord}, $coord[$index][$k]);
	    }
	}
    }
  close GSEQER;
}

#-----------------------------------------------------------------------------#
#-  Test alignment (2).                                                      -#
#-----------------------------------------------------------------------------#
#-  On exclue les séquences :                                                -#
#-      - sans match                                                         -#
#-      - dont la longueur [Start-Stop] != à la longueur de l'alignement     -#
#-----------------------------------------------------------------------------#
sub alignTest
{
  my ($llenSS, $refexcl, @lcoord) = @_;
  my $lenAllEx = 0;
  my $nbExon   = 0;

  for (my $i=0; $i<$#lcoord; $i+=2)
    {
      $lenAllEx += $lcoord[$i+1] - $lcoord[$i] + 1;
      $nbExon++;
    }
  if ($nbExon == 0  ||  $lenAllEx != $llenSS) { $$refexcl = 1; }
}

#-----------------------------------------------------------------------------#
#-  Expert alignment : choice Sim4, GeneSeqer, exclusion or auto process (2) -#
#-----------------------------------------------------------------------------#
sub alignExpert
{
  my ($lid, $lsimOut, $lgseqerOut, $refCdSim, $refCdGS,
      $refsensSim, $refsensGS, $refsens, $refexcl) = @_;
  my $choice = "";
  $$refsens = $$refsensSim;

  print "  - Sim4 / GeneSeqer give different alignments :\n";
  print "\t-> Sim4 :\n";
  for (my $i=0; $i<$#{$refCdSim}; $i+=2)
    {
      my $c = @{$refCdSim}[$i];
      print "\t   ";
      system "grep \"%\" $lsimOut | grep $c";
    }

  print "\t-> GeneSeqer :\n";
  for (my $i=0; $i<$#{$refCdGS}; $i+=2)
    {
      my $c = @{$refCdGS}[$i];
      print "\t  ";
      system "grep \"cDNA\" $lgseqerOut | grep $c";
    }
  print "\n\tFor more details about these alignments see the";
  print "\n\tsim4  ($lsimOut)  and  the";
  print "\n\tGeneSeqer ($lgseqerOut) files.\n";
  while (1)
    {
      print "\n\tEnter  your  choice :  Sim4 (s),  GeneSeqer (g),";
      print "\n\texclusion (e) or choice an automatic process (a)";
      print ".......";

      $choice = <STDIN>;
      chomp $choice;
      if ($choice eq 'e' || $choice eq 'a')
	{
	  print EALOG $lid . " :\t GeneSeqer / Sim4 (excluded by user)\n";
	  $nbExclu++;
	  $$refexcl = 1;
	  if ($choice eq 'a') { $auto = 'a'; }
	  last;
	}
      elsif ($choice eq 's') { last; }
      elsif ($choice eq 'g')
	{
	  @{$refCdSim} = @{$refCdGS};
	  $$refsens    = $$refsensGS;
	  last;
	}
    }
}

#-----------------------------------------------------------------------------#
#-  Extract exon (2).                                                        -#
#-----------------------------------------------------------------------------#
sub ExtractExon
{
    my ($lid, $lcdnafile, $lstart, $lstop) = @_;
    print "\t-> Launch extractseq (exon from cDNA file).......";
    print EXON "$lid ";
    ExtractSeq($lcdnafile, $lstart, $lstop, "s2e_t.txt");
    WriteExtract("s2e_t.txt", *EXON);
    print "done\n";
}

#-----------------------------------------------------------------------------#
#-  Extract exon. (1)                                                        -#
#-----------------------------------------------------------------------------#
sub ExtractExon2
{
    my ($lid, $lfastafile, @lcoord) = @_;

    print "\t-> Launch extractseq (exon from cDNA file)....";
    print EXON "$lid ";
    my $regions = "-regions ";
    my $trash   = " >& /dev/null";
    my $b       = "";
    my $e       = "";
    for (my $i = 1; $i < $#lcoord; $i += 2)
    {
        $b = abs($lcoord[$i]);
        $e = abs($lcoord[$i + 1]);
        $regions .= "$b-$e,";
        my $lenEx = $e - $b + 1;
        if ($lenEx > $exonMaxLen) {$exonMaxLen = $lenEx;}
        if ($lenEx < $exonMinLen) {$exonMinLen = $lenEx;}
    }

    system "$cmd_extractseq $lfastafile $regions -outseq s2e_t.txt $trash";

    if ($lcoord[1] < 0)
    {
        print "R";
        system "$cmd_revseq s2e_t.txt s2e_t2.txt $trash";
        WriteExtract("s2e_t2.txt", *EXON);
        print ".";
    }
    else
    {
        WriteExtract("s2e_t.txt", *EXON);
        print "..";
    }
    print ".done\n";
}

#-----------------------------------------------------------------------------#
#-  Extract utr5/3. (1) & (2)                                                -#
#-----------------------------------------------------------------------------#
sub ExtractUTRs
{
    my ($lid, $lcdnafile, $l5start, $l5stop, $l3start, $l3stop) = @_;

    ## UTR5 ##
    print "\t-> Launch extractseq (UTR5 from cDNA file)....";
    if (abs($l5start) - 1 >= 1)
    {
        print UTR5 "$lid ";
        if (abs($l5start) == abs($l5stop)) {$l5start = 1;}
        ExtractSeq($lcdnafile, abs($l5start), abs($l5stop) - 1, "s2e_t.txt");
        if ($l5stop < 0)
        {
            print "R";
            system "$cmd_revseq s2e_t.txt s2e_t2.txt >& /dev/null";
            WriteExtract("s2e_t2.txt", *UTR5);
            print "..";
        }
        else
        {
            WriteExtract("s2e_t.txt", *UTR5);
            print "...";
        }
        print "done\n";
    }
    else
    {
        print "...no utr5\n";
        $nbNoUTR5++;
    }
    ## UTR3 ##
    print "\t-> Launch extractseq (UTR3 from cDNA file)....";
    if (abs($l3stop) != 0)
    {
        ExtractSeq($lcdnafile, abs($l3start) + 1, abs($l3stop), "s2e_t.txt");
        if ($l3start < 0)
        {
            print "R";
            system "$cmd_revseq s2e_t.txt s2e_t2.txt >& /dev/null";
        }
        else {print ".";}
        open(TMP, "s2e_t.txt") || die "ERROR: Could not open file s2e_t.txt\n";
        my @linesTMP = <TMP>;
        close TMP;
        shift(@linesTMP);
        my $line = "$linesTMP[0]";
        if (length($line) > 2)
        {
            print UTR3 "$lid ";
            foreach $line (@linesTMP) {chomp $line; print UTR3 $line;}
            print UTR3 "\n";
            print "..done\n";
        }
        else
        {
            print "..no utr3\n";
            $nbNoUTR3++;
        }
    }
    else
    {
        print "...no utr3\n";
        $nbNoUTR3++;
    }
}

#-----------------------------------------------------------------------------#
#-  Extract intron. (1) & (2)                                                -#
#-----------------------------------------------------------------------------#
sub ExtractIntr
{
    my ($lid, $lgenofile, $lsens, @lcoord) = @_;

    if ($#lcoord != 3)
    {
        print "\t-> Launch extractseq (intron from genomic file) :\n\t   - ";
        for (my $i = 2; $i < $#lcoord - 1; $i += 2)
        {
            print($i/ 2);
            print INTRON "$lid.[" . ($i / 2) . "] ";
            ExtractSeq($lgenofile, abs($lcoord[$i]) + 1,
		       abs($lcoord[$i + 1]) - 1, "s2e_t.txt");
            my $lenInt = abs($lcoord[$i + 1]) - 1 - abs($lcoord[$i]);
            if ($lenInt > $intronMaxLen) {$intronMaxLen = $lenInt;}
            if ($lenInt < $intronMinLen) {$intronMinLen = $lenInt;}
            if ($lsens eq "-")
            {
                print "R";
                system "$cmd_revseq s2e_t.txt s2e_t2.txt >& /dev/null";
                WriteExtract("s2e_t2.txt", *INTRON);
                print ".";
            }
            else
            {
                WriteExtract("s2e_t.txt", *INTRON);
                print ".";
            }
        }
        print ".done\n";
    }
    else
    {
        print "\t-> Launch extractseq (intron from genomic file)..";
        print "no intron\n";
        $nbNoIntron++;
    }
}

#-----------------------------------------------------------------------------#
#-  Write the extracted sequences (exon, intron, utr5/3) in files. (1) & (2) -#
#-----------------------------------------------------------------------------#
sub WriteExtract
{
    my ($fileIn, $fileOut) = @_;
    open(TMP, $fileIn) || die "ERROR: Could not open file $fileIn\n";
    my @linesTMP = <TMP>;
    close TMP;
    my $line = "";
    shift(@linesTMP);
    foreach $line (@linesTMP) {chomp $line; print $fileOut $line;}
    print $fileOut "\n";
}

#-----------------------------------------------------------------------------#
#-  Extract and write signals with context for WAMBuilder. (1) & (2)         -#
#-----------------------------------------------------------------------------#
sub WAM_FP_TP
{
    my ($lid, $fasta, @coord) = @_;
    my $cmp = 1;
    my $pos;
    my $seq   = "";
    my $start = $coord[1];
    my $stop  = $coord[-2];
    my @don   = ();
    my @acc   = ();

    for (my $i = 2; $i < $#coord - 1; $i += 2)
    {
        push(@don, $coord[$i]);
        push(@acc, $coord[$i + 1]);
    }

    if ($coord[1] < 0)
    {
        system "$cmd_revseq $fasta fastaREV >& /dev/null";
        $fasta = "fastaREV";
    }

    open(FASTA, $fasta) || die "ERROR: Could not open file $fasta\n";
    my @fasta = <FASTA>;
    close FASTA;
    shift(@fasta);
    foreach my $line (@fasta) {chomp $line; $seq .= $line;}
    my $seqLen = length($seq);
    $seq =~ tr/atcg/ATCG/;

    if ($coord[1] < 0)
    {
        my @coordREV = ();
        foreach my $c (@coord) {unshift(@coordREV, $seqLen - abs($c) + 1);}
        $start = $coordREV[1];
        $stop  = $coordREV[-2];
        @don   = ();
        @acc   = ();
        for (my $i = 2; $i < $#coordREV - 1; $i += 2)
        {
            push(@don, $coordREV[$i]);
            push(@acc, $coordREV[$i + 1]);
        }
    }

    print VERIF "$lid\t";

    ## START ##
    my $v      = "ATG:NO";
    my $before = $wamstaB + $worder;
    print "Start...";
    while ($seq =~ /(?=([ATCG]{$before}(ATG)[ATCG]{$wamstaA}))/g)
    {
        $pos = pos($seq) + $before + 1;
        if ($pos != $start)
        {
            print FPSTART ">$lid " . $cmp++ . " ATG:$pos\n$1\n";
        }
        else
        {
            print TPSTART ">$lid " . $cmp++ . " ATG:$pos\n$1\n";
            $v = "ATG:OK";
        }
    }
    print VERIF "$v ";

    ## DONOR ##
    $v = 0;
    if (@don == 0) {push(@don, -1);}
    $cmp = 1;
    my $j     = 0;
    my $donor = $don[$j];
    $before = $wamdonB + $worder;
    print "Donor...";

    while ($seq =~ /(?=([ATCG]{$before}(GT)[ATCG]{$wamdonA}))/g)
    {
        $pos = pos($seq) + $before;
        if ($j < $#don && $pos > $donor) {$donor = $don[++$j];}
        if ($pos != $donor)
        {
            print FPDON ">$lid " . $cmp++ . " GT:$pos\n$1\n";
        }
        else
        {
            print TPDON ">$lid " . $cmp++ . " GT:$pos\n$1\n";
            $v++;
        }
    }
    if ($don[0] != -1)
    {
        print VERIF " GT$v/" . ($#don + 1);
        (($v / ($#don + 1) == 1) ? print VERIF ":OK " : print VERIF ":NO ");
    }

    ## ACCEPTOR ##
    $v = 0;
    if (@acc == 0) {push(@acc, -1);}
    $cmp = 1;
    $j   = 0;
    my $acceptor = $acc[$j];
    $before = $wamaccB + $worder;
    print "Acceptor...";

    while ($seq =~ /(?=([ATCG]{$before}(AG)[ATCG]{$wamaccA}))/g)
    {
        $pos = pos($seq) + $before + 1 + 2;
        if ($j < $#acc && $pos > $acceptor) {$acceptor = $acc[++$j];}
        if ($pos != $acceptor)
        {
            print FPACC ">$lid " . $cmp++ . " AG:$pos\n$1\n";
        }
        else
        {
            print TPACC ">$lid " . $cmp++ . " AG:$pos\n$1\n";
            $v++;
        }
    }
    if ($acc[0] != -1)
    {
        print VERIF " AG$v/" . ($#acc + 1);
        (($v / ($#acc + 1) == 1) ? print VERIF ":OK " : print VERIF ":NO ");
    }

    ## STOP ##
    $v      = "Stop:NO";
    $cmp    = 1;
    $before = $wamstoB + $worder;
    print "Stop...";
    my @st = ("TAA", "TAG", "TGA");
    for (my $i = 0; $i < 3; $i++)
    {

        while ($seq =~ /(?=([ATCG]{$before}($st[$i])[ATCG]{$wamstoA}))/g)
        {
            $pos = pos($seq) + $before + 1 + 2;
            if ($pos != $stop)
            {
                print FPSTOP ">$lid " . $cmp++ . " $st[$i]:$pos\n$1\n";
            }
            else
            {
                print TPSTOP ">$lid " . $cmp++ . " $st[$i]:$pos\n$1\n";
                $v = "$st[$i]:OK";
            }
        }
    }
    print VERIF " $v\n";
}

#-----------------------------------------------------------------------------#
#-  Print (stdout && EALOG) some informations about sequences. (1) & (2)     -#
#-----------------------------------------------------------------------------#
sub PrintInfo
{
    print "----------------------------------"
      . "---------------------------------\n";
    if ($nbExclu != 0)
    {
        print "! WARNING $nbExclu / $nbSeq sequence(s) have been excluded."
          . "\n!         You MUST read the $eugAdaptLog file.\n";
        print "---------------------------------"
	  . "----------------------------------\n";
    }
    print "Informations (for ".($nbSeq-$nbExclu)."/$nbSeq sequence(s)) :\n";
    print "   - Number of :\n";
    print "\t-> Sequences forward : " . ($nbSeq-$nbSeqRev-$nbExclu) . "\n";
    print "\t-> Sequences reverse : $nbSeqRev\n";
    print "\t-> Sequences without utr5   : $nbNoUTR5\n";
    print "\t-> Sequences without utr3   : $nbNoUTR3\n";
    print "\t-> Sequences without intron : $nbNoIntron\n";
    print "   - Length of :\n";
    print "\t-> Exon min.   : $exonMinLen\n";
    print "\t-> Exon max.   : $exonMaxLen\n";
    print "\t-> Intron min. : $intronMinLen\n";
    print "\t-> Intron max. : $intronMaxLen\n";
    system('echo "Finished on "`date`');
    print "-----------------------------------"
      . "--------------------------------\n";
    flush STDOUT;
    print STDERR "\n";

    print EALOG "----------------------------------"
      . "---------------------------------\n";
    print EALOG "Informations (for ".($nbSeq-$nbExclu)."/$nbSeq sequence(s))";
    print EALOG " :\n   - Number of :\n";
    print EALOG "\t-> Sequences forward : ".($nbSeq-$nbSeqRev-$nbExclu)."\n";
    print EALOG "\t-> Sequences reverse : $nbSeqRev\n";
    print EALOG "\t-> Sequences without utr5   : $nbNoUTR5\n";
    print EALOG "\t-> Sequences without utr3   : $nbNoUTR3\n";
    print EALOG "\t-> Sequences without intron : $nbNoIntron\n";
    print EALOG "   - Length of :\n";
    print EALOG "\t-> Exon min.   : $exonMinLen\n";
    print EALOG "\t-> Exon max.   : $exonMaxLen\n";
    print EALOG "\t-> Intron min. : $intronMinLen\n";
    print EALOG "\t-> Intron max. : $intronMaxLen\n";
    print EALOG "-----------------------------------"
      . "--------------------------------\n";
}
