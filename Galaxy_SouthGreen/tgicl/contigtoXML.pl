#!/usr/bin/perl

use Storable;
use Bio::SeqIO;
use strict;

my $working_dir = $ARGV[0];
open Sav,">$working_dir/redundant.fasta";
open Sav_consensus, ">$working_dir/consensus.fasta";


my @rep = ('asm_1',  'asm_2',  'asm_3',  'asm_4',  'asm_5',  'asm_6','asm_7',  'asm_8');
for (@rep) {
	my $rep = $_;
        if ( -e "$working_dir/$rep/ACE" ) {
        	rename( "$working_dir/$rep/ACE", "$working_dir/$rep/ACE.old" )
                          or die "pas possible de renommer\n";
                        open G, "$working_dir/$rep/ACE.old"
                          or die "pas possible de renommer\n";
                        open Sav2, ">$working_dir/$rep/ACE"
                          or die "pas possible de renommer\n";
                        while (<G>) {
                                if (/CO (CL(\d+)Contig\d+) (\d+) (\d+) .+\n/) {
                                        print Sav2 "\n\n$_";
                                }
                                else {
                                        print Sav2 $_;
                                }
                        }
                        close Sav2;
                        close G;
                        open F, "$working_dir/$rep/ACE";

                        my $alignment_temp;
                        my %clust;
			my $ligne    = 0;
                        my $nbreCont = 0;
                        $/ = "\n\n\n\n";

                        while (<F>) {
                                my ( $contigName, $cl, $length, $nbreRead, $seq ) =
                                  /CO (CL(\d+)Contig\d+) (\d+) (\d+) .+\n/;

                                my $titre      = "AS 1 $nbreRead\n";
                                my $paragraphe = $_;
                                if ( $paragraphe =~ /AS \d+ \d+\n/ ) {
                                        $paragraphe =~ s /AS \d+ \d+\n//g;
                                }
                                else {
                                        $titre .= "\n";
                                }
                                $paragraphe = $titre . $paragraphe;
                                if ($nbreRead) {
					
                                        my $in = Bio::SeqIO->new(
                                                -file     => "$working_dir/$rep/contigs",
                                                '-format' => 'fasta'
                                        );
                                        my ($sequence);
                                        while ( my $seq = $in->next_seq() ) {
                                                if ( $seq->id eq "$contigName" ) {
                                                        $sequence = $seq->seq;
                                                }
                                        }
                                        print Sav_consensus ">$contigName\n$sequence\n";
                                        my %hash_cons;
                                        $hash_cons{'length'} = $length;
                                        $hash_cons{'seq'}    = $sequence;
                                        $hash_cons{'cl'}     = $cl;
                                        $hash_cons{'ace'}    = $paragraphe;
                                        my $compt = 0;

                                        while ( $paragraphe =~ /AF (.+?) [CU] .+\n/g ) {
                                                my $est = $1;
                                                $compt++;
						$hash_cons{'sousSeq'}{$compt} = $est;


                                                print Sav ">$est\n" . "todo" . "\n";
                                      

                                        }

                                        my $xmlfilecons = $working_dir . "/$contigName.xml";
                                        store \%hash_cons, $xmlfilecons;
                                }
                        }
                        $/ = "\n";
                        close F;
                }
        }
        close Sav;
        close Sav_consensus;
