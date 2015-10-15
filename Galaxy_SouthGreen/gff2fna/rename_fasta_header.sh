#!/bin/bash

#qsub -q bioinfo.q -b yes -V -N format1 rename_fasta_header.sh

# ARATH
perl -ne 'if($_ =~/^>Chr(\d.+) /){$label=$1;printf (">ARATH%02d\n",$label)} else{print $_}' /bank/arabidopsis_thaliana/TAIR10_whole_chromosomes.fasta > /bank/genfam/ARATH/ARATH-TAIR10-chromosome-genfam.fna
grep '>' /bank/genfam/ARATH/ARATH-TAIR10-chromosome-genfam.fna
perl -ne 'if($_ =~/^Chr(\d+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/){$label=$1;printf ("ARATH%02d\tTAIR10$2\n",$label)}' /bank/genfam/ARATH/TAIR10_GFF3_genes_transposons.gff > /bank/genfam/ARATH/ARATH-TAIR10-sequence_feature-genfam.gff3
awk '{if($0 !~ /^##/){print $1}}' /bank/genfam/ARATH/ARATH-TAIR10-sequence_feature-genfam.gff3 | sort | uniq

# BRADI (gff3 chr + scaff)
#perl -ne 'if($_ =~/^>Bd(\d+)/){$label=$1;printf (">BRADI%02d\n",$label)} elsif ($_ =~/^>scaffold_(\d.+)/) {$label=$1;printf (">BRADI_scaffold%03d\n",$label)} else{print $_}' /bank/brachypodium_distachyon/Bdistachyon_192.fa > /bank/genfam/BRADI/BRADI-JGI1-chromosome-genfam.fna
#grep '>' /bank/genfam/BRADI/BRADI-JGI1-chromosome-genfam.fna
#perl -ne 'if($_ =~/^Bd(.+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/){$label=$1;printf ("BRADI%02d\tJGI1$2\n",$label)} elsif ($_ =~/^scaffold_(.+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/) {$label=$1;printf ("BRADI_scaffold%03d\tJGI1$2\n",$label)}else{print $_}' /bank/brachypodium_distachyon/Bdistachyon_192_gene_exons.gff3 > /bank/genfam/BRADI/BRADI-JGI1-sequence_feature-genfam.gff3
#awk '{if($0 !~ /^##/){print $1}}' /bank/genfam/BRADI/BRADI-JGI1-sequence_feature-genfam.gff3 | sort | uniq

#GLYMA (présence de scdcaffold 4 digits dans génome)
#perl -ne 'if($_ =~/^>Gm(\d.+)/){$label=$1;printf (">GLYMA%02d\n",$label)} elsif ($_ =~/^>scaffold_(\d.+)/) {$label=$1;printf (">GLYMA_scaffold%04d\n",$label)} else{print $_}' /bank/glycine_max/Gmax_v1.1_189.fa > /bank/genfam/GLYMA/GLYMA-JGI1-chromosome-genfam.fna
#grep '>' /bank/genfam/GLYMA/GLYMA-JGI1-chromosome-genfam.fna
#perl -ne 'if($_ =~/^Gm(.+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/){$label=$1;printf ("GLYMA%02d\tJGI1$2\n",$label)} elsif ($_ =~/^scaffold_(.+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/) {$label=$1;printf ("GLYMA_scaffold%04d\tJGI1$2\n",$label)}else{print $_}' /bank/genfam/GLYMA/GLYMA.gff3 > /bank/genfam/GLYMA/GLYMA-JGI1-sequence_feature-genfam.gff3
#awk '{if($0 !~ /^##/){print $1}}' /bank/genfam/GLYMA/GLYMA-JGI1-sequence_feature-genfam.gff3 | sort | uniq

#GOSRA (présence de scaffold 4 digits dans génome)
#perl -ne 'if($_ =~/^>Chr(\d.+)/){$label=$1;printf (">GOSRA%02d\n",$label)} elsif ($_ =~/^>scaffold_(\d.+)/) {$label=$1;printf (">GOSRA_scaffold%04d\n",$label)} else{print $_}' /bank/gossypium_raimondii/JGI_pseudochromosome > /bank/genfam/GOSRA/GOSRA-JGI1-chromosome-genfam.fna
#grep '>' /bank/genfam/GOSRA/GOSRA-JGI1-chromosome-genfam.fna
#perl -ne 'if($_ =~/^Chr(.+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/){$label=$1;printf ("GOSRA%02d\tJGI1$2\n",$label)} elsif ($_ =~/^scaffold_(.+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/) {$label=$1;printf ("GOSRA_scaffold%04d\tJGI1$2\n",$label)}else{print $_}' /bank/genfam/GOSRA/GOSRA.gff3 > /bank/genfam/GOSRA/GOSRA-JGI1-sequence_feature-genfam.gff3
#awk '{if($0 !~ /^##/){print $1}}' /bank/genfam/GOSRA/GOSRA-JGI1-sequence_feature-genfam.gff3 | sort | uniq

# LOTJA (gff3 hétérogène avec fin avec n°accession et mettre sur 2 digit et gff3 commence a 0)
#perl -ne 'if($_ =~/^>LjChr(\d+)/){$label=$1;printf (">LOTJA%02d\n",$label)} elsif ($_ =~/^>.+_(\d.+)/) {$label=$1;printf (">LOTJA_contig%06d\n",$label)}else{print $_}' /bank/genfam/LOTJA/LOTJA-Kazusa2.5-chromosome-contig-genfam.fna > /bank/genfam/LOTJA/LOTJA-Kazusa2.5-chromosome-contig2-genfam.fna
#grep '>' /bank/genfam/LOTJA/LOTJA-Kazusa2.5-chromosome-contig2-genfam.fna
#perl -ne 'if($_ =~/^Chr(.+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/){$label=$1+1;printf ("LOTJA%02d\tKazusa2.5$2\n",$label)} elsif ($_ =~/^.+_(.+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/) {$label=$1;printf ("LOTJA_contig%06d\tKazusa2.5$2\n",$label)}else{print $_}' /bank/genfam/LOTJA/Lj2.5_gene_models.gff3 > /bank/genfam/LOTJA/LOTJA-Kazusa2.5-sequence_feature-genfam.gff3
#awk '{if($0 !~ /^##/){print $1}}' /bank/genfam/LOTJA/LOTJA-Kazusa2.5-sequence_feature-genfam.gff3 | sort | uniq

# MEDTR (gff3 hétérogène avec n°accession dans colone ID: AC130806 en début et chromosome à la suite)(-)
#perl -ne 'if($_ =~/^>chr(\d+)/){$label=$1;printf (">MEDTR%02d\n",$label)} elsif ($_ =~/^>AC(\d+)(\.\d+)|^>CU(\d+)(\.\d+)/) {$label=$1;printf (">MEDTR_scaffold%06d$2\n",$label)} else{print $_}' /bank/medicago_truncatula/Mtruncatula_198.fa > /bank/genfam/MEDTR/MEDTR-Mt3.5v5-chromosome-genfam.fna
#grep '>' /bank/genfam/MEDTR/MEDTR-Mt3.5v5-chromosome-genfam.fna
#perl -ne 'if($_ =~/^chr(.+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/){$label=$1;printf ("MEDTR%02d\tMt3.5v5$2\n",$label)} elsif ($_ =~/^AC(\d+)(\.\d+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/) {$label=$1;printf ("MEDTR_scaffold%06d$2\tMt3.5v5$3\n",$label)}elsif ($_ =~/^CU(\d+)(\.\d+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/) {$label=$1;printf ("MEDTR_scaffold%06d$2\tMt3.5v5$3\n",$label)}else{print $_}' /bank/medicago_truncatula/Mtruncatula_198_gene_exons.gff3 > /bank/genfam/MEDTR/MEDTR-Mt3.5v5-sequence_feature-genfam.gff3
#awk '{if($0 !~ /^##/){print $1}}' /bank/genfam/MEDTR/MEDTR-Mt3.5v5-sequence_feature-genfam.gff3 | sort | uniq

# MUSAC
#perl -ne 'if($_ =~/^>chr(.+)/){$label=$1;printf (">MUSAC%02d\n",$label)} else{print $_}' /bank/musa_acuminata/MUSA_pseudochromosome.fna > /bank/genfam/MUSAC/MUSAC-Musa1-chromosome-genfam.fna
#grep '>' /bank/genfam/MUSAC/MUSAC-Musa1-chromosome-genfam.fna
#perl -ne 'if($_ =~/^chr(.+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/){$label=$1;printf ("MUSAC%02d\tMusa1$2\n",$label)} else{print $_}' /bank/genfam/MUSAC/MUSAC.gff3 > /bank/genfam/MUSAC/MUSAC-Musa1-sequence_feature-genfam.gff3
#awk '{if($0 !~ /^##/){print $1}}' /bank/genfam/MUSAC/MUSAC-Musa1-sequence_feature-genfam.gff3 | sort | uniq

# ORYSI
#perl -ne 'if($_ =~/^>Osi(\d.+)/){$label=$1;printf (">ORYSI%02d\n",$label)} else{print $_}' /bank/oryza_sativa_indica/BGI2_pseudochromosome > /bank/genfam/ORYSI/ORYSI-BGI2-chromosome-genfam.fna
#grep '>' /bank/genfam/ORYSI/ORYSI-BGI2-chromosome-genfam.fna
#perl -ne 'if($_ =~/^Osi(.+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/){$label=$1;printf ("ORYSI%02d\tBGI2$2\n",$label)} else{print $_}' /bank/gff3/ORYSI1.gff3 > /bank/genfam/ORYSI/ORYSI-BGI2-sequence_feature-genfam.gff3
#perl -ne 'if($_ !~ /^##/){s/Osi/ORYSI/; print $_}' /bank/gff3/ORYSI1.gff3 > /bank/genfam/ORYSI/ORYSI-BGI2-sequence_feature-genfam.gff3
#awk '{if($0 !~ /^##/){print $1}}' /bank/genfam/ORYSI/ORYSI-BGI2-sequence_feature-genfam.gff3 | sort | uniq

# ORYSJ
#perl -ne 'if($_ =~/^>Os(\d.+)/){$label=$1;printf (">ORYSJ%02d\n",$label)} else{print $_}' /bank/oryza_sativa_japonica/MSU7_pseudochromosome.fna > /bank/genfam/ORYSJ/ORYSJ-MSU7-chromosome-genfam.fna
#grep '>' /bank/genfam/ORYSJ/ORYSJ-MSU7-chromosome-genfam.fna
#perl -ne 'if($_ =~/^Os(.+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/){$label=$1;printf ("ORYSJ%02d\tMSU7$2\n",$label)} else{print $_}' /bank/gff3/ORYSJ1.gff3 > /bank/genfam/ORYSJ/ORYSJ-MSU7-sequence_feature-genfam.gff3
#perl -ne 'if($_ !~ /^##/){s/Os/ORYSJ/; print $_}' /bank/gff3/ORYSJ1.gff3 > /bank/genfam/ORYSJ/ORYSJ-MSU7-sequence_feature-genfam.gff3
#awk '{if($0 !~ /^##/){print $1}}' /bank/genfam/ORYSJ/ORYSJ-MSU7-sequence_feature-genfam.gff3 | sort | uniq

# POPTR (gff3: chr + scaff) (-)
#perl -ne 'if($_ =~/^>Chr(\d+)/){$label=$1;printf (">POPTR%02d\n",$label)} elsif ($_ =~/^>scaffold_(\d.+)/) {$label=$1;printf (">POPTR_scaffold%04d\n",$label)} else{print $_}' /bank/populus_trichocarpa/Ptrichocarpa_210.fa > /bank/genfam/POPTR/POPTR-JGI2-chromosome-genfam.fna
#grep '>' /bank/genfam/POPTR/POPTR-JGI2-chromosome-genfam.fna
#perl -ne 'if($_ =~/^Chr(\d+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/){$label=$1;printf ("POPTR%02d\tJGI2$2\n",$label)} elsif ($_ =~/^scaffold_(\d+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/) {$label=$1;printf ("POPTR_scaffold%04d\tJGI2$2\n",$label)} else{print $_}' /bank/populus_trichocarpa/Ptrichocarpa_210_gene_exons.gff3 > /bank/genfam/POPTR/POPTR-JGI2-sequence_feature-genfam.gff3
#awk '{if($0 !~ /^##/){print $1}}' /bank/genfam/POPTR/POPTR-JGI2-sequence_feature-genfam.gff3 | sort | uniq

# SOLLC
#perl -ne 'if($_ =~/^>SL2.40ch(\d.+)/){$label=$1;printf (">SOLLC%02d\n",$label)} else{print $_}' /bank/solanum_lycopersicum/S_lycopersicum_chromosomes.2.40.fa > /bank/genfam/SOLLC/SOLLC-ITAG2.40-chromosome-genfam.fna
#grep '>' /bank/genfam/SOLLC/SOLLC-ITAG2.40-chromosome-genfam.fna
#perl -ne 'if($_ =~/^SL2.40ch(.+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/){$label=$1;printf ("SOLLC%02d\tITAG2.40$2\n",$label)} else{print $_}' /bank/genfam/SOLLC/SOLLC.gff3 > /bank/genfam/SOLLC/SOLLC-ITAG2.40-sequence_feature-genfam.gff3
#perl -ne 'if($_ !~ /^##/){s/SL2.40ch/SOLLC/; print $_}' /bank/genfam/SOLLC/SOLLC.gff3 > /bank/genfam/SOLLC/SOLLC-ITAG2.40-sequence_feature-genfam.gff3
#awk '{if($0 !~ /^##/){print $1}}' /bank/genfam/SOLLC/SOLLC-ITAG2.40-sequence_feature-genfam.gff3 | sort | uniq

# SORBI (tail du gff3 différent: super_) (-)
#perl -ne 'if($_ =~/^>Sb(\d.+)/){$label=$1;printf (">SORBI%02d\n",$label)} else{print $_}' /bank/genfam/SORBI/Sbicolor_79.fa > /bank/genfam/SORBI/SORBI-JGI1.4-chromosome-genfam.fna
#grep '>' /bank/genfam/SORBI/SORBI-JGI1.4-chromosome-genfam.fna
#perl -ne 'if($_ =~/^chromosome_(.+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/){$label=$1;printf ("SORBI%02d\tJGI1.4$2\n",$label)} elsif ($_ =~/^super_(.+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/) {$label=$1;printf ("SORBI_scaffold%05d\tJGI1.4$2\n",$label)}else{print $_}' /bank/genfam/SORBI/Sbicolor_79_gene_exons.gff3 > /bank/genfam/SORBI/SORBI-JGI1.4-sequence_feature-genfam.gff3
#awk '{if($0 !~ /^##/){print $1}}' /bank/genfam/SORBI/SORBI-JGI1.4-sequence_feature-genfam.gff3 | sort | uniq

# THECC (Un_random => 00)
#perl -ne 'if($_ =~/^>Tc(\d.+)/){$label=$1;printf (">THECC%02d\n",$label)} else{print $_}' /bank/theobroma_cacao/TC_pseudochromosome_tot.fna > /bank/genfam/THECC/THECC-COCOA1-chromosome-genfam.fna
#grep '>' /bank/genfam/THECC/THECC-COCOA1-chromosome-genfam.fna
#perl -ne 'if($_ =~/^Tc(.+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/){$label=$1;printf ("THECC%02d\tCOCOA1$2\n",$label)} else{print $_}' /bank/genfam/THECC/THECC.gff3 > /bank/genfam/THECC/THECC-COCOA1-sequence_feature-genfam.gff3
#perl -ne 'if($_ !~ /^##/){s/Tc/THECC/; print $_}' /bank/genfam/THECC/THECC.gff3 > /bank/genfam/THECC/THECC-COCOA1-sequence_feature-genfam.gff3
#awk '{if($0 !~ /^##/){print $1}}' /bank/genfam/THECC/THECC-COCOA1-sequence_feature-genfam.gff3 | sort | uniq

# VITVI
#perl -ne 'if($_ =~/^>chr(.+)/){$label=$1;printf (">VITVI%02d\n",$label)} else{print $_}' /bank/vitis_vinifera/GENOSCOPE1_pseudochromosome > /bank/genfam/VITVI/VITVI-GENOSCOPE1-chromosome-genfam.fna
#grep '>' /bank/genfam/VITVI/VITVI-GENOSCOPE1-chromosome-genfam.fna
#perl -ne 'if($_ =~/^chr(.+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/){$label=$1;printf ("VITVI%02d\tGENOSCOPE1$2\n",$label)} else{print $_}' /bank/genfam/VITVI/VITVI.gff3 > /bank/genfam/VITVI/VITVI-GENOSCOPE1-sequence_feature-genfam.gff3
#awk '{if($0 !~ /^##/){print $1}}' /bank/genfam/VITVI/VITVI-GENOSCOPE1-sequence_feature-genfam.gff3 | sort | uniq

# MAIZE
#perl -ne 'if($_ =~/^>(\d+)/){$label=$1;printf (">MAIZE%02d\n",$label)} elsif ($_ =~/^>(UNKNOWN)/){$label=$1;printf ("MAIZE%02d\n",$label)}else{print $_}' /bank/zea_maize/Zmays_181.fa > /bank/genfam/MAIZE/MAIZE-MGDB5b60-chromosome-genfam.fna
#grep '>' /bank/genfam/MAIZE/MAIZE-MGDB5b60-chromosome-genfam.fna
#perl -ne 'if($_ =~/^(\d+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/){$label=$1;printf ("MAIZE%02d\tMGDB5b60$2\n",$label)} elsif ($_ =~/^(UNKNOWN)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/){$label=$1;printf ("MAIZE%02d\tMGDB5b60$2\n",$label)}else{print $_}' /bank/zea_maize/Zmays_181_gene_exons.gff3 > /bank/genfam/MAIZE/MAIZE-MGDB5b60-sequence_feature-genfam.gff3
#awk '{if($0 !~ /^##/){print $1}}' /bank/genfam/MAIZE/MAIZE-MGDB5b60-sequence_feature-genfam.gff3 | sort | uniq

# MALDO
#perl -ne 'if($_ =~/^>MDC(\d+\.\d+)/){$label=$1;printf (">MALDO_scaffold$label\n")} else{print $_}' /bank/genfam/MALDO/Mdomestica-chr.fna > /bank/genfam/MALDO/MALDO-JGI1-chromosome-genfam.fna
#grep '>' /bank/genfam/MALDO/MALDO-JGI1-chromosome-genfam.fna
#perl -ne 'if($_ =~/^MDC(\d+\.\d+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/){$label=$1;printf ("MALDO_scaffold$label\tJGI1$2\n")} else{print $_}' /bank/genfam/MALDO/MALDO.gff3 > /bank/genfam/MALDO/MALDO-JGI1-sequence_feature-genfam.gff3
#awk '{if($0 !~ /^##/){print $1}}' /bank/genfam/MALDO/MALDO-JGI1-sequence_feature-genfam.gff3 | sort | uniq

# MANES
#perl -ne 'if($_ =~/^>scaffold(\d.+)/){$label=$1;printf (">MANES_scaffold%05d\n",$label)} else{print $_}' /bank/genfam/MANES/Mesculenta-chr.fna > /bank/genfam/MANES/MANES-JGI4.1-chromosome-genfam.fna
#grep '>' /bank/genfam/MANES/MANES-JGI4.1-chromosome-genfam.fna
#perl -ne 'if($_ =~/^scaffold(.+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/){$label=$1;printf ("MANES_scaffold%02d\tJGI4.1$2\n",$label)} else{print $_}' /bank/genfam/MANES/MANES.gff3 > /bank/genfam/MANES/MANES-JGI4.1-sequence_feature-genfam.gff3
#perl -ne 'if($_ !~ /^##/){s/scaffold/MANES_scaffold/; print $_}' /bank/genfam/MANES/MANES.gff3 > /bank/genfam/MANES/MANES-JGI4.1-sequence_feature-genfam.gff3
#awk '{if($0 !~ /^##/){print $1}}' /bank/genfam/MANES/MANES-JGI4.1-sequence_feature-genfam.gff3 | sort | uniq

# RICCO
#perl -ne 'if($_ =~/^>(\d.+)/){$label=$1;printf (">RICCO_scaffold%05d\n",$label)} else{print $_}' /bank/genfam/RICCO/Rcommunis-chr.fna > /bank/genfam/RICCO/RICCO-JGI0.1-chromosome-genfam.fna
#grep '>' /bank/genfam/RICCO/RICCO-JGI0.1-chromosome-genfam.fna
#perl -ne 'if($_ =~/^(.+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/){$label=$1;printf ("RICCO_scaffold%05d\tJGI0.1$2\n",$label)} else{print $_}' /bank/genfam/RICCO/RICCO.gff3 > /bank/genfam/RICCO/RICCO-JGI0.1-sequence_feature-genfam.gff3
#awk '{if($0 !~ /^##/){print $1}}' /bank/genfam/RICCO/RICCO-JGI0.1-sequence_feature-genfam.gff3 | sort | uniq

# SETIT
#perl -ne 'if($_ =~/^>scaffold_(\d+)/){$label=$1;printf (">SETIT_scaffold%03d\n",$label)} else{print $_}' /bank/genfam/SETIT/Sitalica_164_v2.fa > /bank/genfam/SETIT/SETIT-JGI1-chromosome-genfam.fna
#grep '>' /bank/genfam/SETIT/SETIT-JGI1-chromosome-genfam.fna
#perl -ne 'if($_ =~/^scaffold_(.+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/){$label=$1;printf ("SETIT_scaffold%03d\tJGI1$2\n",$label)} else{print $_}' /bank/genfam/SETIT/SETIT.gff3 > /bank/genfam/SETIT/SETIT-JGI1-sequence_feature-genfam.gff3
#awk '{if($0 !~ /^##/){print $1}}' /bank/genfam/SETIT/SETIT-JGI1-sequence_feature-genfam.gff3 | sort | uniq

# SOLTU
#perl -ne 'if($_ =~/^>chr(\d.+)/){$label=$1;printf (">SOLTU%02d\n",$label)} else{print $_}' /bank/genfam/SOLTU/Stuberosum_206_v3.fna > /bank/genfam/SOLTU/SOLTU-JGI3-chromosome-genfam.fna
#grep '>' /bank/genfam/SOLTU/SOLTU-JGI3-chromosome-genfam.fna
#perl -ne 'if($_ =~/^chr(.+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/){$label=$1;printf ("SOLTU%02d\tJGI3$2\n",$label)} else{print $_}' /bank/genfam/SOLTU/Stuberosum_206_v3.4.gff3 > /bank/genfam/SOLTU/SOLTU-JGI3-sequence_feature-genfam.gff3
#perl -ne 'if($_ !~ /^##/){s/chr/SOLTU/; print $_}' /bank/genfam/SOLTU/Stuberosum_206_v3.4.gff3 > /bank/genfam/SOLTU/SOLTU-JGI3-sequence_feature-genfam.gff3
#awk '{if($0 !~ /^##/){print $1}}' /bank/genfam/SOLTU/SOLTU-JGI3-sequence_feature-genfam.gff3 | sort | uniq

# COFCA
#perl -ne 'if($_ =~/^>chr(.+)/){$label=$1;printf (">COFCA%02d\n",$label)} else{print $_}' /bank/genfam/COFCA/pseudomolecules.fa > /bank/genfam/COFCA/COFCA-GENOSCOPE1-chromosome-genfam.fna
#grep '>' /bank/genfam/COFCA/COFCA-GENOSCOPE1-chromosome-genfam.fna
#perl -ne 'if($_ =~/^chr(.+)\t.+\ttransposable_element(\t.+\t.+\t.+\t.+\t.+\t.+)/){$label=$1;printf ("COFCA%02d\tGENOSCOPE1\ttransp_element$2\n",$label)} elsif ($_ =~/^chr(.+)\t.+(\t.+\t.+\t.+\t.+\t.+\t.+\t.+)/){$label=$1;printf ("COFCA%02d\tGENOSCOPE1$2\n",$label)} else{print $_}' /bank/genfam/COFCA/COFCA.gff3 > /bank/genfam/COFCA/COFCA-GENOSCOPE1-sequence_feature-genfam.gff3
#awk '{if($0 !~ /^##/){print $1}}' /bank/genfam/COFCA/COFCA-GENOSCOPE1-sequence_feature-genfam.gff3 | sort | uniq