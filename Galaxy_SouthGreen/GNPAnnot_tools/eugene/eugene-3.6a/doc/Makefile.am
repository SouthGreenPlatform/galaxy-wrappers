all : EuGeneDoc.pdf

EuGeneDoc.pdf : GenerateDoc.tcl Doc.tex ParaOptimization.fig ParaOptimization.pdf Sequences/SYNO_ARATH.fasta  Sequences/SYNO_ARATH.fasta.blast0 Sequences/SYNO_ARATH.fasta.blast1 Sequences/SYNO_ARATH.fasta.blast2 Sequences/SYNO_ARATH.fasta.est Sequences/SYNO_ARATH.fasta.spliceP Sequences/SYNO_ARATH.fasta.splicePR Sequences/SYNO_ARATH.fasta.splices Sequences/SYNO_ARATH.fasta.splicesR Sequences/SYNO_ARATH.fasta.starts  Sequences/SYNO_ARATH.fasta.startsR
	tclsh GenerateDoc.tcl

cleandoc :
	rm -f CmdFlags.* EuGeneDoc.* *~

