package FIT; #First In-data Treatment

use strict;
use Bio::SeqIO;
use moduls::DMC;

use Exporter;

our @ISA= qw( Exporter );

# these CAN be exported.
our @EXPORT_OK = qw( readTypeFile createAllelesLocus identifyAlleles createAllelesDist);

# these are exported by default.
our @EXPORT = qw( readTypeFile createAllelesLocus identifyAlleles createAllelesDist);

=for comment

	FILE		FIT.pm - Include functions for load a type in mlst test (1st step)
	AUTHOR		Verdier Axel
	CREATE_DATE	04-18-2014
	LAST_DATE	04-25-2014
	FUNCTIONS	readTypeFile createAllelesLocus identifyAlleles createAllelesDist
			readTypeFile - Load a fasta file, and return an array of the locus (name and sequence)
			createAllelesLocus - create the alleles of each Locus matrix and return it in a table
			identifyAlleles - identify the allele of each locus of the new type and and this allele in the array of locus alleles if necessary (and this distance with the other in the array of alleles dist per locus)
			createAllelesDist - read distances matrix files and create the big alleles distances matrix
=cut


#Functions
#====================
=head1 readTypeFile

	Overview	Load a fasta file, and return an array of the locus (name and sequence)
	
	parameters	0:inputPath 1:nbLocus
			inputPath : the path to the fasta file that contains the locus of the type to test
			nbLocus : the number of locus for this disease files
			
	return		An array (for all locus) of hash table (name and sequence).
	
	example		@test=FIT::readTypeFile(path/file.fas); #get the array from the fasta file
			print "$test[0]{name}\n";#display the name of the first locus.
=cut
sub readTypeFile{
	my $inputPath=\@_[0];
	my $nbLocus=\@_[1];
	my @sequenceIn=();
	if (! -e $$inputPath){
		print "$$inputPath n'existe pas ou n'est pas un fichier.\n";
		exit;
	}
	else {
		my $numLocus=0;
		my $seqio = Bio::SeqIO->new(-file => $$inputPath, '-format' => 'Fasta');
		while(my $seq = $seqio->next_seq){
			if($numLocus>$$nbLocus-1){
				print "\ntoo much sequence in $$inputPath. Expected $$nbLocus, 1 per locus\n";
				exit;
			}
			$sequenceIn[$numLocus]={"name" => $seq->display_id, "sequence" => $seq -> seq};
	 		$numLocus++;
		}
		if ($numLocus<$$nbLocus-1){
			print "not enough sequence in $$inputPath. Expected $$nbLocus, 1 per locus\n";
		}
	}
	return @sequenceIn;
}

#====================
=head1 createAllelesLocus

	Overview	create the alleles of each Locus matrix and return it in a table
	
	parameters	0:pathLocusFiles 1:nbLocus
			pathLocusFiles: the path to the 'model' fasta file that contains the locus alleles => for example/file1.fa example/file2.fa [...], it's example/file.fa
			nbLocus : the number of locus for this disease files
			
	return		An array (locus) of an array (alleles) of hash tablea (id, name, sequence) 
	
	example		my @allelesLocus=FIT::createAllelesLocus(path/file.fa, NBLOCUS);
			print $allelesLocus[0][0]{sequence}."\n";
=cut
sub createAllelesLocus{
	my $pathLocusFiles=\@_[0];
	my $nbLocus=\@_[1];
	
	my @allelesLocus;
	
	#Load alleles for each locus
	for (my $n=1; $n <= $$nbLocus; $n++){
		my ($extension)=$$pathLocusFiles=~m/^.*(\.\w+)$/;
		my ($path)=$$pathLocusFiles=~m/^(.*)$extension$/;
		$path.=$n.$extension;
		
		my $seqio = Bio::SeqIO->new(-file => $path, -format => 'Fasta');
		my $id=0;
		my @allelesInLocus;
		while (my $seq = $seqio->next_seq){
			@allelesInLocus[$id]={'id'=>$id, 'name' => $seq -> display_id, 'sequence' => $seq -> seq};
			if ($allelesInLocus[$id]{sequence}=~m/\r/){print "\\r found !! Be carefull !!! \n\n";}
			$id++;
		}
		$allelesLocus[$n-1]=\@allelesInLocus;
	}
	return @allelesLocus
}

#====================
=head1 identifyAlleles

	Overview	identify the allele of each locus of the new type and and this allele in the array of locus alleles if necessary (and this distance with the other in the array of alleles dist per locus)
	
	parameters	0:typeIn 1:allelesLocus 2:allelesDist 3:DISTANCE_VALUE 4:nbLocus 5:minGapSize
			typeIn : an array of sequence (1 for each locus)
			allelesLocus : the array of locus allelesInLocus
			allelesDist: the array of alleles distance per locus
			DISTANCE_VALUE: a hash table of the value distance for align nucleotides type (match, mismatch, indel, gap, vntr)
			nbLocus : the number of locus for this disease
			minGapSize : the minimal size of a gap
			
			
	return		N/A
	
	example		my @inType=FIT::readTypeFile($inputPath, NBLOCUS);
			my @allelesLocus=FIT::createAllelesLocus($path, NBLOCUS);
			my @allelesDist=FIT::createAllelesDist($path, $NBLOCUS);
			FIT::identifyAlleles(\@inType,\@allelesLocus, \@allelesDist, $NBLOCUS);
=cut
sub identifyAlleles{
	my $typeIn=\@{@_[0]};
	my $allelesLocus=\@{@_[1]};
	my $allelesDist=\@{@_[2]};
	my $DISTANCE_VALUE=\%{@_[3]};
	my $nbLocus=\@_[4];
	my $minGapSize=\@_[5];
		
	#Assign its allele for each locus of the typeIn
	my $numLocus=0;
	while($numLocus<$$nbLocus){
		my $nameType=$$typeIn[$numLocus]{name};
		my $sequenceType=$$typeIn[$numLocus]{sequence};
		my $idAllele=0;
		my $testLocusNb=$numLocus+1;
		
		
		my $alleleFind=0;
		my $test=@$allelesLocus[$numLocus];
		foreach (@{@$allelesLocus[$numLocus]}){
			my $allele=$$allelesLocus[$numLocus][$idAllele]{sequence};
			if(simplyAlign($sequenceType, $allele)){
				$alleleFind=1;
				last;
				}
			$idAllele++;
		}
		if($alleleFind==1){$$typeIn[$numLocus]={'id'=>$idAllele, 'name' => $nameType, 'sequence' => $sequenceType};}
		else{
			$$allelesLocus[$numLocus][$idAllele]={'id'=>$idAllele, 'name' => $nameType, 'sequence' => $sequenceType};
			$$typeIn[$numLocus]={'id'=>$idAllele, 'name' => $nameType, 'sequence' => $sequenceType};
			
			#calc the distance between this allele and the other and add them in the matrix
			for (my $i=0; $i<=$idAllele;$i++){
				my %align=DMC::align2Seq($$allelesLocus[$numLocus][$idAllele]{sequence}, $$allelesLocus[$numLocus][$i]{sequence}, $$minGapSize);
				$$allelesDist[$numLocus][$idAllele][$i]=$align{mismatches}*$$DISTANCE_VALUE{mismatch} + $align{indels}*$$DISTANCE_VALUE{indel} + $align{gaps}*$$DISTANCE_VALUE{gap} + $align{vntr}*$$DISTANCE_VALUE{vntr} + $align{matches}*$$DISTANCE_VALUE{match};
				#Check if there is an error in the in-type sequences
				if ($$allelesDist[$numLocus][$idAllele][$i]>=50){die("ERROR: distance >=50 : improbable sequence for locus number ".($idAllele+1)."\n");}
			}
		}

		
		$numLocus++;
	}
}

#====================
=head1 createAllelesDist

	Overview	read distances matrix files and create the big alleles distances matrix
		
	parameters	0:pathFiles 1:nbLocus 2:distAllelesMatrix 3:MINGAPSIZE 4:DISTANCE_VALUE
			pathFiles: the path to the matrix files + the base of files names (test/mccp-1-matrix.csv, test/mccp-2-matrix.csv, ... => test/mccp.csv)
			
	return		an array of alleles distance matrix per locus @distAllelesMatrix[$locus][$i][$j]
	
	example		-
=cut
sub createAllelesDist{
	my $pathFiles=@_[0];
	my $nbLocus=@_[1];
	my $MINGAPSIZE=@_[2];
	my %DISTANCE_VALUE=%{@_[3]};
	my @distAllelesMatrix;
	
	
	for (my $loc=0; $loc<$nbLocus; $loc ++){
		my ($file, $ext)=$pathFiles=~m/^(.*)\.(\w+)$/;
		my @alleles;
		my $numSeq=0;
		my $nLoc=$loc+1;
		my $inputPath=$file.$nLoc."\.".$ext;
		my $seqio = Bio::SeqIO->new(-file => $inputPath, '-format' => 'Fasta');
		while(my $seq = $seqio->next_seq){
			$alleles[$numSeq]{name}=$seq->display_id;
			$alleles[$numSeq]{sequence}=$seq->seq;
	 		$numSeq++;
		}
		$numSeq--; #Last incrementation correction
		#Do the matrix
		$|=1;
		my @matrix;
		for(my $i=0; $i<=$numSeq; $i++){
			for(my $j=0; $j<=$i; $j++){
				my %align=DMC::align2Seq($alleles[$i]{sequence}, $alleles[$j]{sequence}, $MINGAPSIZE);
				$matrix[$i][$j]=$align{mismatches}*$DISTANCE_VALUE{mismatch} + $align{indels}*$DISTANCE_VALUE{indel} + $align{gaps}*$DISTANCE_VALUE{gap} + $align{vntr}*$DISTANCE_VALUE{vntr} + $align{matches}*$DISTANCE_VALUE{match};
			}
		}
		$distAllelesMatrix[$loc]=\@matrix;
	}
	return @distAllelesMatrix;
}

#====================
=head1 simplyAlign

	Overview	compare 2 sequence and return a boolean : true if they're the same
		
	parameters	0:seq1 1:seq2
			seq1: the first sequence
			seq2: the seconde sequence
			
	return		an boolean value of the correspondance
	
	example		simplyAlign('ATGCGT','ATGGGT');
=cut

sub simplyAlign{
	return @_[0]=~m/^(@_[1])$/;
}
1;
