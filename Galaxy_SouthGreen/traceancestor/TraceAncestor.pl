#!/usr/bin/perl

# Version 1 : 01/06/18 

use warnings;
use strict;

#Command line processing.
use Getopt::Long;

# Initialization of files and variables that can be chosen by the user to launch the script.
my $input;
my $vcf;
my $color;
my $ploidy;
my $window_size = 10;
my $window_cut = 100;
my $LOD = 3;
my $freq = 0.99;
my $help = "";
my $hybrid;
my $indivs;
my $merge;
Getopt::Long::Configure ('bundling');
GetOptions ('t|input=s' => \$input, ## reference matrix
	    'v|vcf=s' => \$vcf, ## vcf of hybrids. Hybrids should not have space in their names and no special characters.
	    'c|color=s'=> \$color, ## color file
	    'p|ploidy=i' => \$ploidy, ## ploidy of the hybrid population (2, 3 or 4)
	    'w|window=i' =>\$window_size, ## number of ancestral markers by windows on a chromosome
	    'l|lod=i' =>\$LOD, ## LOD value to conclude for one hypothesis or an other
	    's|threshold=f' =>\$freq, ## Theoretical frequency used to calcul the LOD
	    'k|cut=i' => \$window_cut, ## number of K bases in one window
	    'i|ind=s' => \$hybrid, ## particular hybrid you want to focus on (optional). If -f not indicated.
	    'f|focus=s'=> \$indivs, ## file with several hybrids to focus on (optional). If -i not indicated.
	    'm|merge!'=> \$merge,	
	    'h|help!'=> \$help,)
or die("Error in command line arguments\n");
if($help) { ## display this help when -help is called
	print "\nusage: TraceAncestor.pl [-t matrix file] [-v vcf file] [-c color file] [-p ploidy] [-w number of markers by window] [-s threshold for LOD] [-k window size in K-bases] [-i hybrid name to focus on] [-f file of hybrids to focus on]\n\n";
	print "-t | --input : reference matrice. \n";
	print "-v | --vcf : vcf of the hybrid population\n";
	print "-c | --color : color file\n";
	print "-p | --ploidy : ploidy of the hybrid population\n";
	print "-w | --window : number of markers by window\n";
	print "-l | --lod : LOD value to conclude for one hypothesis\n";
	print "-s | --freq : theoretical frequency used to calcul the LOD\n";
	print "-k | --cut : number of K bases in one window\n";
	print "-i | --ind : particular hybrid you want to focus on. If -f not indicated (optional).\n";
	print "-f | --focus : file containing several hybrids to focus on (more than one). Each one asociated with a unique number.\n";
	print "-h | --help : display this help\n\n";
	exit;
}

## check to presence of necessary parameters
if (!$input){
	print "You need a matrix file (-t)\n";
	exit;
}   
if (!$vcf){
	print "You need a vcf (-v)\n";
	exit;
} 
if (!$color){
	print "You need a color file (-c)\n";
	exit;
}  
if (!$ploidy){
	print "You need to indicate the ploidy number of your population (-p)\n";
	exit;
} 
if ($indivs && $hybrid){
	print "You can't have a focus file (-f) and focus on a single hybrid (-i) at the same time\n";
	exit;
}

#------------------------------------------------------------------------------------------------------------------------------------------------------
# STEP 0 : If we work on WGS --> we redo an ancestor matrix with only snp that are also in the vcf file.

if ($merge){
	my@val;
	my%mergeAncestor;
	open(FA,"$input") or die ("Error: matrix of ancestors wont open\n"); #opening of the ancestor matrix
	while (my $li = <FA>){ # for each marker
		chomp($li);
		if ($li !~ m/^ancestor/){
			my@ligne = split("\t", $li);
			my$ancetre = $ligne[0];
			my$chromosome = $ligne[1];
			my$position = $ligne[2];
			my$allele = $ligne[3];
			$mergeAncestor{$ancetre}{$chromosome}{$position} = $allele;
		}
	}
	close FA;
	open(FX, ">new_AncestorMatrix");
	print FX "ancestor\tchromosome\tposition\tallele\n";
	open(FV,"$vcf") or die ("Error: can't open the vcf file\n");
	while (my $li = <FV>){ # for each marker
		chomp($li);
		if ($li !~ m/^#.+$/){
			@val = split("\t",$li);
			my$chr = $val[0]; # chromosome
			my$pos = $val[1]; # position
			for my$anc (sort keys(%mergeAncestor)){
				for my$c(sort keys(%{$mergeAncestor{$anc}})){
					for my$p(sort {$a<=>$b} keys(%{$mergeAncestor{$anc}{$c}})){
						if(($c eq $chr) && ($p eq $pos)){
							print FX "$anc\t$c\t$p\t$mergeAncestor{$anc}{$c}{$p}\n";
						}
					}
				}
			}
		}
	}
	close FV;
	close FX;
}

#------------------------------------------------------------------------------------------------------------------------------------------------------

# STEP 1: creation of the not overlapping windows of $window_size (parameter chose by the user) markers

my%bloc_position; # key 1 = ancestor / key 2 = chromosome / key 3 = position of the end of the window => value = every marker positions in one window separated by a -
my%endChro; # key = chromosome => value = position of the last marker of every chromosome. We will assume this is the length of the chromosome.
my@ref; # list of the name of each ancestral genome.
my%allele_refalt; # key 1 = ancestor / key 2 = chromosome / key3 = position of the marker => value = base of the ancestral allele
my%catPos; # key 1 = ancestor / key 2 = chromosome => value = concatenation of every position of a windows separated by a -
my%hashRef; # marker counter by window
my@chromosomList; # list of chromosomes in the ancestor matrix
if(!$merge){
	open(F1,"$input") or die ("Error: matrix of ancestors wont open\n"); #opening of the ancestor matrix
}
else{
	open(F1,"new_AncestorMatrix") or die ("Error: new matrix of ancestors wont open\n"); #opening of the ancestor matrix
}
while (my $li = <F1>){ # for each marker
	chomp($li);
	# if line is not the header
	if ($li !~ m/^ancestor/){
		my@ligne = split("\t", $li);
		my$ancetre = $ligne[0];
		my$chromosome = $ligne[1];
		my$position = $ligne[2];
		my$allele = $ligne[3];
		# list of last position for every chromosome
		$endChro{$chromosome} = $position;
		if (!grep { $_ eq $ancetre } @ref){
			push (@ref, $ancetre);
		}
		# list of chromosomes
		if (!grep { $_ eq $chromosome } @chromosomList){
			push (@chromosomList, $chromosome);
		}
		# We make windows of 10 consecutives maarkers -> %bloc_position
		if (!defined($catPos{$ancetre}{$chromosome})){
			$catPos{$ancetre}{$chromosome} = $position; #initialization of concatenation of positions at the begining of each chromosome.
			$hashRef{$ancetre}{$chromosome} = 1; #initialization of the marker counter
		}
		if ($hashRef{$ancetre}{$chromosome} != $window_size-1){ # if the counter is less than the wanted number of markers by windows
			if ($catPos{$ancetre}{$chromosome} eq ""){
				$catPos{$ancetre}{$chromosome} = $position;
				$hashRef{$ancetre}{$chromosome} = 1;
			}
			elsif ($catPos{$ancetre}{$chromosome} ne $position){
				$catPos{$ancetre}{$chromosome} = $catPos{$ancetre}{$chromosome}."-".$position;
				$hashRef{$ancetre}{$chromosome}++;
			}
		}
		elsif($hashRef{$ancetre}{$chromosome} == $window_size-1){ # if the counter equal the wanted number of markers by windows
			$bloc_position{$ancetre}{$chromosome}{$position} = $catPos{$ancetre}{$chromosome}."-".$position;
			$hashRef{$ancetre}{$chromosome} = 0; # reinitialization
			$catPos{$ancetre}{$chromosome} = "";
		}
		$allele_refalt{$ancetre}{$chromosome}{$position} = $allele;
	}
}
close F1;
# We are numbering each windows in a hash %blocs.			
my %blocs;
my %blocspos;
my$end;
for my $parent (sort keys %bloc_position) {
    	for my $chro (sort keys %{$bloc_position{$parent}}) {
		my$boo=0; #if boo==0, it's the begining of a chromosome
		my$num = 0;
		for my $pos (sort {$a<=>$b} keys %{$bloc_position{$parent}{$chro}}) {
			$num++;
			$blocs{$parent}{$chro}{$num} = $bloc_position{$parent}{$chro}{$pos};
			my@marqueur = split("-",$blocs{$parent}{$chro}{$num});
			if ($boo == 0){
				$blocspos{$parent}{$chro}{$num}= "1-".$marqueur[scalar(@marqueur)-1];
			}
			else{
				$blocspos{$parent}{$chro}{$num}=$end."-".$marqueur[scalar(@marqueur)-1];
			}
			$end = $marqueur[scalar(@marqueur)-1];
			$boo = 1;
		}
	}
}

#--------------------------------------------------------------------------------------------------------

# STEP 2: Calculation of the frequency of specifics reads for each ancestors by window.

open(F3,"$vcf") or die ("Error: can't open the vcf file\n");
my@individus; # list of samples. Samples names begin at $individu[9]
my@valBlock; # On line of the vcf
my%nbAncestralOrNot; # key 1: ancestor / key 2 : hybrid / key 3 : chromosome / key 4: window number / key 5 : "NOT" or "ANCESTRAL" / value : If "NOT" = sums of non-ancestral reads. If "ANCESTRAL" = sum of ancestral reads
my$chr;
my$pos="";
my$numbloc; 
my$AlleleRef; # Ancestral allele for a given position and a given ancestor
my$Ancestral = "ANCESTRAL";
my$Not = "NOT";
my$base1; # first base for one position in the vcf
my$base2; # second base for one position in the vcf (can be multiple)
my@bases2; # list for tue multiple alternative base in the vcf alternative base
my%frequences; # Frequency of ancestral alleles by marker
my%NM; # hash counting the number of marker really present in the vcf for each window (could be usefull later).
my%realEndChro; #real length of each chromosome. We will use these values if they are available.
my$booEndChro = 0; # if = 0 we use %endChro else we use %realEndChro as end of chromosome
while (my$li = <F3>){
	chomp($li);
	if ($li =~ m/^#.+$/){
		if ($li =~ m/(#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT.+)$/){
			@individus = split(/\s/,$li); # List of samples
		}
		if ($li =~ m/##contig=<ID=(\S+)\,length=(\d+)\>/){
			$realEndChro{$1} = $2; # length ($2) for every chromosomes ($1)
			$booEndChro = 1;
		}
	}
	else{ # for each position / line in the vcf
		@valBlock = split("\t",$li);
		$chr = $valBlock[0]; # chromosome
		if(!$chr){
			next;
		}
		if(!grep { $_ eq $chr } @chromosomList){ # we work only on the chromosomes that are in the matrix of ancestors.
			next;
		};
		$pos = $valBlock[1]; # position
		$base1 = $valBlock[3];
		$base2 = $valBlock[4];
		my@base2;
		if ($base2 =~ /[\w|\*\,]+/){ # if base2 is multiple, we consider the first one (the more frequent)
			@bases2 = split(",", $base2);
		}
		for my $parent (sort keys %blocs){
			foreach my $num(sort {$a<=>$b} keys(%{$blocs{$parent}{$chr}})){
				my @allpos = split("-",$blocs{$parent}{$chr}{$num}); # list of the markers in this window
				if (grep { $_ == $pos } @allpos){ # We check that the snp is in the reference marker list.
					$NM{$parent}{$chr}{$num}++; # incrementation of the number of existing markers by windows.
					$AlleleRef = $allele_refalt{$parent}{$chr}{$pos};
					for (my$ind = 9; $ind < scalar(@individus); $ind++){
						if ($hybrid){ # if we have a focus on only one hybrid :
							if($individus[$ind] eq $hybrid){
								if ($valBlock[$ind] =~ m/^[\d+\.\/]+\:(\-\d+|\d+)\,([\d+\,\-]+)\:(\d+)\:.+/){
									my$refer = $1; # number of reads for the base1
									my@alters = split(",",$2);
									my$sumreads = $3; # number of reads for the for a marqueur (base1 + base2)
									if($sumreads > 0){ # if base1 + base2 != 0
										if ($AlleleRef eq $base1){ # if base1 is the ancestral allele
											$frequences{$parent}{$individus[$ind]}{$chr}{$num}{$pos} = $refer / $sumreads; # frequency of the ancestral allele for this marker
											if (defined($nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Ancestral})){
												$nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Ancestral} = $nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Ancestral} + $refer;
												$nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Not} = $nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Not} + ($sumreads-$refer);
											}
											else{
												$nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Ancestral} = $refer;
												$nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Not} = ($sumreads-$refer);
											}
										}
										else{
											for(my$a=0; $a < scalar(@bases2); $a++){
												if ($AlleleRef eq $bases2[$a]){ # if one of the base2 is the ancestral allele

													$frequences{$parent}{$individus[$ind]}{$chr}{$num}{$pos} = $alters[$a] / $sumreads; # frequency of the ancestral allele for this marker		
													if (defined($nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Ancestral})){
														$nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Ancestral} = $nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Ancestral} + $alters[$a];
														$nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Not} = $nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Not} + ($sumreads-$alters[$a]);
													}
													else{
														$nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Ancestral} = $alters[$a];
														$nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Not} = ($sumreads-$alters[$a]);
													}
												}
											}
										}
									}
								}
							}
						}
						else{ # if we have no focus on one hybrid but on all hybrid or on a list of hybrids
							if ($valBlock[$ind] =~ m/^[\d+\.\/]+\:(\-\d+|\d+)\,([\d+\,\-]+)\:(\d+)\:.+/){
								my$refer = $1; # number of reads for the base1
								my@alters = split(",",$2);
								my$sumreads = $3; # number of reads for the for a marqueur (base1 + base2)
								if($sumreads > 0){ # if base1 + base2 != 0 
									if ($AlleleRef eq $base1){ # if base1 if the ancestral allele
										$frequences{$parent}{$individus[$ind]}{$chr}{$num}{$pos} = $refer / $sumreads; # frequency of the ancestral allele for this marker
										if (defined($nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Ancestral})){
											$nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Ancestral} = $nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Ancestral} + $refer;
											$nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Not} = $nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Not} + ($sumreads-$refer);
										}
										else{
											$nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Ancestral} = $refer;
											$nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Not} = ($sumreads-$refer);
										}
									}
									else{
										for(my$a=0; $a < scalar(@bases2); $a++){
											if ($AlleleRef eq $bases2[$a]){ # if one of the base2 if the ancestral allele
												$frequences{$parent}{$individus[$ind]}{$chr}{$num}{$pos} = $alters[$a] / $sumreads; # frequency of the ancestral allele for this marker		
												if (defined($nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Ancestral})){
													$nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Ancestral} = $nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Ancestral} + $alters[$a];
													$nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Not} = $nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Not} + ($sumreads-$alters[$a]);
												}
												else{
													$nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Ancestral} = $alters[$a];
													$nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Not} = ($sumreads-$alters[$a]);
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
close F3;


my%percentage; # frequency of ancestral allele in a given window (ancestral reads / sum of ancestral and non ancestral reads). If there is no marker of a given widow in the vcf --> window is NA (will become NE later).
for my $parent(sort keys(%blocs)){
	if ($hybrid){ # if focus on one hybrid
		for my $chr(sort keys(%{$blocs{$parent}})){
			for my $num(sort {$a<=>$b} keys(%{$blocs{$parent}{$chr}})){
				if (!exists ($nbAncestralOrNot{$parent}{$hybrid}{$chr}{$num}{$Ancestral})){
					$percentage{$parent}{$hybrid}{$chr}{$num} = "NA";
				}
				else{
					my $sum = $nbAncestralOrNot{$parent}{$hybrid}{$chr}{$num}{$Ancestral} + $nbAncestralOrNot{$parent}{$hybrid}{$chr}{$num}{$Not};
					if($sum == 0){
						$percentage{$parent}{$hybrid}{$chr}{$num} = "NA";
					}
					else{
						$percentage{$parent}{$hybrid}{$chr}{$num} = $nbAncestralOrNot{$parent}{$hybrid}{$chr}{$num}{$Ancestral}/$sum;
					}
				}
			}
		}
		
	}
	else{ 
		for (my$ind = 9; $ind < scalar(@individus); $ind++){
			for my $chr(sort keys(%{$blocs{$parent}})){
				for my $num(sort {$a<=>$b} keys(%{$blocs{$parent}{$chr}})){					
					if (!exists ($nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Ancestral})){
						$percentage{$parent}{$individus[$ind]}{$chr}{$num} = "NA";
					}
					else{
						my $sum = $nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Ancestral} + $nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Not};
						if($sum == 0){
							$percentage{$parent}{$individus[$ind]}{$chr}{$num} = "NA";
						}
						else{
							$percentage{$parent}{$individus[$ind]}{$chr}{$num} = $nbAncestralOrNot{$parent}{$individus[$ind]}{$chr}{$num}{$Ancestral}/$sum;
						}
					}
				}
			}
		}
	}
}

open (Ffreq, ">ancestorFreq") or die ("cant write ancestorFreq\n"); # output => frequency of ancestral alleles along chromosomes
print Ffreq "Hybrid\tAncestry\tChromosom\tPosition_Start\tPosition_End\tFrequence\n"; # header of ancestorfreq
for my$parent (sort keys(%blocspos)){
	if ($hybrid){ # focus on one particular hybrid
		for my$chro (sort keys(%{$blocspos{$parent}})){
			for my$num (sort {$a<=>$b} keys(%{$blocspos{$parent}{$chro}})){
				my@position = split("-",$blocspos{$parent}{$chro}{$num});
				my$start = $position[0];
				# Correction of the begining of each bloc. -> +1 postition.
				my$start_cor = 1;
				if ($start != 1){
					$start_cor = $start + 1;
				}
				my$end = $position[1];
				print Ffreq "$hybrid\t$parent\t$chro\t$start_cor\t$end\t$percentage{$parent}{$hybrid}{$chro}{$num}\n";
			}
		}
	}
	elsif ($indivs){ # focus on all the hybrids in the focus file
		open(Findivs, "$indivs") or die ("can't open hybrids file"); # read the focus file to extract a list of hybrids to analyze
		my%indivs; #hash of focus file : key = individu name , value = number
		while (my $li = <Findivs>){
			chomp($li);
			if ($li =~ m/\d+\s+\S+/){
				my@i=split(/\s+/,$li);
				$indivs{$i[1]} = $i[0];
			}
		}
		close Findivs;
		for my$ind (sort keys(%indivs)){
			for my$chro (sort keys(%{$blocspos{$parent}})){
				for my$num (sort {$a<=>$b} keys(%{$blocspos{$parent}{$chro}})){
					my@position = split("-",$blocspos{$parent}{$chro}{$num});
					my$start = $position[0];
					# Correction of the begining of each bloc. -> +1 postition.
					my$start_cor = 1;
					if ($start != 1){
						$start_cor = $start + 1;
					}
					my$end = $position[1];
					print Ffreq "$ind\t$parent\t$chro\t$start_cor\t$end\t$percentage{$parent}{$ind}{$chro}{$num}\n";
				}
			}
		}
	}
	else{ # print all
		for (my$ind = 9; $ind < scalar(@individus); $ind++){
			for my$chro (sort keys(%{$blocspos{$parent}})){
				for my$num (sort {$a<=>$b} keys(%{$blocspos{$parent}{$chro}})){
					my@position = split("-",$blocspos{$parent}{$chro}{$num});
					my$start = $position[0];
					# Correction of the begining of each bloc. -> +1 postition.
					my$start_cor = 1;
					if ($start != 1){
						$start_cor = $start + 1;
					}
					my$end = $position[1];
					print Ffreq "$individus[$ind]\t$parent\t$chro\t$start_cor\t$end\t$percentage{$parent}{$individus[$ind]}{$chro}{$num}\n"
				}
			}
		}
	}
}
close Ffreq;

#--------------------------------------------------------------------------------------------------------------------
# STEP 3 : Estimate allelic dosage of each ancestors by windows. We use the Maximum of likelihood compared to theoretical frequency (LOD). 

my@dosage; # estimated frequencies.
if($ploidy == 2){
	$dosage[0] = 1-$freq;
	$dosage[1] = 0.50;
	$dosage[2] = $freq;
}
if($ploidy == 3){
	$dosage[0] = 1-$freq;
	$dosage[1] = 0.33;
	$dosage[2] = 0.66;
	$dosage[3] = $freq;
}
if($ploidy == 4){
	$dosage[0] = 1-$freq;
	$dosage[1] = 0.25;
	$dosage[2] = 0.50;
	$dosage[3] = 0.75;
	$dosage[4] = $freq;
}

# fonction of log10
sub log10 {
     my $n = shift;
     return log($n)/log(10);
}

# calculation of LOD
my%dosageAllelique; # key 1 = hybrid / key 2 = parent / key 3 = chromosome / key 4 = window number => value = allelic dosage of an ancestor for a window.
my$dos; 
my$ref; # Number of ancestral reads
my$alt; # Number of reads non ancestral 
for my $parent (sort keys %percentage) {
    for my $ind (sort keys %{$percentage{$parent}}) {
        for my $chr (sort keys %{$percentage{$parent}{$ind}}) {
     		for my $num (sort {$a<=>$b} keys %{$percentage{$parent}{$ind}{$chr}}) {
			if ($percentage{$parent}{$ind}{$chr}{$num} eq "NA"){
				$dosageAllelique{$ind}{$parent}{$chr}{$num} = "NE";
			}
			else{
				$ref = $nbAncestralOrNot{$parent}{$ind}{$chr}{$num}{$Ancestral}; # ancestral reads
				$alt = $nbAncestralOrNot{$parent}{$ind}{$chr}{$num}{$Not}; # non ancestral reads

				if($ploidy == 2){ # simplified version
					# hyp 0 vs hyp 1 :	
					my$lod0vs1 = $ref*log10($dosage[0]/$dosage[1]) + $alt*log10($dosage[2]/$dosage[1]);
					# If lod > 3 --> hyp 0 kept
					if ($lod0vs1 > $LOD){ # hyp 0 > hyp 1
						$dos = 0;
						$dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
					}
					elsif($lod0vs1 < -$LOD){ # hyp 1 > hyp 0
						# hyp 1 vs hyp 2 :
						my$lod1vs2 = $ref*log10($dosage[1]/$dosage[2]) + $alt*log10($dosage[1]/$dosage[0]);
						# If lod 3 --> hyp 1 kept
						if ($lod1vs2 > $LOD){ # hyp 1 > hyp 2
							$dos = 1;
							$dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
						}
						elsif($lod1vs2 < -$LOD){ # hyp 1 < hyp 2
							$dos = 2;
							$dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
						}
						else{ # indetermination between hyp1 and hyp2
							$dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; # Indetermination 
						}
					}
					else{ # indetermination between hyp1 and hyp0
						$dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; # Indetermination 
					}
				}
				elsif($ploidy == 4){ # simplified version
					# hyp 0 vs hyp 1 :
					my$lod0vs1 = $ref*log10($dosage[0]/$dosage[1]) + $alt*log10($dosage[4]/$dosage[3]);
					# If lod > 3 --> hyp 0 kept
					if ($lod0vs1 > $LOD){ # hyp 0 > hyp 1
						$dos = 0;
						$dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
					}
					elsif($lod0vs1 < -$LOD){ # hyp 1 > hyp 0
						# hyp 1 vs hyp 2 :
						my$lod1vs2 = $ref*log10($dosage[1]/$dosage[2]) + $alt*log10($dosage[3]/$dosage[2]);
						# si lod > 3 --> hyp 1 kept
						if ($lod1vs2 > $LOD){ # hyp 1 > hyp 2
							$dos = 1;
							$dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
						}
						elsif($lod1vs2 < -$LOD){ # hyp 2 > hyp 1
							# hyp 2 vs hyp 3 :
							my$lod2vs3 = $ref*log10($dosage[2]/$dosage[3]) + $alt*log10($dosage[2]/$dosage[1]);
							if ($lod2vs3 > $LOD){ # hyp 2 > hyp 3
								$dos = 2;
								$dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
							}
							elsif ($lod2vs3 < -$LOD){ # hyp 3 > hyp 2
								# hyp 3 vs hyp 4 :
								my$lod3vs4 = $ref*log10($dosage[3]/$dosage[4]) + $alt*log10($dosage[1]/$dosage[0]);
								if ($lod3vs4 > $LOD){ # hyp 3 > hyp 4
									$dos = 3;
									$dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
								}
								elsif ($lod3vs4 < -$LOD){ # hyp 4 > hyp 3
									$dos = 4;
									$dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
								}
								else{ #indetermination between hyp4 and hyp3
									$dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; # Indetermination 
								}
							}
							else{ #indetermination between hyp3 and hyp2
								$dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; # Indetermination 
							}
						}
						else{ #indetermination between hyp1 and hyp2
							$dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; # Indetermination 
						}
					}
					else{ #indetermination between hyp1 et hyp0
						$dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; # Indetermination 
					}


				}
				elsif($ploidy == 3){ # normal version
				# hyp 0 vs hyp 1 :	
					my$lod0vs1 = $ref*log10($dosage[0]/$dosage[1]) + $alt*log10($dosage[3]/$dosage[2]);
					# If lod > 3 --> hyp 0 kept
					if ($lod0vs1 > $LOD){ # hyp 0 > hyp 1
						$dos = 0;
						$dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
					}
					elsif($lod0vs1 < -$LOD){ # hyp 1 > hyp 0
						# hyp 1 vs hyp 2 :
						my$lod1vs2 = $ref*log10($dosage[1]/$dosage[2]) + $alt*log10($dosage[2]/$dosage[1]);
						# If lod > 3 --> hyp 1 kept
						if ($lod1vs2 > $LOD){ # hyp 1 > hyp 2
							$dos = 1;
							$dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
						}
						# else, hyp 2 vs hyp 3 :
						elsif($lod1vs2 < -$LOD){ # hyp 2 < hyp 1
							my$lod2vs3 = $ref*log10($dosage[2]/$dosage[3]) + $alt*log10($dosage[1]/$dosage[0]);
							# si lod > 3 --> hyp 2 kept
							if ($lod2vs3 > $LOD){ # hyp 2 > hyp 3
								$dos = 2;
								$dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
							}
							# si lod < -3  --> hyp 3 kept
							elsif ($lod2vs3 < -$LOD){ # hyp 3 > hyp 2
								$dos = 3;
								$dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
							}
							else{ 
								$dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
							}
						}
						else{ #indetermination between hyp1 and hyp2
							# hyp1 vs hyp3
							my$lod1vs3 = $ref*log10($dosage[1]/$dosage[3]) + $alt*log10($dosage[2]/$dosage[0]);
							if ($lod1vs3 > $LOD){ # hyp 1 > hyp 3;
								$dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
							}
							else { # hyp 3 > hyp 1 ou indetermination
								# hyp 2 vs hyp 3
								my$lod2vs3 = $ref*log10($dosage[2]/$dosage[3]) + $alt*log10($dosage[1]/$dosage[0]);
								# si lod > 3 --> hyp2 kept
								if ($lod2vs3 > $LOD){ # hyp 2 > hyp 3
									$dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
								}
								# si lod < -3  --> hyp3 kept
								elsif ($lod2vs3 < -$LOD){ # hyp 3 > hyp 2
									$dos = 3;
									$dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
								}
								else{ 
									$dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
								}
							}
						}
					}
					else {  #indetermination between hyp0 and hyp1
						my$lod0vs2 = $ref*log10($dosage[0]/$dosage[2]) + $alt*log10($dosage[3]/$dosage[1]);
					# hyp 0 vs hyp 2 :
						if ($lod0vs2 > $LOD){ # hyp 0 > hyp 2
							$dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
						}
						# lod < 3 --> hyp 2 kept but we have to test others hypothesis.
						elsif($lod0vs2 < -$LOD){ # hyp 2 > hyp 0
							# we test hyp 1 vs hyp2 :
							my$lod1vs2 = $ref*log10($dosage[1]/$dosage[2]) + $alt*log10($dosage[2]/$dosage[1]);
							# if lod > 3 --> hyp 1 kept
							if ($lod1vs2 > $LOD){ # hyp 1 > hyp 2
								$dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
							}
							# else, hyp 2 vs hyp 3 :
							elsif($lod1vs2 < -$LOD){ # hyp 2 < hyp 1
								my$lod2vs3 = $ref*log10($dosage[2]/$dosage[3]) + $alt*log10($dosage[1]/$dosage[0]);
							
								# if lod > 3 --> hyp 2 kept
								if ($lod2vs3 > $LOD){ # hyp 2 > hyp 3
									$dos = 2;
									$dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
								}
								# si lod < -3  --> hyp 3 kept
								elsif ($lod2vs3 < -$LOD){ # hyp 3 > hyp 2
									$dos = 3;
									$dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
								}
								else{ 
									$dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
								}
							}
							else{ #indetermination between hyp1 and hyp2
								# then hyp1 vs hyp3
								my$lod1vs3 = $ref*log10($dosage[1]/$dosage[3]) + $alt*log10($dosage[2]/$dosage[0]);
								if ($lod1vs3 > $LOD){ # hyp 1 > hyp 3
									$dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
								}
								elsif ($lod1vs3 < -$LOD){
									# hyp 2 and hyp 3
									my$lod2vs3 = $ref*log10($dosage[2]/$dosage[3]) + $alt*log10($dosage[1]/$dosage[0]);
									# if lod > 3 --> hyp 2 kept
									if ($lod2vs3 > $LOD){ # hyp 2 > hyp 3
										$dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
									}
									# if lod < -3  --> hyp 3 kept
									elsif ($lod2vs3 < -$LOD){ # hyp 3 > hyp 2
										$dos = 3;
										$dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
									}
									else{ 
										$dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
									}
								}
								else { # indetermination 1v3
									$dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA";	
								}
							}
						}
						else{ # indetermination between 0 and 2
							# hyp 0 vs hyp 3
							my$lod0vs3 = $ref*log10($dosage[0]/$dosage[3]) + $alt*log10($dosage[3]/$dosage[0]);
							if ($lod0vs3 > 3){ # hyp 0 > hyp 3
								$dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
							}
							elsif ($lod0vs3 < -$LOD){ # hyp 3 > hyp 0 --> we have to test other hypothesis
								# hyp 1 vs hyp 3 :
								my$lod1vs3 = $ref*log10($dosage[1]/$dosage[3]) + $alt*log10($dosage[2]/$dosage[0]);
								# if lod > 3 --> ID
								if ($lod1vs3 > $LOD){ # hyp 1 > hyp 3
									$dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
								}
								# else, hyp 2 vs hyp 3 :
								elsif($lod1vs3 < -$LOD){ # hyp 2 < hyp 1
									my$lod2vs3 = $ref*log10($dosage[2]/$dosage[3]) + $alt*log10($dosage[1]/$dosage[0]);
									# if lod > 3 --> hyp 2 kept
									if ($lod2vs3 > $LOD){ # hyp 2 > hyp 3
										$dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
									}
									# if lod < -3  --> hyp 3 kept
									elsif ($lod2vs3 < -$LOD){ # hyp 3 > hyp 2
										$dos = 3;
										$dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
									}
									else{ 
										$dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
									}
								}
								else{ #indetermination between hyp1 and hyp3
									$dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
								}
							}
							else{ # indétermintation between 0 and 3
								$dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
							}
						}
					}
				}
			}
		}
	}
    }
}

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Lenght of chromosome determination.
my%theEndChr; 
my%theNumChr; 
if ($booEndChro == 0){
	my$n = 0;
	for my$x(sort keys %endChro){
		if (grep { $_ eq $x } @chromosomList){
			$n++;
			$theEndChr{$x} = $endChro{$x};
			$theNumChr{$x} = $n;
		}
	}
	%endChro = %theEndChr;
}
elsif ($booEndChro == 1){
	my$n = 0;
	for my$x(sort keys %realEndChro){
		if (grep { $_ eq $x } @chromosomList){
			$n++;
			$theEndChr{$x} = $realEndChro{$x};
			$theNumChr{$x} = $n;
		}
	}
	%endChro = %theEndChr;
}

# step of keys inversion (parent/chrom/num) of %blocs to (chrom/parent/num) %blocsInv 
my%blocsInv;
for my$parent(sort keys %blocspos){
	for my $CHROMOSOMES(sort keys %{$blocspos{$parent}}){
		for my $num (sort {$a<=>$b} keys %{$blocspos{$parent}{$CHROMOSOMES}}){
			$blocsInv{$CHROMOSOMES}{$parent}{$num} = $blocspos{$parent}{$CHROMOSOMES}{$num};
		}
	}
}

# fill little windows with allelic dosages. The sum of every ancestors for one little window should be equal to the ploidy. Else it's an indetermination -> NA 

my$cut=$window_cut*1000;
my$startslice=1;
my$endslice=$cut;
my%dosbywindows;
my$sum;
for my$ind (sort keys %dosageAllelique){
	for my$CHROMOSOMES(sort keys %blocsInv){
		$startslice=1;
		$endslice=$cut;
		while ($startslice < $endChro{$CHROMOSOMES}){
			for my $parent(sort keys %{$blocsInv{$CHROMOSOMES}}){
				$sum = 0;
				for my$num (sort {$a<=>$b} keys %{$blocsInv{$CHROMOSOMES}{$parent}}){
					my@markers = split("-",$blocsInv{$CHROMOSOMES}{$parent}{$num});
					my$start = $markers[0];
					my$end = $markers[1];
					if($startslice >= $start && $startslice < $end){
						$dosbywindows{$ind}{$CHROMOSOMES}{$startslice}{$parent} = $dosageAllelique{$ind}{$parent}{$CHROMOSOMES}{$num};
					}
				}
			}
			$startslice = $endslice;
			$endslice = $startslice + $cut;
		}
	}

}

for my$ind (sort keys %dosbywindows){
	for my$CHROMOSOMES(sort keys %{$dosbywindows{$ind}}){
		for my$startslice(sort {$a<=>$b} keys %{$dosbywindows{$ind}{$CHROMOSOMES}}){
			foreach my$parent(@ref){
				if (!exists($dosbywindows{$ind}{$CHROMOSOMES}{$startslice}{$parent})){
					$dosbywindows{$ind}{$CHROMOSOMES}{$startslice}{$parent} = "NE"; 
				} 
			}
		}
	}
}

my%color; # 1 color by ancestor
open(Fcolor, "$color") or die ("Erreur d'ouverture de Fcolor\n");
my$count_col = 0; 
while (my $li = <Fcolor>){
	chomp($li);
	if ($li =~ m/\S+\s+\S+/){
		my@paint = split(/\s+/, $li);
		$color{$paint[0]}=$paint[1];
		$count_col++;
	}
}
#$color{"NA"} = "#B9B9B9"; # grey for NA
$color{"separator"} = "#000000"; # black for zones inter-chromosomes
if (scalar(@ref)+1 != $count_col){
	print "the NA color will be grey\n";
	if(scalar(@ref) == $count_col){
		$color{"NA"} = "#B9B9B9"; # grey for NA --> if the user have not chosen the NA color in the color file
	}
	else{
		print "there is not enougth colors indicated for the number of ancestors\n";
	}
}
close Fcolor;

my%hasher;
if($hybrid){ ## If focus on hybrid
	my$ind = $hybrid;
	for my $CHROMOSOMES (sort keys %{$dosbywindows{$ind}}){
		for my $slice (sort {$a<=>$b} keys %{$dosbywindows{$ind}{$CHROMOSOMES}}){
			my$booNa = 0; # boolean if NA for at least one ancestor in the slice
			my%dosparent; # allelic dosage by parent for this slice
			my$dosBySlice = 0; # sum of allelic dosage for a slice
			for my $parent (sort keys %{$dosbywindows{$ind}{$CHROMOSOMES}{$slice}}){
				if ($dosbywindows{$ind}{$CHROMOSOMES}{$slice}{$parent} eq "NE"){ # NE = blocs with no corresponding markers in vcf => Dosage = 0
					$dosbywindows{$ind}{$CHROMOSOMES}{$slice}{$parent} = 0;

					if ($booNa == 0){
						$dosparent{$parent} = $dosbywindows{$ind}{$CHROMOSOMES}{$slice}{$parent};
						$dosBySlice = $dosBySlice + $dosparent{$parent};
					}
					elsif($booNa == 1){
						$dosBySlice = "NA";
					}			
				}
				elsif($dosbywindows{$ind}{$CHROMOSOMES}{$slice}{$parent} eq "NA"){ # NA = bloc with LOD indetermintation -> NA
					$dosparent{$parent} = "NA";
					$booNa = 1;
					$dosBySlice = "NA";
				}
				else{ # real dosage
					if ($booNa == 0){
						$dosparent{$parent} = $dosbywindows{$ind}{$CHROMOSOMES}{$slice}{$parent};
						$dosBySlice = $dosBySlice + $dosparent{$parent};
					}
					elsif($booNa == 1){
						$dosBySlice = "NA";
					}
				}
			}
			# fill $hasher with color for each slice
			# If NA of != ploidy --> NA. Can't conclude anything.
			if ($dosBySlice eq "NA"){
				for (my$p = 0; $p < $ploidy; $p++){
					my$slice2 = $slice + $cut;
					$hasher{$CHROMOSOMES}{$p}{$slice} = $color{"NA"};
				}
			}
			elsif ($dosBySlice != $ploidy){
				for (my$p=0; $p < $ploidy; $p++){
					my$slice2 = $slice + $cut;
					$hasher{$CHROMOSOMES}{$p}{$slice} = $color{"NA"};
				}
			}
			# Else, if = ploidy. We can paint.
			elsif ($dosBySlice == $ploidy){
				my$slice2 = $slice + $cut;
				my$p = 0; # incrementation of p until the poidy each time we paint an haplotype
				for my$parent (sort keys %dosparent){
					my$x = $dosparent{$parent}; # We paint as much haplotype than allelic dosage for this parent
					if ($x != 0){
						for(my$i = 1; $i <= $x; $i++){
							if ($p == 0){
								$hasher{$CHROMOSOMES}{$p}{$slice} = $color{$parent};
							}
							elsif ($p == 1){
								$hasher{$CHROMOSOMES}{$p}{$slice} = $color{$parent};
							}
							elsif ($p == 2){
								$hasher{$CHROMOSOMES}{$p}{$slice} = $color{$parent};
							}
							elsif ($p == 3){
								$hasher{$CHROMOSOMES}{$p}{$slice} = $color{$parent};
							}
							$p++;
						}
					}
				}
				
			}
		}
	}
	open(F11, ">ideogram_$hybrid") or die ("Error: can't open output for ideogram\n"); # for a particular hybrid. Output = ideogram format.
	for my$CHROMOSOMES (sort keys %hasher){	
		my$endC = $endChro{$CHROMOSOMES};
		my$numchr = $theNumChr{$CHROMOSOMES};
		for my$p (sort {$a<=>$b} keys %{$hasher{$CHROMOSOMES}}){
			my$color = "";
			my$begin = 1;
			my$endS = 0;
			my$sli;
			for my$slice (sort {$a<=>$b} keys %{$hasher{$CHROMOSOMES}{$p}}){
				if ($color eq ""){
					$color = $hasher{$CHROMOSOMES}{$p}{$slice};
				}
				# we join slice of the same color together before we print
				if ($color ne $hasher{$CHROMOSOMES}{$p}{$slice}){
					$endS = $slice;
					if ($endS <= $endC){
						print F11 "$numchr $p $begin $endS $color\n";
					}
					else{
						print F11 "$numchr $p $begin $endC $color\n";
					}
					$begin = $slice + 1;
				}
				$color = $hasher{$CHROMOSOMES}{$p}{$slice};
				$sli = $slice;
			}
			if ($endS < $endC){
				if ($sli + $cut > $endS){
					if ($sli + $cut >= $endC){ # If the chromosome longer than the end of the last slice -> grey 'til the end 
						$endS =	$endC;
						print F11 "$numchr $p $begin $endS $color\n"; #color of the last slices
					}	
					else{
						if($color ne $color{"NA"}){
							$endS =	$sli + $cut;
							print F11 "$numchr $p $begin $endS $color\n"; #color of the last slices
							$begin = $endS + 1;
							$color = $color{"NA"};
							print F11 "$numchr $p $begin $endC $color\n"; #end of chro grey
						}
						else{ #if color = color of NA
							print F11 "$numchr $p $begin $endC $color\n";
						}
					}
				}
				else {
					$begin = $endS + 1;
					$color = $color{"NA"};
					print F11 "$numchr $p $begin $endC $color\n"; #end of chro grey
				}
			}
			
		}
	}
	close F11;
	# Print the length file.
	# chr length_chromosome haplotypes_names
	open(FIL, ">len_ideogram_$hybrid") or die ("Error : can't open length output for ideogram\n");
	for my$ch(sort keys %endChro){
		if ($ploidy == 2){
			print FIL "$ch $endChro{$ch} 01\n";
		}
		if ($ploidy == 3){
			print FIL "$ch $endChro{$ch} 012\n";
		}
		if ($ploidy == 4){
			print FIL "$ch $endChro{$ch} 0123\n";
		}
	}
}
else{ # If there is no focus on an hybrid: same things but on every individuals
	for my$ind (sort keys %dosbywindows){
		for my $CHROMOSOMES (sort keys %{$dosbywindows{$ind}}){
			for my $slice (sort {$a<=>$b} keys %{$dosbywindows{$ind}{$CHROMOSOMES}}){
				my$booNa=0;
				my%dosparent;
				my$dosBySlice = 0;
				for my $parent (sort keys %{$dosbywindows{$ind}{$CHROMOSOMES}{$slice}}){
					if ($dosbywindows{$ind}{$CHROMOSOMES}{$slice}{$parent} eq "NE"){ #NE = allelic dosage -> 0
						$dosbywindows{$ind}{$CHROMOSOMES}{$slice}{$parent} = 0;
						if ($booNa == 0){
							$dosBySlice = $dosBySlice + $dosbywindows{$ind}{$CHROMOSOMES}{$slice}{$parent};
							$dosparent{$parent} = $dosbywindows{$ind}{$CHROMOSOMES}{$slice}{$parent};
						}
						elsif($booNa == 1){
							$dosBySlice = "NA";
						}			
					}
					elsif($dosbywindows{$ind}{$CHROMOSOMES}{$slice}{$parent} eq "NA"){ #NA = indetermination of LOD -> NA
						$dosparent{$parent} = "NA";
						$booNa = 1;
						$dosBySlice = "NA";
					}
					else{ # real dosage
						if ($booNa == 0){
							$dosBySlice = $dosBySlice + $dosbywindows{$ind}{$CHROMOSOMES}{$slice}{$parent};
							$dosparent{$parent} = $dosbywindows{$ind}{$CHROMOSOMES}{$slice}{$parent};
						}
						elsif($booNa == 1){
							$dosBySlice = "NA";
						}
					}
				}
				# If NA of != ploidy --> NA. Can't conclude anything.
				if ($dosBySlice eq "NA"){
					for (my$p=0; $p < $ploidy; $p++){
						$hasher{$ind}{$CHROMOSOMES}{$p}{$slice} = $color{"NA"};
					}
				}
				elsif ($dosBySlice != $ploidy){
					for (my$p=0; $p < $ploidy; $p++){
						$hasher{$ind}{$CHROMOSOMES}{$p}{$slice} = $color{"NA"};
					}
				}
				# else -> painting
				elsif ($dosBySlice == $ploidy){
					my$p = 0;
					for my$parent (sort keys %dosparent){
						my$x = $dosparent{$parent};
						if ($x != 0){
							for(my$i = 1; $i <= $x; $i++){
								if ($p == 0){
									$hasher{$ind}{$CHROMOSOMES}{$p}{$slice} = $color{$parent};
								}
								elsif ($p == 1){
									$hasher{$ind}{$CHROMOSOMES}{$p}{$slice} = $color{$parent};
								}
								elsif ($p == 2){
									$hasher{$ind}{$CHROMOSOMES}{$p}{$slice} = $color{$parent};
								}
								elsif ($p == 3){
									$hasher{$ind}{$CHROMOSOMES}{$p}{$slice} = $color{$parent};
								}
								$p++;
							}
						}
					}
				}
			}
		}
	}
	# Write file compatible with circos.js
	open(F12, ">circos") or die ("Error : can't write on circos output file\n");
	if ($indivs){ #If there is a focus file, the output will be on the individuals focused
		open(Findivs, "$indivs") or die ("can't open hybrids file");
		my%indivs; #hash of focus file : key = individu name , value = number
		while (my $li = <Findivs>){
			chomp($li);
			my@i=split(/\s+/,$li);
			$indivs{$i[1]} = $i[0];
		}
		close Findivs;
		for my$ind (sort keys %indivs){
			for my$CHROMOSOMES (sort keys %{$hasher{$ind}}){
				my$endC = $endChro{$CHROMOSOMES};
				for my$p (sort {$a<=>$b} keys %{$hasher{$ind}{$CHROMOSOMES}}){
					my$sli;
					my$endS = 0;
					my$color = "";
					my$begin = 1;
					for my$slice (sort {$a<=>$b} keys %{$hasher{$ind}{$CHROMOSOMES}{$p}}){
						if ($color eq ""){
							$color = $hasher{$ind}{$CHROMOSOMES}{$p}{$slice};
						}
						# Join slices of the same color before printing
						if ($color ne $hasher{$ind}{$CHROMOSOMES}{$p}{$slice}){
							$endS = $slice;
							if ($endS <= $endC){
								print F12 "$CHROMOSOMES $begin $endS $color $ind\n";
							}
							else{
								print F12 "$CHROMOSOMES $begin $endC $color $ind\n";
							}
							$begin = $slice + 1;
						}
						$color = $hasher{$ind}{$CHROMOSOMES}{$p}{$slice};
						$sli = $slice;
					}
					if ($endS < $endC){ # end of chromosome = grey
						if ($sli + $cut > $endS){
							if ($sli + $cut >= $endC){ # If the chromosome longer than the end of the last slice -> grey 'til the end 
								$endS =	$endC;
								print F12 "$CHROMOSOMES $begin $endS $color $ind\n"; #color of the last slices
							}	
							else{
								if($color ne $color{"NA"}){
									$endS =	$sli + $cut;
									print F12 "$CHROMOSOMES $begin $endS $color $ind\n"; #color of the last slices
									$begin = $endS + 1;
									$color = $color{"NA"};
									print F12 "$CHROMOSOMES $begin $endC $color $ind\n"; #end of chro grey
								}
								else{
									print F12 "$CHROMOSOMES $begin $endC $color $ind\n"; #end of chro grey
								}
							}
						}
						else {
							$begin = $endS + 1;
							$color = $color{"NA"};
							print F12 "$CHROMOSOMES $begin $endC $color $ind\n"; #end of chro grey
						}
					}
				}
			}
		}
	}
	else{ # If there is no focus file, the output will be on all the individuals : same thing.
		for my$ind (sort keys %hasher){
			for my$CHROMOSOMES (sort keys %{$hasher{$ind}}){
				my$endC = $endChro{$CHROMOSOMES};
				for my$p (sort {$a<=>$b} keys %{$hasher{$ind}{$CHROMOSOMES}}){
					my$sli;
					my$endS = 0;
					my$color = "";
					my$begin = 1;
					for my$slice (sort {$a<=>$b} keys %{$hasher{$ind}{$CHROMOSOMES}{$p}}){
						if ($color eq ""){
							$color = $hasher{$ind}{$CHROMOSOMES}{$p}{$slice};
						}
						if ($color ne $hasher{$ind}{$CHROMOSOMES}{$p}{$slice}){
							$endS = $slice ;
							if ($endS <= $endC){
								print F12 "$CHROMOSOMES $begin $endS $color $ind\n";
							}
							else{
								print F12 "$CHROMOSOMES $begin $endC $color $ind\n";
							}
							$begin = $slice + 1;
						}
						$color = $hasher{$ind}{$CHROMOSOMES}{$p}{$slice};
						$sli = $slice;
					}
					if ($endS < $endC){ # end of chromosome = grey
						if ($sli + $cut > $endS){
							if ($sli + $cut >= $endC){ # If the chromosome longer than the end of the last slice -> grey 'til the end 
								$endS =	$endC;
								print F12 "$CHROMOSOMES $begin $endS $color $ind\n"; #color of the last slices
							}	
							else{
								if($color ne $color{"NA"}){
									$endS =	$sli + $cut;
									print F12 "$CHROMOSOMES $begin $endS $color $ind\n"; #color of the last slices
									$begin = $endS + 1;
									$color = $color{"NA"};
									print F12 "$CHROMOSOMES $begin $endC $color $ind\n"; #end of chro grey
								}
								else{
									print F12 "$CHROMOSOMES $begin $endC $color $ind\n"; #end of chro grey
								}
							}
						}
						else {
							$begin = $endS + 1;
							$color = $color{"NA"};
							print F12 "$CHROMOSOMES $begin $endC $color $ind\n"; #end of chro grey
						}
					}
				}
			}
		}
	}
	# Output of ciros file lenght
	open(FLCirc, ">len_allInd_Circos") or die ("Error : can't write circos file length\n");
	for my$ch(sort keys %endChro){
		print FLCirc "$ch $endChro{$ch}\n";
	}
	close FLCirc;
	close F12;
	# Output = ideogram file but for several hybrid. All chromosomes are put side by side separated by a black zone.
	open(F13, ">ideogram_allInd") or die ("Error : can't write ideogram output file\n");
	open(F14, ">len_allInd_Ideogram") or die ("Error : can't write ideogram output length file\n");
	# If there is a focus file, the output will be on the individuals focused
	if ($indivs){
		open(Findivs, "$indivs") or die ("can't open hybrids file");
		my%indivs; #hash of focus file : key = individu name , value = number
		while (my $li = <Findivs>){
			chomp($li);
			my@i=split(/\s+/,$li);
			$indivs{$i[1]} = $i[0];
		}
		close Findivs;
		for my$ind (sort keys %indivs){
			my$indNum = $indivs{$ind};
			my$endC;
			my$endS;	
			for (my$p = 0; $p < $ploidy; $p++){
				my$sli;
				my$begin = 1;
				$endS = 0;
				$endC = 0;
				for my$CHROMOSOMES (sort keys %{$hasher{$ind}}){
					my$color = "";
					for my$slice (sort {$a<=>$b} keys %{$hasher{$ind}{$CHROMOSOMES}{$p}}){
						if ($color eq ""){
							$color = $hasher{$ind}{$CHROMOSOMES}{$p}{$slice};
						}
						# join slice of the same color before printing
						if ($color ne $hasher{$ind}{$CHROMOSOMES}{$p}{$slice}){
							$endS = $endC + $slice;
							if ($endS <= $endC + $endChro{$CHROMOSOMES}){	
								print F13 "$indNum $p $begin $endS $color\n";
							}
							else{
								print F13 "$indNum $p $begin $endC $color\n";
							}
							$begin = $endC + $slice + 1;
						}
						$color = $hasher{$ind}{$CHROMOSOMES}{$p}{$slice};
						$sli = $endC + $slice;
					}
					$endC = $endC + $endChro{$CHROMOSOMES};
					if ($endS <= $endC){ # NA until the end of chromosome (grey) if no information
						if ($sli + $cut > $endS){
							if ($sli + $cut >= $endC){ 
								$endS = $endC;
								print F13 "$indNum $p $begin $endS $color\n";
								$begin = $endC + 1;
								$endS = $begin + $cut * 20;
								$endC = $endC + $cut * 20; # The total size of the big "chromosome" = all chromosomes length + black zones
								$color = $color{"separator"}; # little black zone to separate different chromosomes
								print F13 "$indNum $p $begin $endS $color\n"; 
								$begin = $endS + 1;
							}
							else{ # If the chromosome longer than the end of the last slice -> grey 'til the end 
								$endS = $sli + $cut;
								print F13 "$indNum $p $begin $endS $color\n";
								$begin = $endS + 1;
								$color = $color{"NA"};
								print F13 "$indNum $p $begin $endC $color\n";
								$begin = $endC + 1;
								$endS = $begin + $cut * 20;
								$endC = $endC + $cut * 20; # The total size of the big "chromosome" = all chromosomes length + black zones
								$color = $color{"separator"}; # little black zone to separate different chromosomes
								print F13 "$indNum $p $begin $endS $color\n"; 
								$begin = $endS + 1;
							}
						}
						else{
							$begin = $endS + 1;
							$color = $color{"NA"};
							print F13 "$indNum $p $begin $endC $color\n";
							$begin = $endC + 1;
							$endS = $begin + $cut * 20;
							$endC = $endC + $cut * 20; # The total size of the big "chromosome" = all chromosomes length + black zones
							$color = $color{"separator"}; # little black zone to separate different chromosomes
							print F13 "$indNum $p $begin $endS $color\n"; 
							$begin = $endS + 1;
						}
					}
				}
			}
			if ($ploidy == 2){
				print F14 "$indNum $endS 01\n";
			}
			if ($ploidy == 3){
				print F14 "$indNum $endS 012\n";
			}
			if ($ploidy == 4){
				print F14 "$indNum $endS 0123\n";
			}
		}
	}
	#If there is no focus file, the output will be on all the individuals
	else{
		my$blop = 0; #indice of each hybrid
		for my$ind (sort keys %hasher){
			$blop++;
			my$endC;
			my$endS;	
			for (my$p = 0; $p < $ploidy; $p++){
				my$sli;
				my$begin = 1;
				$endS = 0;
				$endC = 0;
				for my$CHROMOSOMES (sort keys %{$hasher{$ind}}){
					my$color = "";
					for my$slice (sort {$a<=>$b} keys %{$hasher{$ind}{$CHROMOSOMES}{$p}}){
						if ($color eq ""){
							$color = $hasher{$ind}{$CHROMOSOMES}{$p}{$slice};
						}
						# join slice of the same color before printing
						if ($color ne $hasher{$ind}{$CHROMOSOMES}{$p}{$slice}){
							$endS = $endC + $slice;
							if ($endS <= $endC + $endChro{$CHROMOSOMES}){
								print F13 "$blop $p $begin $endS $color\n";
							}
							else{
								print F13 "$blop $p $begin $endC $color\n";
							}
							$begin = $endC + $slice + 1;
						}
						$color = $hasher{$ind}{$CHROMOSOMES}{$p}{$slice};
						$sli = $endC + $slice;
					}
					$endC = $endC + $endChro{$CHROMOSOMES};
					if ($endS <= $endC){ # NA until the end of chromosome (grey) if no information
						if ($sli + $cut > $endS){
							if ($sli + $cut >= $endC){ 
								$endS = $endC;
								print F13 "$blop $p $begin $endS $color\n";
								$begin = $endC + 1;
								$endS = $begin + $cut * 20;
								$endC = $endC + $cut * 20; # The total size of the big "chromosome" = all chromosomes length + black zones
								$color = $color{"separator"}; # little black zone to separate different chromosomes
								print F13 "$blop $p $begin $endS $color\n"; 
								$begin = $endS + 1;
							}
							else{ # If the chromosome longer than the end of the last slice -> grey 'til the end 
								$endS = $sli + $cut;
								print F13 "$blop $p $begin $endS $color\n";
								$begin = $endS + 1;
								$color = $color{"NA"};
								print F13 "$blop $p $begin $endC $color\n";
								$begin = $endC + 1;
								$endS = $begin + $cut * 20;
								$endC = $endC + $cut * 20; # The total size of the big "chromosome" = all chromosomes length + black zones
								$color = $color{"separator"}; # little black zone to separate different chromosomes
								print F13 "$blop $p $begin $endS $color\n"; 
								$begin = $endS + 1;
							}
						}
						else{
							$begin = $endS + 1;
							$color = $color{"NA"};
							print F13 "$blop $p $begin $endC $color\n";
							$begin = $endC + 1;
							$endS = $begin + $cut * 20;
							$endC = $endC + $cut * 20; # The total size of the big "chromosome" = all chromosomes length + black zones
							$color = $color{"separator"}; # little black zone to separate different chromosomes
							print F13 "$blop $p $begin $endS $color\n"; 
							$begin = $endS + 1;
						}
					}
				}
			}
			if ($ploidy == 2){
				print F14 "$blop $endS 01\n";
			}
			if ($ploidy == 3){
				print F14 "$blop $endS 012\n";
			}
			if ($ploidy == 4){
				print F14 "$blop $endS 0123\n";
			}
		}
	}
	close F13;
	close F14;
}
print("finish\n");
