#!/usr/bin/perl 

=for comment

	FILE		mlst.pl - for create a tree of distances between types
	AUTHOR 		Verdier Axel
	CREATE_DATE	04-22-2014
	LAST_DATE	05-27-2014
=cut

use strict;
use File::Basename;
use FindBin;              # locate this script
use lib "$FindBin::Bin";  # use the current directory
my $abs_path = $FindBin::Bin;
use moduls::FIT;#1st step
use moduls::DMC;#2nd step
use moduls::CTO;#3st step
use moduls::readOptions;

use Getopt::Long qw(:config bundling);
use Pod::Usage;

print "\t[== PRE-STEP ==]\n";

#Read Options
print "Reading Options";
$|=1;
my $DATAS_MCCP;
my $DATAS_CBPP;
my $NBLOCUS;
my $MINGAPSIZE;
my %DISTANCE_VALUE;
readOptions::readOptions($abs_path, $DATAS_MCCP, $DATAS_CBPP, $NBLOCUS, $MINGAPSIZE, \%DISTANCE_VALUE);
print " - DONE\n";

#Parameters
print "Reading Parameters";
$|=1;

my $man = 0;
my $help = 0;
my $inputPath='';
my $outputMatrix='';
my $outputTree='';
my @inSequences=(); 
my $name='';
my $disease=0;#Default 0: mccp
GetOptions(	'help|h|?' => \$help,
		man => \$man,
		'file|f=s' => \$inputPath,
		'sequences|s=s' => \@inSequences,
		'name|n=s' => \$name,
		'disease|d=s' => \$disease,
		'oMatrix|matrix|m=s'=> \$outputMatrix,
		'oTree|tree|t=s'=>\$outputTree)
	or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitval => 0, -verbose => 2) if $man;
	#Inputs sequence or file
if($inputPath && @inSequences or (!$inputPath && !@inSequences)){die("\nERROR : Must have an input file OR inputs sequences\n");}
my $inputName="";
if (@inSequences != 0){
	@inSequences = split(/,/,join(',',@inSequences));
	if(@inSequences>$NBLOCUS){die("\nERROR : need $NBLOCUS sequences (1 per locus), there are ".@inSequences."\n");}
	elsif(@inSequences<$NBLOCUS){for (my $i=@inSequences; $i<$NBLOCUS; $i++){$inSequences[$i]='';}}#if there is not enough seq in the inputs (no amplification ?), we add empty sequences
	for (my $i=0; $i<@inSequences; $i++){
		$inSequences[$i]=~s/\s//g; #To delete all space, newlines and tabulations
	}
	if(!$name){$name="yours";}
}
else{
	$inputName=basename $inputPath;
#	if ($inputName!~m/^.*\.fa(s)?$/ and $inputName!~m/^.*\.fasta$/){die("\nfichier $inputName invalide ! Il doit être de format fasta (.fa, .fas ou .fasta)\n");}
	if(!$name){$name=$inputName;}
}
	#Disease
if($disease eq "1" or $disease eq "cbpp" or $disease eq "CBPP3"){$disease=1;}
else{$disease=0;}

print " - DONE\n";

print "\n\t[== STEP 1 ==]\n";
#Read file
print "Reading File !";
$|=1;
my @inType=();
#if it's the file
if($inputPath){@inType=FIT::readTypeFile($inputPath, $NBLOCUS);}
#else it's sequences
else{
	for(my $i=0; $i<@inSequences; $i++){
		$inType[$i]={"name" => "loc_$i", "sequence" => $inSequences[$i]};
	}
}
print " - DONE\n";

#Create the matrix of alleles per locus
print "Alleles matrix creation !";
$|=1;
my $path;
if($disease==0){$path=$DATAS_MCCP;}
elsif($disease==1){$path=$DATAS_CBPP;}
else{
	print "INTERN_ERROR : assign \$path, please report it\n";
	exit;
}
my @allelesLocus=FIT::createAllelesLocus($path, $NBLOCUS);
print " - DONE\n";

#Create the matrix of alleles distances per locus
print "Alleles distances matrix creation !";
$|=1;
my @allelesDist=FIT::createAllelesDist($path, $NBLOCUS, $MINGAPSIZE, \%DISTANCE_VALUE);
print " - DONE\n";

#Search the inType Locus allele id
print "Search alleles of file entry correspondences !";
$|=1;
FIT::identifyAlleles(\@inType,\@allelesLocus, \@allelesDist, \%DISTANCE_VALUE, $NBLOCUS, $MINGAPSIZE);
print " - DONE\n";

print "\n\t[== STEP 2 ==]\n";
#Create the types alleles matrix
print "Create type alleles matrix !";
$|=1;
$path=~s/\..*$//;
$path.="-matrix.csv";
my @nameTypes;
my @typesMatrix = DMC::createTypesMatrix($path, $NBLOCUS, \@nameTypes);
my $sizeMatrix=@typesMatrix;
my $newType=1;
my $inTypeId=0;
my $n=0;
	#Search the id of the in-type and add it in the types matrix if necessary
while ($newType==1 && $n<$sizeMatrix){
	for (my $i=0; $i<$NBLOCUS; $i++){
		if ($typesMatrix[$n][$i] != $inType[$i]{id}){$i=$NBLOCUS;}
		else{
			if ($i == $NBLOCUS-1){
				$newType=0;
				$inTypeId=$n;
				$nameTypes[$n]=$nameTypes[$n]."_".$name;
			}  
		}
	}
	$n++;
}
if ($newType == 1){
		for (my $i=0; $i<$NBLOCUS; $i++){
			$typesMatrix[$sizeMatrix][$i]=$inType[$i]{id};
			
		}
		$nameTypes[$sizeMatrix]=$name;
		$inTypeId=$sizeMatrix;
}
print " - DONE\n";
if ($newType){print "It's an unknow type.\tPlease contact\n";}
else{print "This type is similar to an other type : id=$inTypeId\n";}
#write the matrix if necessary
if($outputMatrix){
	print "Write alleles matrix !";
	$|=1;
	my $fileContent="";
	for(my $loc=0; $loc< $NBLOCUS; $loc++){
		my ($nameLoc)=$allelesLocus[$loc][0]{name}=~m/^(.*)[\-_][Aa]l.*$/;
		$fileContent.="\t$nameLoc";
	}
	for (my $i=0; $i<@typesMatrix; $i++){
		$fileContent.="\n".$nameTypes[$i];
		for(my $loc=0; $loc< $NBLOCUS; $loc++){
			$fileContent.="\t".($typesMatrix[$i][$loc]+1);
		}
	}
	open(MATRIXFILE, '>', $outputMatrix);
	print MATRIXFILE $fileContent;
	close(MATRIXFILE);
	print " - DONE\n";
}

#Create the distances matrix of types
print "Create types distance matrix !";
$|=1;
my @matrixDist=DMC::createDistMatrix(\@typesMatrix, \@allelesDist, $NBLOCUS);
print " - DONE\n";

print "\n\t[== STEP 3 ==]\n";
#Create the tree
print "Calculate the tree !";
$|=1;
my $tree=CTO::calcTree(\@matrixDist, \@nameTypes, $inTypeId);
print " - DONE\n";
$|=1;
#write the tree
if($outputTree){
	print "Create newick file !";
	open(FILE, '>', $outputTree);
	print FILE $tree->as_text('newick');
	close(FILE);
	print " - DONE\n";
}

print "\n\t[== END ==]\n";

__END__

=head1 NAME

mlst.pl - Compare an in-type with the know type of the diseases mccp or cbpp. Return a matrix and a newick tree

=head1 SYNOPSIS

mlst.pl [options]

 Options:
   --help	brief help message
   --man	full documentation
   --file	the input path of the multi-fasta file
   --sequences	the alleles sequences of the type to test
   --oMatrix	the output path for the types alleles id matrix (.csv)
   --oTree	the output path for the newick tree (.nwk)

=head1 OPTIONS

=over 8

=item B<--help | -h | -?>

Print a brief help message and exits.

=item B<--man>

Prints the manual page and exits.

=item B<--file | -f>

Path of the input multi-fasta file

=item B<--sequences | -s>

The sequences separate with ','. If there is an empty sequence ('no amplification'), you should have ',,'.

=item B<--oMatrix | --matrix | -m>

Path for the ouput matrix

=item B<--oTree | --tree | -t>

Path for the output tree

=back

=head1 DESCRIPTION

This program will read the given input multi-fasta file of the sequences B<or> the sequences of each locus (8) of the type and compare the type to all known. You can see a distances tree of this results.

=head1 EXAMPLES

=over 8

=item with a Multi-fasta file

C<./mlst --file path/to/the/multifasta_file.fa(s) --disease cbpp --matrix path/to/where/you/want/the/matrix --tree path/to/where/you/want/the/tree>

=item directly the sequences with no amplification for the 7th locus

C<./mlst --sequences ATGCGTA[...]GTCATGA,GGATA[...],[...],[...],[...],[...],,[...] -d mccp -m path/for/the/matrix -t path/for/the/tree>

=back

=cut



