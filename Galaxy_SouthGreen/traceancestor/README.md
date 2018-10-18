
TraceAncestor.pl is a script that permits to estimate the allelic dosage of ancestral alleles in hybrid individuals and then to perform chromosom painting.


# TUTORIEL

## TAprefilter.pl

### Usage

This script is used to define a matrix of ancestry informative markers from an other matrix containing GST (inter-population differentiation parameter) information.

	perl TAprefilter.pl [-t matrix file] [-g threshold for gst] [-m threshold for missing data]

### Inputs

	-t | --input : reference matrice with GST and alternative allele frequence (F) informations.

The first column names schould be : 

| #CHROM | POS | REF | ALT | %Nref | GST1 | GST2 | GST3 | F1 | F2 | F3 |
- #CHROM = chromosome name
- POS = position of DSNP
- REF = Base of the reference allele of this DSNP
- ALT = Base of the alternative allele of this DSNP
- %Nref = Percentage of maximal missing data for this DSNP
- GST = value of GST (inter-population differentiation parameter) (With 1,2,3 the ancestors names)
- F = Alternative allele frequency for each ancestor (With 1,2,3 the ancestors names)

example : <http://montreuillois.cirad.fr/acomte/gbnomique/blob/master/data/MatriceDSNPs.csv>

	-g | --gst : maximum value of the GST of a DSNP to be a marker (default = 0.9). 

	-m | --missing : maximum value of missing data frequence for a DSNP to be a marker (default = 0.3) 

### Ouputs

A matrix containing all the ancestry informative markers for every ancestors.

| ancestor | chromosome | position | allele |
 
- ancestor = Ancestor names
- chromosome = Chromosome numbers
- position = Position of the DSNP marker
- allele = Base of the ancestral allele

## TraceAncestor.pl

### Usage

	perl TraceAncestor.pl [-t matrix file] [-v vcf file] [-c color file] [-p ploidy] [-w number of markers by window] [-l LOD value] [-s threshold for LOD] [-k window size in K-bases] [-i particular hybrid to focus on] [-f focus file with several hybrids to focus on]

### Inputs

	-t | --input : reference matrice. 

The first column names schould be : 

| ancestor | chromosome | position | allele |
 
- ancestor = Ancestor names
- chromosome = Chromosome numbers
- position = Position of the DSNP marker
- allele = Base of the ancestral allele

example : <http://montreuillois.cirad.fr/acomte/gbnomique/blob/master/data/Diagnosis_matrix>

	-v | --vcf : vcf of the hybrid population

example : <http://montreuillois.cirad.fr/acomte/gbnomique/blob/master/data/Citrus.vcf>

	-c | --color : color file

Color file where the user can choose one color for each ancestor for the painting.
Indeterminations (NA) will always be grey.
If you use the ideogram output of all the hybrids, the bands between chromosoms will always be black.

example : <http://montreuillois.cirad.fr/acomte/gbnomique/blob/master/data/color>

	-p | --ploidy : ploidy of the hybrid population (2, 3 or 4).

	-w | --window : number of markers by window (default = 10)

	-l | --lod : lod value to conclude for one hypothesis or an other (default = 3)

	-s | --threshold : threshold for the calcul of LOD score (default = 0.99)

	-k | --cut : number of K bases for one window (default = 100)

	-i | --ind : particular hybrid you want to focus on (optional). If -f not indicated. 

	-f | --focus : file containing several hybrids to focus on (more than one). If -i not indicated (optional). If -i and -f are blank, focus on all hybrids.

	-h | --help : display an help

### Ouputs

- If -i is filled: 
	
	1- *ideogram_hybridname* : the painting data. An Ideogram output compatible with ideogram.js : <http://genomeharvest.southgreen.fr/visu/ideogram/newindex.php>
	
	2- *len_ideogram_hybridname* : the chromosomes data. An Ideogram output compatible with ideogram.js
	
	3- *ancestorFreq* : frequency of ancestors alleles along chromosome for the particular hybrid focused.


- If -f is filled: 

	1- *ideogram_allInd* : the painting data. An output compatible with ideogram.js. Instead of a normal ideogram, chromosomes are replaces by hybrids. Each chromosomes is put side to side to represent the whole genome.  : <http://genomeharvest.southgreen.fr/visu/ideogram/newindex.php>
	
	2- *len_allInd_ideogram* : the chromosomes data. Length of the representation of genomes of each hybrids. An output compatible with ideogram.js
	
	3- *ancestorFreq* : frequency of ancestors alleles along chromosome for all the hybrids of the focus_file.
	
	4- *circos* : the painting data. An output compatible with circos.js <http://genomeharvest.southgreen.fr/visu/circosJS/demo/index.php>

	5- *len_allInd_Circos* : the chromosomes data. An output compatible with circos.js <http://genomeharvest.southgreen.fr/visu/circosJS/demo/index.php>


- If none is filled : 

	1- *ideogram_allInd* : the painting data. An Ideogram output compatible with ideogram.js. Instead of a normal ideogram, chromosomes are replaces by hybrids. Each chromosomes is put side to side to represent the whole genome : <http://genomeharvest.southgreen.fr/visu/ideogram/newindex.php>. 
	
	2- *len_allInd_ideogram* : the chromosomes data. Length of the representation of genomes of each hybrids. An Ideogram output compatible with ideogram.js
	
	3- *ancestorFreq* : frequency of ancestors alleles along chromosome for all hybrids.

	4- *circos* : the painting data. An output compatible with circos.js <http://genomeharvest.southgreen.fr/visu/circosJS/demo/index.php>

	5- *len_allInd_Circos* : the chromosomes data. An output compatible with circos.js <http://genomeharvest.southgreen.fr/visu/circosJS/demo/index.php>

