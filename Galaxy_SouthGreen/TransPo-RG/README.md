# TransPo-RG

----------
TransPo-RG, *Transfer of Position to Resequenced Genome* (or TPRG), is a python
script for transfering position of annotation from one genome assembly to a another
(e.g. newer version for a reference genome). It extract sequences with flanking
nucleotide (50bp in each side by default) and map them on the target genome.
Then we generate the output file in the same format of the entry file with only
changing positions. 

### Usage :

It's a command-line program. It needs three files in entry : two genome file in
FASTA format and one tabbed file containing position of annotation. The tabbed
file can be VCF, GFF3 or BED. 

Example of a simple command for a reference gemone fasta1, a target genome fasta2 and a
VCF file v1.vcf :

```
python transpo-rg.py -f1 fasta1.fasta -f2 fasta2.fasta -ti v1.vcf
```

The output file will only be a VCF file named v1\_out.vcf. It will go in directory
named result/ in your current directory. 
Both name are configurable with option so you can redirect output file in whatever
directory you want.
If you want to keep all the intermediate file like alignment file (in SAM), you can
use the option -n (--notempfile).


### Requires :  

* bwa mem (0.7.15)
* bedtools (2.24.0)
* lib python :
    * biopython (1.65+) 
    * bedtools (2.24.0)
    * pybedtools (0.7.10)

----------
```
Version : 0.6.2
Last update : 13 Aug 2018

usage: transpo-rg.py [-f1 FASTA1] [-f2 FASTA2] [-ti TABINPUT] [-b FLANK] [-c]
                     [-d DIRECTORY] [-h] [-i] [-l] [-n] [-o OUT] [-t TYPEA]
                     [-u] [-v {0,1,2}] [-ver] [-w]

Required arguments :
  -f1 FASTA1, --fasta1 FASTA1
                        Input of genome fasta1 (reference) <fasta>.
  -f2 FASTA2, --fasta2 FASTA2
                        Input of genome fasta2 (target) <fasta>.
  -ti TABINPUT, --tabinput TABINPUT
                        Input of tabbed file related to fasta1 <bed/gff/vcf>

Optional arguments :
  -b FLANK, --flank FLANK
                        Size of flank region to extract from each side of the
                        annotation (default : 50).
  -c, --cds             Enable control of postions of CDS inside mRNA (or
                        gene).
  -d DIRECTORY, --directory DIRECTORY
                        Name of the directory where files are generated
                        (default : result).
  -h, --help            Show this help message then exit.
  -i, --index           Create the index of fasta2 if it doesn't exist (-ii
                        for forcing).
  -l, --loss            Enable the creation of 'StatsLoss' file which contain
                        percentage of loss.
  -n, --notempfile      Create all file instead of using temporary file.
  -o OUT, --output OUT  Output file (same format of tabbed file input).
  -t TYPEA, --type TYPEA
                        Only extract annotation of specified type (gff file
                        only) (example : "mRNA exon cds" (case insensitive)).
  -u, --update          Update the .version file by checking the git log
  -v {0,1,2}, --verbose {0,1,2}
                        Change verbosity level.
  -ver, --version       Show version and date of last update then exit.
  -w, --warning         Disable warnings.
```
