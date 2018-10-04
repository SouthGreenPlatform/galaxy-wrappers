#!/usr/bin/env python3
from Bio import SeqIO
from subprocess import call, Popen, PIPE
from pybedtools import BedTool
from version import getVersion, updateTag
from timeit import Timer
import sys, os, re, tempfile, argparse, warnings, datetime, collections

parser = argparse.ArgumentParser(add_help=False)

required = parser.add_argument_group("Required arguments ")
optional = parser.add_argument_group("Optional arguments ")
required.add_argument("-f1", "--fasta1", dest="fasta1",default=None, help="Input of genome fasta1 (reference) <fasta>.")
required.add_argument("-f2", "--fasta2", dest="fasta2", default=None, help="Input of genome fasta2 (target) <fasta>.")
required.add_argument("-ti", "--tabinput", dest="tabinput", default=None, help="Input of tabbed file related to fasta1 <bed/gff/vcf>")
optional.add_argument("-b", "--flank", dest="flank", type=int, default=50, help="Size of flank region to extract from each side of the annotation (default : 50).")
optional.add_argument("-c", "--cds",dest="cds", action="store_true", help="Enable control of postions of CDS inside mRNA (or gene).")
optional.add_argument("-d", "--directory", dest="directory", default="result", help="Name of the directory where files are generated (default : result).")
optional.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS,help="Show this help message then exit.")
optional.add_argument("-i", "--index",dest="index", default=None, action="count", help="Create the index of fasta2 if it doesn't exist (-ii for forcing).")
optional.add_argument("-l", "--loss",dest="loss", action="store_true", help="Enable the creation of 'StatsLoss' file which contain percentage of loss.")
optional.add_argument("-n", "--notempfile",dest="notempf", action="store_true", help="Create all file instead of using temporary file.")
optional.add_argument("-o", "--output", dest="out", default=None, help="Output file (same format of tabbed file input).")
optional.add_argument("-t", "--type",dest="typeA", default=None, help="Only extract annotation of specified type (gff file only) (example : \"mRNA exon cds\" (case insensitive)).")
optional.add_argument("-u", "--update", dest="update", action="store_true", help="Update the .version file by checking the git log")
optional.add_argument("-v", "--verbose",dest="verbose", default=1, type=int, choices=[0,1,2], help="Change verbosity level.")
optional.add_argument("-ver", "--version",dest="version", action="store_true", help="Show version and date of last update then exit.")
optional.add_argument("-w", "--warning",dest="warn", action="store_true", help="Disable warnings.")
if __name__ == '__main__':
    """
        Some configuration for the argument parser like
        checking require argument : --fasta1, --fasta2
        and --tabinput. We setup the output file name,
        create output directory (default : 'result/'), setup the
        tempfile or file depending --notempfile option
        and signal useless argument (e.g. --type if the
        format of tabinput isn't GFF).
    """
    args = parser.parse_args()

    if len(sys.argv)<2:
        getVersion()
        parser.print_help()
        sys.exit(1)
    if args.update :
        updateTag()
        sys.exit(1)
    if args.version :
        getVersion()
        sys.exit(1)
    if args.warn :
        warnings.filterwarnings("ignore")

    if args.fasta1 == None:
        sys.exit("ERROR : Argument --fasta1 (-f1) is missing.")
    if args.fasta2 == None:
        sys.exit("ERROR : Argument --fasta2 (-f2) is missing.")
    if args.tabinput == None:
        sys.exit("ERROR : Argument --tabinput (-ti) is missing.")

    global ext
    ext=BedTool(args.tabinput).file_type
    if ext == "gff" :
        ext ="gff3"
    if args.directory[-1]=="/" :
        args.directory=args.directory[:-1]
    try:
        os.mkdir(args.directory)
        if args.verbose !=0 :
            print("\n ----- Creating directory '"+args.directory+"/'. -----")
    except :
        pass

    # Creating tempfile (or file if args --notempfile is enable)
    tabOut=tempfile.NamedTemporaryFile()
    tabO=tabOut.name
    if not args.notempf :
        tabOPut=tempfile.NamedTemporaryFile()
        tabOP=tabOPut.name
        aln=tempfile.NamedTemporaryFile()
        alnN=aln.name
        selectS = tempfile.NamedTemporaryFile()
        selectedSeq = selectS.name
    else :
        tab=".".join((args.tabinput).split(".")[:-1])
        tabOP=args.directory+"/"+tab.split("/")[-1]+"_out."+ext
        alnN=args.directory+"/aln_out_"+ext+".sam"
        fasta1Out=".".join((args.fasta1).split(".")[:-1])
        selectedSeq=args.directory+"/"+fasta1Out.split("/")[-1]+"_selected_"+ext+".fasta"

    if args.typeA != None and ext != "gff3" :
        print("")
        warnings.warn("Argument --typeA is ignored because --tabinput format isn't GFF.",Warning)

    if args.cds and ext != "gff3" :
        print("")
        warnings.warn("Argument --cds is ignored because --tabinput format isn't GFF.",Warning)

    if args.typeA != None and ext == "gff3" :
        typ=(re.findall("[a-zA-Z0-9]+",args.typeA))
        typeAclean=[l.lower() for l in typ]
    else :
        typeAclean=""
    if typeAclean ==[] :
        sys.exit("ERROR : Argument --type (-t) is empty.")
    for i in typeAclean :
        if i not in ["cds","mrna","gene"] and args.cds :
            print("")
            warnings.warn("Argument --cds is ignored because there is unsupported type in --type argument when -c is passed.",Warning)
            args.cds=None
            break
    if ext == "gff3" and args.cds and ("gene" not in typeAclean and "mrna" not in typeAclean or "cds" not in typeAclean) :
        print("")
        warnings.warn("Argument --cds is ignored because there is missing type in --type argument (required : 'CDS' and mRNA or gene).",Warning)
 
    res2 = re.search("/([^/]+)(\.\w+)",args.fasta2)
    if args.out == None:
        args.out = args.directory+"/"+res2.group(1)+"_out."+ext
    else :
        args.out = args.directory+"/"+args.out.split(".")[0]+"."+ext
        if os.path.isfile(args.out) :
            print("")
            warnings.warn("Argument --output is ignored because file '"+args.out+"' already exist.",Warning)
            args.out = args.directory+"/"+res2.group(1)+"_out."+ext

def checkDependency() :
    """
        Check if all third-party program are installed by
        making a system call with "module list" then
        comparing with a required module list. If at least
        one thrird-party program isn't installed, it quit.
        Currently, this function check for environment
        modules but it need further test or a more generic
        test.
    """
    p=Popen(["module list"],shell=True, stderr=PIPE)
    outs,err=p.communicate()
    listM=re.findall("\w+/\w+/(?:\d+\.?)+",str(err))
#    mod=tempfile.NamedTemporaryFile()
#    listMod=mod.name
#    with open(listMod,"r+") as out :
#        call(["module list"],shell=True,stderr=out)
#        out.seek(0)
#        for line in out :
#            res=re.findall("\)\s([^\s]+)\s+",line)
#            if "listM" in locals() :
#                for i in res :
#                    listM.append(i)
#            else :
#                listM=res
    mandatoryMod=["bioinfo/bwa/0.7.15","bioinfo/bedtools/2.24.0"]
    goInstall=""
    for i in mandatoryMod:
        if i not in listM :
            goInstall += ("/".join(i.split("/")[1:]))+"  "
    if goInstall :
        sys.exit("ERROR ! Please, install following tools : " + goInstall)
    return

def parseFa() :
    """
        Later, in order to get the flank region and not
        going over last position of the chromosome, we need
        to index the length of them.
        This function parse the inputted file fasta1 and
        write in a tempfile all the chromosome ID and length.
    """
    lenChr = tempfile.NamedTemporaryFile()
    with open(lenChr.name,"w") as lenC :
        for seqF in SeqIO.parse(args.fasta1,"fasta") :
            lenC.write(seqF.id+"\t"+str(len(seqF))+"\n")
    return (lenChr)

def prefix() :
    """
        In order to link the tabbed file and the fasta
        file, we need to get the same chromosome ID. This
        function check for this and change ID prefix in the
        tabbed file to the same as the ID prefix in the
        fasta file.
    """
    global change
    change=False
    fileTab=BedTool(args.tabinput)
    fasta=SeqIO.parse(args.fasta1,"fasta")
    first_seq=next(fasta)
    if (fileTab[0][0][0:-1])!= (first_seq.id[0:-1]) :
        print("\n ----- Changing prefix of chromosome ID in '"+ args.tabinput+"' to match with chromosome ID in '"+args.fasta1+"' ----- ")
        s=""
        with open(args.tabinput) as tabHash :
            for hasht in tabHash :
                if hasht[0] == "#" :
                    s+=hasht
                else :
                    break
            for feat in fileTab :
                numc=re.search("[^\d]?(\d+)$",feat[0])
                if numc :
                    s+= first_seq.id[0:-1]+numc.group(1)+"\s"+("\s".join(feat[1:])+"\n")
        BedTool(s, from_string=True, deli="\s").saveas(tabO)
        change=True
    return

def cutGff() :
    """
        If the argument --typeA is specified, we need to cut
        the gff in order to keep only specified type. We use
        lower() function on each type in order to be case
        insensitive.
    """
    global typeAclean
    if args.verbose != 0:
        print("\n ----- Extracting GFF's annotation which match '"+",".join(typeAclean)+"'. -----")
    with open(args.tabinput,"r") as inpTab, open(tabO,"w") as out:
        for line in inpTab :
            fields = line.strip().split("\t")
            if fields[2].lower() in typeAclean :
                out.write(("\t".join(fields))+"\n")
    return

def getFlank() :
    """
        Change the position in tabbed file in order to
        include flank region. For vcf file, we need to
        generate a new temporary tabbed file which contain
        the ID, start position and stop position. If the
        argument --notempfile is passed, it create the file
        [tabbed_file_name]_out.[ext] in the directory
        ./result.
    """
    lenChr=parseFa()
    if (args.typeA != None and ext== "gff3") or change :
        fileTab2=BedTool(tabO)
    else :
        fileTab2=BedTool(args.tabinput)
    if args.verbose!=0 and args.notempf:
        print("\n ----- Creating file '"+tabOP +"'. ----- ")
    if (fileTab2.file_type == "bed"):
        with open(lenChr.name,"r") as lenC :
            res=""
            for feature in fileTab2 :
                fjoin=("__".join(str(feature).split("\t"))).replace(" ","\\s")
                for line in lenC :
                    lenghtC=re.search(feature.chrom+"\t(\d+)",line)
                    if lenghtC :
                        break
                lenC.seek(0)
                if feature.stop+args.flank > int(lenghtC.group(1)): # TODO : unreadable, gotta change that
                    stop=int(lenghtC.group(1))
                else :
                    stop=feature.stop+args.flank
                if feature.start-args.flank < 0 :
                    start=0
                else :
                    start=feature.start-args.flank
                res += feature.chrom +"\t"+str(start)+"\t"+ str(stop)+"\t"+fjoin+"\n"
        BedTool(res, from_string=True, deli="\t").saveas(tabOP)
        #fileTab2.slop(b=args.flank, g=lenChr.name ,output=tabOP, header=True)
    if (fileTab2.file_type == "gff"):
        with open(lenChr.name,"r") as lenC :
            res=""
            countCds=0
            for feature in fileTab2 :
                taline=str(feature).split("\t")
                if taline[2]=="CDS" :
                    countCds+=1
                    taline[-1]=taline[-1][0:-1]+";Note="+str(countCds)
                else :
                    countCds=0
                fjoin=("__".join(taline)).replace(" ","\\s")
                for line in lenC :
                    lenghtC=re.search(feature.chrom+"\t(\d+)",line)
                    if lenghtC :
                        break
                lenC.seek(0)
                if feature.stop+args.flank > int(lenghtC.group(1)): # TODO : unreadable, gotta change that
                    stop=int(lenghtC.group(1))
                else :
                    stop=feature.stop+args.flank
                if feature.start-args.flank < 0 :
                    start=0
                else :
                    start=feature.start+1-args.flank
                res += feature.chrom +"\t"+str(start)+"\t"+ str(stop)+"\t"+fjoin+"\n"
        BedTool(res, from_string=True, deli="\t").saveas(tabOP)
        #fileTab2.slop(b=args.flank, g=lenChr.name ,output=tabOP, header=True)
    elif (fileTab2.file_type == "vcf") :
        with open(lenChr.name,"r") as lenC :
            res=""
            for feature in fileTab2 :
                fjoin=("__".join(str(feature).split("\t"))).replace(" ","\\s")
                for line in lenC :
                    lenghtC=re.search(feature.chrom+"\t(\d+)",line)
                    if lenghtC :
                        break
                lenC.seek(0)
                if feature.stop+args.flank-1+(len(feature[3])-1) > int(lenghtC.group(1)): # TODO : unreadable, gotta change that
                    stop=int(lenghtC.group(1))
                else :
                    stop=feature.stop+args.flank-1+(len(feature[3])-1)
                if feature.start-args.flank-1 < 0 :
                    start=0
                else :
                    start=feature.start-args.flank-1
                res += feature.chrom +"\t"+str(start)+"\t"+ str(stop)+"\t"+fjoin+"\n"
        BedTool(res, from_string=True, deli="\t").saveas(tabOP)
    return (tabOP)

def getFasta(tabOP):
    """
        This function just extract the sequence defined in
        the tabbed file create by getFlank(). The output is
        a fasta file. If --notempfile is enable, it create the
        file [fasta1_name]_selected_[ext].fasta in the
        directory ./result.
    """
    if args.verbose!=0 :
        print("\n ----- Extracting flanking sequences in '"+args.fasta1+"'. -----")
    if not args.notempf :
        #BedTool(tabOP.name).sequence(fi=args.fasta1).save_seqs(selectedSeq)
        BedTool(tabOP).sequence(name=True,fi=args.fasta1).save_seqs(selectedSeq)
    else :
        if args.verbose !=0 :
            print("\n ----- Creating file '"+selectedSeq+"'. ----- ")
        #BedTool(tabOP).sequence(fi=args.fasta1).save_seqs(selectedSeq)
        BedTool(tabOP).sequence(name=True,fi=args.fasta1).save_seqs(selectedSeq)
    return

def index() :
    """
        This function generate the index of inputted fasta2
        file with a system call to bwa index.
        Only called if call of bwa mem failed or if argument
        -ii is specified.
    """
    if args.verbose !=0:
        print("\n ----- Generating index of '"+args.fasta2+"' using 'bwa index'. -----")
    #TODO : remove system call (possible ? can't find any wrapper)
    if args.verbose ==2:
        print("")
        call(["bwa","index",args.fasta2])
    else :
        p=Popen(["bwa","index",args.fasta2], stderr=PIPE)
        err=p.communicate()
    return

def align(tabOp):
    """
        This function make a system call to bwa mem to align
        fasta2 file with selected sequence of fasta1. The
        output file is in sam format.
        If argument --notempfile is enable it create the
        file "aln_out_[ext].sam".
    """
    getFasta(tabOp)
    if args.notempf and args.verbose !=0 :
          print("\n ----- Creating file '"+alnN+"'. ----- ")
    with open(alnN,"w") as alout:
        while True :
            if args.verbose !=0:
                print("\n ----- Mapping selected sequence to fasta2 using 'bwa mem'. ----- \n")
            p=Popen(["bwa","mem", args.fasta2, selectedSeq],stdout=alout, stderr=PIPE)
            err=p.communicate()
            if err[-1].decode("utf-8")[1]=="E" :
                if args.index == None :
                    print("")
                    sys.exit("ERROR : Fail to locate the index files.")
                if args.index==1 :
                    index()
                else :
                    sys.exit(1)
            else :
                break
        if args.verbose==2:
            print(err[-1].decode("utf-8"))
    return

def parseCigar(sam) :
    """
        This function parse the CIGAR code in alignment file
        in order to get the length of the match.
    """
    lenCig=[]
    for i in sam :
        if int(i[1]) != 0: # Ignoring complementary match (flag 2048)
            continue
        leng=0
        res =re.findall("\d+\w",i[5])
        for i in res :
            if i[-1] in ["M","=","X","I","S"] :
                leng+=int(i[:-1])
        lenCig.append(leng)
    return lenCig

#def sorting(aln) :
#    fil=[]
#    if args.verbose != 0 :
#        print(" ----- Sorting file '"+aln+"'. ----- ")
#    with open(aln, "r") as samsoul, open(sortSam,"w") as soulsam :
#        for i in samsoul :
#            if i[0]!="@" :
#                fil.append(i.split("\t"))
#            else :
#                soulsam.write(i)
#        #for i in sorted(samsoul) :
#        #    soulsam.write(i)
#        for i in sorted(fil, key=lambda ok: (ok[0][0],int(ok[0].split("_")[1].split(":")[0]),int(ok[0].split(":")[1].split("-")[0]),int(ok[0].split(":")[1].split("-")[1]))) :
#            soulsam.write("\t".join(i))
#        alnN=sortSam

def samToTab() :
    """
        In order to get the result of alignment usable for
        further  use, we create a new tabbed file in the
        same format as the inputted tabbed file. This file
        is the final output of this program.
        First of all, we extract the start and stop position
        and remove the flank region to get the initial
        annotation length. Then we generate the output file
        with these new position and keep all the information
        of the original tabbed file.
        The name of output file is by defaut
        "[fasta2_name]_out.[ext]" or the name that user
        specify with argument --output (but still with [ext]
        as extension). It return the name of the file created.
    """
    tabou=""
    if ext !="vcf" :
        lengh=parseCigar(BedTool(alnN))
    if (ext=="gff3" and args.typeA != None) or change:
        args.tabinput=tabO
    with open(args.tabinput,"r") as tabi, open(alnN,"r") as sam :
        for i in tabi : # tabi = tabbed file after all modification 
            line=i.split("\t")
            if line[0][0] == "#" :
                tabou += i
                continue
            tabou += "# File generated the "+datetime.datetime.now().strftime("%d %b %Y") + " with following command line : \n"+"# "+" ".join(sys.argv)+"\n"
            break
        for f in sam : # samf = alignment file inside Bedtools object
            samf=f.split("\t")
            if samf[0][0]=="@" :
                continue
            if int(samf[1])!=0 : # Ignoring complementary match (flag 2048)
                continue
            if ext != "vcf" :
                leng=0
                res=re.findall("\d+\w",samf[5])
                for i in res :
                    if i[-1] in ["M","=","X","I","S"] :
                        leng+=int(i[:-1])
            tab=samf[0].replace("\\s"," ").split("__")
            if ext == "gff3" :
                tab[3]=int(samf[3])+args.flank
                tab[4]=int(tab[3])+leng-args.flank
                tabou+="\s".join(map(str,tab))+"\n"
            elif ext == "bed" :
                tab[1]=int(samf[3])+args.flank
                tab[2]=int(tab[1])+leng-args.flank
                tabou+="\s".join(map(str,tab))+"\n"
            elif ext == "vcf" and samf[5]==str(len(tab[3])+(args.flank*2))+"M" :
                tab[1]=int(samf[3])+args.flank
                tabou+="\s".join(map(str,tab))+"\n"
            #if f[11][-1]!="0" and f[5]=="101M":     # show ID of sequence which contain a missmatch
            #    print(f[0].split(":")[1]+"\t"+f[12]) 
            #if ext == "vcf" and samf[5]=="101M" and samf[5]==samf[12].split(":")[-1]+"M" and okw : # perfect match only for snp
            #if ext == "vcf" and samf[5]==str(len(samf[0].split("__")[3])+(args.flank*2))+"M" and okw :
        if args.verbose!=0 :
            print(" ----- Creating file '"+args.out+"'. ----- \n")
        BedTool(tabou,from_string=True, deli="\s").saveas(args.out)
    return (args.out)

def getPosCds(tab) :
    """
        This function extract position of mRNA (by defaut,
        else it's position of gene) and the relative
        position of CDS within.
        It take a file in gff format and return a
        dictionary with tuple of position of mRNA (or gene)
        and his number in keys and list of list of position
        of CDS in value.
    """
    dicoPos={}
    posGene=()
    global typeAclean
    if args.typeA == None :
        typeD = True
    else :
        if ("gene" in typeAclean and "mrna" in typeAclean) :
            typeD = True
        else :
            typeD=False
    with open(tab,"r") as out :
        numGene=0
        for line in out :
            if line[0]=="#" :
                continue
            lineSplit=line.split("\t")
            typeA = lineSplit[2].lower()
            if typeD and typeA == "gene" :
                continue
            start = lineSplit[3]
            stop = lineSplit[4]
            if typeA == "gene" or typeA == "mrna" :
                numGene+=1
                getags=lineSplit[-1]
                posGene=(getags,numGene,int(start),int(stop))
                if posGene not in dicoPos.keys():
                    dicoPos[posGene]=[]
            if typeA == "cds":
                cdsStart=int(start)-int(posGene[1])
                cdsStop=int(stop)-int(posGene[1])
                if cdsStart > cdsStop :
                    warnings.resetwarnings()
                    warnings.filterwarnings("error")
                    warnings.warn("Start > stop",Warning)
                resTag=re.search("ID=(\w+((\.?\w+)?)+)",getags)
                resTagCds=re.search("Parent=(\w+((\.?\w+)?)+)",lineSplit[-1])
                if resTag and resTagCds:
                    if str(resTag.group(1))==str(resTagCds.group(1)) :
                        dicoPos[posGene].append([cdsStart,cdsStop])
    return dicoPos

def isComplete(samtotabOut) :
    """
        This function call the function getPosCds(tab) with
        both inputted tabbed file and newly generated tabbed
        file and check if CDS in the newly generated file
        are in the same position within the mRNA (or gene).
        It take the name of the tabbed file newly generated
        in argument.
    """
    if ext == "gff3" :
        dicoPos1=getPosCds(args.tabinput)
        dicoPos2=getPosCds(samtotabOut)
        outTab = samtotabOut.split("/")[-1]
        geneInt=[]
        #lastG=0
        geneOk=0
        ok=0
        countG=0
        selectable=False
        filtered = "# File generated the "+datetime.datetime.now().strftime("%d %b %Y") + " with following command line : \n"+"# "+" ".join(sys.argv)+"\n"
        for key1 in dicoPos1.keys() :
            for key2 in dicoPos2.keys() :
                if key2[0]==key1[0] :
                    if len(dicoPos1[key1]) == len(dicoPos2[key2]) :
                        geneInt.append(key1[1])
                       # for v in range (0,len(dicoPos1[key1])) :
                       #     print(dicoPos1[key1][v])
                       #     if dicoPos1[key1][v] == dicoPos2[key2][v] : # TODO : c'est de la merde.
                       #         geneOk+=1
                       # print(geneOk)
                       # if geneOk >= len(dicoPos1[key1]) : # here we can add/rm condition to accept or not the mRNA/gene
                       #    # add the mRNA/gene number to the list of "acceptable mRNA/gene to select"
                       #     geneOk = 0
                       # else :
                       #     geneOk = 0
        if "gene" in typeAclean :
            typeC="gene"
        elif "mrna" in typeAclean :
            typeC="mrna"
        with open(samtotabOut,"r") as tabou :
            for line in tabou :
                if line[0]=="#":
                    continue
                lineS=line.strip().split("\t")
                if lineS[2].lower() == typeC : # TODO : unreadable
                    resTag=re.search("ID=(\w+((\.?\w+)?)+)",lineS[-1])
                    if resTag :
                        geneId=resTag.group(1)
                if lineS[2] =="CDS" :
                    resTagCds=re.search("Parent=(\w+((\.?\w+)?)+)",lineS[-1])
                    if resTagCds :
                        cdsId=resTagCds.group(1)
                if lineS[2].lower() ==typeC :
                    countG+=1
                    if countG in geneInt :
                        selectable=True
                    else :
                        selectable=False
                if lineS[2] =="CDS" and geneId!=cdsId :
                    selectable=False
                if selectable :
                    filtered+=("\s".join(lineS))+"\n"
            countG=0
        if args.verbose != 0 :
            print(" ----- Generating filtered GFF file '"+(args.directory+"/filtered_"+outTab)+"'. -----\n")
        BedTool(filtered, from_string=True, deli="\s").saveas(args.directory+"/filtered_"+outTab)
    return

def indexTab(fil) :
    dico={}
    with open(fil,"r") as ok :
        for i in ok :
            tab=i.split("\t")
            if tab[0][0] == "#" :
                continue
            if tab[0] not in dico.keys() :
                dico[tab[0]]=1
            else :
                dico[tab[0]]+=1
    return(dico)

def checkLoss() :
    dicOri=collections.OrderedDict(sorted(indexTab(args.tabinput).items(), key=lambda jm: (jm[0].split("_")[0],int(jm[0].split("_")[-1]))))
    dicNew=collections.OrderedDict(sorted(indexTab(args.out).items(), key=lambda jm: (jm[0].split("_")[0],int(jm[0].split("_")[-1]))))
    statFile=args.directory+"/"+".".join(args.fasta1.split(".")[:-1]).split("/")[-1]+"_statsLoss"
    with open(statFile,"w") as lossOut :
        print(" ----- Creating file '"+statFile+"'. ----- \n")
        lossOut.write("#ID\tREF\tTARGET\tLOSS\t%\n")
        for key in dicOri.keys() :
            if key in dicNew.keys():
                loss=(int(dicOri[key])-int(dicNew[key]))
                lossOut.write(key+"\t"+str(dicOri[key])+"\t"+str(dicNew[key])+"\t"+str(loss)+"\t"+str("%.2f" % round(loss*100/dicOri[key],2))+"%"+"\n")
    return

if __name__ == "__main__":
    tm = datetime.datetime.now()
    checkDependency()
    if args.typeA != None and ext == "gff3" :
        cutGff()
    prefix()
    tabOp = getFlank()
    if args.index==2:
        index()
    align(tabOp)
    #sorting(alnN)
    samtotabOut=samToTab()
    if (args.cds and ext=="gff3" and (args.typeA==None or (("gene" in typeAclean or "mrna" in typeAclean) and "cds" in typeAclean))) :
        isComplete(samtotabOut)
    if args.loss :
        checkLoss()
