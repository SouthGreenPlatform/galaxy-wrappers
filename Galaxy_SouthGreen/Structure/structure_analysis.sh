#!/bin/bash

mkdir tmpdir$$

directory=`dirname $0`

MAINPARAM="tmpdir$$/mainparam.txt"
EXTRAPARAM="tmpdir$$/extraparam.txt"
ISUSERMAINPARAM=1
ISUSEREXTRAPARAM=1
TMPOUT="none"
DATAFILE=$1
OUTFILE=$2
LOGFILE=$3
kmin=$4
runs=$5
simsum_file=$6
kmax=$7
Rfile=$8
NUMINDS=0
shift; shift; shift; shift; shift; shift; shift; shift;

for param in $*
do
	if [ "none" != "$TMPOUT" ]
	then
		name=`echo $param | awk -F"," '{print $1}'`
		value=`echo $param | awk -F"," '{print $2}'`
		
		if [ -n "$value" ]
		then
			if [ $value == "true" ]
			then
				value=1
			elif [ $value == "false" ]
			then
				value=0
			fi
			
			echo "#define $name  $value" >> "${TMPOUT}"
		fi
	fi
	
	if [ "mainparam" == $param ]
	then
		TMPOUT=$MAINPARAM
		ISUSERMAINPARAM=0
	fi
	
	if [ "extraparam" == $param ]
	then
		TMPOUT=$EXTRAPARAM
		ISUSEREXTRAPARAM=0
	fi
done
#fi


if [ $ISUSERMAINPARAM -eq 0 ]
then
	MAINPARAM=$MAINPARAM
fi

if [ $ISUSERMAINPARAM -eq 1 ]
then
	MAINPARAM=$1
fi

if [ $ISUSEREXTRAPARAM -eq 0 ]
then
	EXTRAPARAM=$EXTRAPARAM
fi

if [ $ISUSEREXTRAPARAM -eq 1 ]
then
	EXTRAPARAM=${!#} #get last arg
fi

#get number of induviduals from mainparam file
NUMINDS=`more $MAINPARAM | perl -lne 'if( /^#\w+\s+NUMINDS\s+(\d+)/i ){ print "$1" }'`
let NUMINDS+=1



#cp $MAINPARAM /usr/local/bioinfo/galaxy_dev/galaxy_dist/tools/structure/newXml/test/mainparam
#cp $EXTRAPARAM /usr/local/bioinfo/galaxy_dev/galaxy_dist/tools/structure/newXml/test/extraparam
#cp $DATAFILE /usr/local/bioinfo/galaxy_dev/galaxy_dist/tools/structure/newXml/test/datafile

# Inclure le nombre de groupe
echo -e "#!/bin/bash\n" > run.sh
#echo -e "sleep 120 " >> run.sh
for K in `seq $kmin $kmax`
do
	echo -e "$directory/structure -m $MAINPARAM -e $EXTRAPARAM -i $DATAFILE -K $K -o tmpdir$$/resultat_${K}_\${SGE_TASK_ID}" >> run.sh
	echo -e "SEP $K " >> run.sh
done

qsub -t 1-$runs -tc 16 -o $LOGFILE -e $LOGFILE -N structure_galaxy -q web.q -sync y "run.sh"

ls -R > $LOGFILE
cat run.sh >> $LOGFILE

#exit;

# Creating the simsum file
if [ -e "tmpdir$$/resultat_${kmin}_1_f" ]
then
	title="File\tRun\tK\tLnP(D)\tVar\tmean_Ln_likhd"
	for fst_no in $(seq $kmin $kmax)
	do
		title=$title"\tFst_$fst_no"
	done
	
	echo -e "$title" >> "$simsum_file";
	
	for fi in `ls tmpdir$$/resultat*`
	do
		nbRun=$(expr match $fi '.*resultat_[0-9]\+_\(.*\)_f$')
		
		current_k=$(expr match $fi '.*resultat_\([0-9]\)\+_.*_f$')
		
		txt=$(grep 'Estimated Ln Prob of Data' $fi)
		lnP=$(expr match "$txt" '.*= \(.*\)$')
		
		txt=$(grep 'Mean value of ln likelihood' $fi)
		mean_ln=$(expr match "$txt" '.*= \([0-9\.\-]\+\).*$')
		
		txt=$(grep 'Variance of ln likelihood' $fi)
		varlhd=$(expr match "$txt" '.*= \([0-9\.\-]\+\).*$')
		
		txt=$(grep 'Mean value of alpha' $fi)
		mean_alpha=$(expr match "$txt" '.*= \([0-9\.\-]\+\).*$')
		
		fst=""
		for ptitK in $(seq 1 $kmax)
		do
			txt=$(grep "Mean value of Fst_$ptitK" $fi)
			if [ -z "$txt" ]; then
				fst=$fst"NA\t"
			else
				fst=$fst$(expr match "$txt" '.*= \([0-9\.\-]\+\).*$')"\t"
			fi
		done
	
		echo -e "$fi\t$nbRun\t$K\t$lnP\t$varlhd\t$mean_ln\t$fst" >> "$simsum_file";
		echo -e  "$current_k\t$fi" >> "$Rfile";
	done

	tar -zcvf $OUTFILE tmpdir$$/resultat*
fi

# AFAC
if [ -e "tmpdir$$/resultat_f" ]
then
	#mv tmpdir$$/resultat_f $OUTFILE
	`grep 'Estimated Ln Prob of Data' tmpdir$$/resultat_f > $OUTFILE`
	`grep -A $NUMINDS 'Inferred ancestry of individuals:' tmpdir$$/resultat_f >> $OUTFILE`
fi
