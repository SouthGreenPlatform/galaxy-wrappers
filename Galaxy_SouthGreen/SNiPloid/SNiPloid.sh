#!/bin/bash


directory=`dirname $0`

if [ -d "$HOME/galaxy/static/images" ]; then
 cp -rf $directory/img_sniploid $HOME/galaxy/static/images
fi
if [ -d "$HOME/galaxy/static/images" ]; then
 cp -rf $directory/img_sniploid $HOME/galaxy/static/images
fi


if [ "$1" = "poly" ] # poly analysis
then
	perl $directory/SNiPloid.pl --dp $2 --dp2	$3 --vp $4 --cpp $5 --vp2 $6 --cpp2 $7 --oc $8 --oh $9 --ohs ${10} --ocs ${11} --elq ${12}  --vfp1 ${14} --vfp2 ${15} --img ../../../static/images/img_sniploid/ 2>&1 >>${13} 
	if [ "${16}" = "yes" ]
	then
		perl $directory/DrawMapOfOccurences.pl -c ${11} -a ${17} -o ${18} -m 50 -t polyploid_polyploid
	fi
elif [ $2 -eq 0 ] # Ref ext
then
	perl $directory/SNiPloid.pl  --dp $3 --dg1 $4 --vp $5 --cpp $6 --vg1 $7 --cg1 $8 --dg2 $9 --vg2 $10 --cg2 $11 -- ref 0 --oc ${12} --oh ${13} --ohs ${14} --ocs ${15} --elq ${16} --img ../../../static/images/img_sniploid/ 2>&1 >>${17}
else	# Ref int == 1
	perl $directory/SNiPloid.pl --gn2 $3 --dp $4 --dg1 $5 --vp $6 --cpp $7 --vg1 $8 --cg1 $9 --ref 1 --oc ${10} --oh ${11} --ohs ${12} --ocs ${13} --elq ${14}  --img ../../../static/images/img_sniploid/ 2>&1 >>${15}
	if [ "${16}" = "yes" ]
	then
		echo ${16} >>${15}
		perl $directory/DrawMapOfOccurences.pl -c ${13} -a ${17} -o ${18} -m 50 -t polyploid_diploid 2>&1 >>${15}
	fi
fi
