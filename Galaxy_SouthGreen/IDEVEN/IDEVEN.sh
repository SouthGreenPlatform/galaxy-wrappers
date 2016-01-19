#!/bin/bash

spe1=$1;
spe2=$2;
familyfile=$3;
output=$4;

module load system/java/jre8

if [ -d "/bank/genfam/IDEVEN/${spe1}_${spe2}/" ]; then
ls /bank/genfam/IDEVEN/${spe1}_${spe2}/*.aligncoords | \
(while read files
do
	all=$files;
done 

ls /bank/genfam/IDEVEN/${spe1}_${spe2}/*.aligncoords.qac*.40 | \
(while read ofiles
do
	ortho=$ofiles;
done 

ls /bank/genfam/IDEVEN/${spe1}_${spe2}/*.ks | grep -v 'old' | \
(while read afiles
    do
        allks=$afiles;
    done 

ls /bank/genfam/IDEVEN/${spe1}_${spe1}/*.aligncoords | \
(while read p1files
    do
        para1=$p1files;
    done 

ls /bank/genfam/IDEVEN/${spe1}_${spe1}/*.ks | grep -v 'old'  | \
(while read pk1files
    do
        para1ks=$pk1files;
    done 

ls /bank/genfam/IDEVEN/${spe2}_${spe2}/*.aligncoords | \
(while read p2files
    do
        para2=$p2files;
    done 

ls /bank/genfam/IDEVEN/${spe2}_${spe2}/*.ks | grep -v 'old' | \
(while read pk2files
    do
        para2ks=$pk2files;
    done 
	
ls /bank/genfam/IDEVEN/${spe1}_${spe2}/*.q.localdups  | \
(while read tan1files
    do
        tan1file=$tan1files;
         echo $tan1file;
    done 
	
ls /bank/genfam/IDEVEN/${spe1}_${spe2}/*.s.localdups  | \
(while read tan2files
    do
        tan2file=$tan2files;
        echo $tan2file;
    done 
	
echo java  -jar ~/galaxy/tools/SouthGreen/IDEVEN/IDEVEN.jar $allks $para1ks $para2ks $para2 $para1 $all $ortho $tan1file $tan2file $familyfile $output
# if [ -z "$tan1file" ]; then
java  -jar ~/galaxy/tools/SouthGreen/IDEVEN/IDEVEN.jar $allks $para1ks $para2ks $para2 $para1 $all $ortho $tan1file $tan2file $familyfile $output
)
)
)
)
)
)
)
)
)
elif [ -d "/bank/genfam/IDEVEN/${spe2}_${spe1}/" ]; then
ls /bank/genfam/IDEVEN/${spe2}_${spe1}/*.aligncoords | \
(while read files
do
	all=$files;
done 

ls /bank/genfam/IDEVEN/${spe2}_${spe1}/*.aligncoords.qac*.40 | \
(while read ofiles
do
	ortho=$ofiles;
done 

ls /bank/genfam/IDEVEN/${spe2}_${spe1}/*.ks | grep -v 'old' | \
(while read afiles
    do
        allks=$afiles;
    done 

ls /bank/genfam/IDEVEN/${spe1}_${spe1}/*.aligncoords | \
(while read p1files
    do
        para1=$p1files;
    done 

ls /bank/genfam/IDEVEN/${spe1}_${spe1}/*.ks | grep -v 'old' | \
(while read pk1files
    do
        para1ks=$pk1files;
    done 

ls /bank/genfam/IDEVEN/${spe2}_${spe2}/*.aligncoords | \
(while read p2files
    do
        para2=$p2files;
    done 

ls /bank/genfam/IDEVEN/${spe2}_${spe2}/*.ks | grep -v 'old' | \
(while read pk2files
    do
        para2ks=$pk2files;
    done 

	ls /bank/genfam/IDEVEN/${spe2}_${spe1}/*.q.localdups  | \
(while read tan1files
    do
        tan1file=$tan1files;
    done 
	
ls /bank/genfam/IDEVEN/${spe2}_${spe1}/*.s.localdups  | \
(while read tan2files
    do
        tan2file=$tan2files;
    done 
	
echo java  -jar ~/galaxy/tools/SouthGreen/IDEVEN/IDEVEN.jar $allks $para1ks $para2ks $para2 $para1 $all $ortho $tan1file $tan2file $familyfile $output
# if [ -z "$tan1file" ]; then  
java  -jar ~/galaxy/tools/SouthGreen/IDEVEN/IDEVEN.jar $allks $para1ks $para2ks $para2 $para1 $all $ortho $tan1file $tan2file $familyfile $output
)
)
)
)
)
)
)
)
)
fi

# elif [ -d "/bank/genfam/IDEVEN/${spe2}_${spe1}/" ]; then
#   all= "$(ls /bank/genfam/IDEVEN/${spe2}_${spe1}/*.aligncoords)";
#   ortho="$(ls /bank/genfam/IDEVEN/${spe2}_${spe1}/*.aligncoords.qac*.40)";
#   allks="$(ls /bank/genfam/IDEVEN/${spe2}_${spe1}/*.ks)"
# ls /bank/genfam/IDEVEN/${spe2}_${spe1}/*.aligncoords | \
# (while read files
# do
# 	all=$files;
# done 
# echo $all)
# 
# 
# fi
# 
# 
# if [ -d "/bank/genfam/IDEVEN/${spe2}_${spe2}/" ]; then
#   para2= "$(ls /bank/genfam/IDEVEN/${spe2}_${spe2}/*.aligncoords)";
#   para2ks="$(ls /bank/genfam/IDEVEN/${spe2}_${spe2}/*.ks)"
# fi 

# ./IDEVE.jar $allks $para1ks $para2ks $para2 $para1 $all $familyfile $ortho $output


# ls /bank/genfam/IDEVEN/MEDTR_BRADI/*.aligncoords | \
# (while read files
# do
# 	all=$files;
# done 
# ls /bank/genfam/IDEVEN/MEDTR_BRADI/*.aligncoords.qac*.40 | \
# (while read ofiles
# do
# 	ortho=$ofiles;
# done 
# ls /bank/genfam/IDEVEN/MEDTR_BRADI/*.ks | \
# (while read afiles
# do
# 	allks=$afiles;
# done 
# echo $all
# echo $ortho
# echo $allks
# )
# )
# )
