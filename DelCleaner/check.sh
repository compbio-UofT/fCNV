#!/bin/bash

path=/dupa-filer/laci/I1/chr1/fcnv_data 
threshToDelete=0.0
threshToAdd=0.9

rm mfile
rm pfile
rm addfile

lenn=$1

for file in $path/* 
do 
	name=`basename $file`
	#echo $name
	type=`echo $name | awk '{split($0,a,"-"); print a[5]}'`
	if [ "$type" != "delete.target.txt" ]
	then
		continue
	fi
	begin=`echo $name | awk '{split($0,a,"-"); print a[3]}'`
	end=`echo $name | awk '{split($0,a,"-"); print a[4]}'`
	size=$(( $end - $begin ))
	if (( $size != $lenn ))
	then
		continue
	fi	
	res=`pypy delOriginCleaner.py /dupa-filer/laci/I1/chr1/trio.phase.vcf $file`
	prob=`echo $res | awk '{print $2}'`
	echo $prob
	if (( `echo "$prob > $threshToDelete" | bc` == 1 ))
	then
		echo $name
		source=`echo $res | awk '{print $1}'`
		if [ $source == "M" ]
		then
			echo M
			echo $prob" "$name >> mfile
		fi
		if [ $source == "P" ]
		then
			echo P
			echo $prob" "$name >> pfile
		fi
	fi
	
	#pypy delOriginCleaner.py /dupa-filer/laci/I1/chr1/trio.phase.vcf A-chr1-192925846-193025846-delete.alleles_doc.txt 
done 

num=`wc -l mfile | awk '{print $1}'`
if (( $num < 10 ))
then
	dif=$(( 10 - $num ))
	for i in $(seq 1 $dif)
	do
		while [ 1 ]
		do
			file=`pypy getRandomName.py $lenn`
			echo $file
			res=`pypy delOriginCleaner.py /dupa-filer/laci/I1/chr1/trio.phase.vcf $file`
			echo $res
			echo $file
			source=`echo $res | awk '{print $1}'`
			prob=`echo $res | awk '{print $2}'`
			if [ $source == "M" ] && (( `echo "$prob >= $threshToAdd" | bc` == 1 ))
			then
				echo $file >> addfile
				break
			fi
		done
	done
fi

num=`wc -l pfile | awk '{print $1}'`
if (( $num < 10 ))
then
	dif=$(( 10 - $num ))
	for i in $(seq 1 $dif)
	do
		while [ 1 ]
		do
			file=`pypy getRandomName.py $lenn`
			echo $file
			res=`pypy delOriginCleaner.py /dupa-filer/laci/I1/chr1/trio.phase.vcf $file`
			echo $res
			echo $file
			source=`echo $res | awk '{print $1}'`
			prob=`echo $res | awk '{print $2}'`
			if [ $source == "P" ] && (( `echo "$prob >= $threshToAdd" | bc` == 1 ))
			then
				echo $file >> addfile
				break
			fi
		done
	done
fi

