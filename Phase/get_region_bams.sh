#!/bin/bash

# $1 - chromosome to extract

#logfile=log_$0.$$.log
#exec > $logfile 2>&1

date

if [ $# -ne 1 ]; then
    echo "Wrong number of arguments: got $#, expected 1"
    exit
fi

bindir="/dupa-filer/laci/bin"

#parse out chromosome ID
region=$1
chr=`echo $region | awk 'BEGIN {FS=":"}{print $1}'`

echo "merging and removing PCR duplicates:"

res_dir="/dupa-filer/laci/I1/$chr"

samtools merge -h /dupa-filer/laci/I1/chr20/header.sam -R "$chr" - /dupa-filer/laci/phs000500.v1.p1/reads-src/SRP017783/SRS387404/*/*/*.bam | samtools rmdup - $res_dir/__M.allreads.part.bam &
sleep 200
samtools merge -R "$chr" - /dupa-filer/laci/phs000500.v1.p1/reads-src/SRP017783/SRS387404/SRX21977*/*/*.bam | samtools rmdup - "$res_dir/__M.part.bam" &
sleep 200
samtools merge -R "$chr" - /dupa-filer/laci/phs000500.v1.p1/reads-src/SRP017783/SRS387408/*/*/*.bam | samtools rmdup - "$res_dir/__P.part.bam" &
echo "wait1"
wait
samtools merge -R "$chr" - /dupa-filer/laci/phs000500.v1.p1/reads-src/SRP017783/SRS467990/*/*/*.bam | samtools rmdup - "$res_dir/__F.part.bam" &
sleep 200
samtools merge -R "$chr" - /dupa-filer/laci/phs000500.v1.p1/reads-src/SRP017783/SRS387409/*/*/*.bam | samtools rmdup - "$res_dir/__plasma.part.bam" &
echo "wait2"
wait
samtools index "$res_dir/__M.allreads.part.bam" &
sleep 100
samtools index "$res_dir/__M.part.bam" &
samtools index "$res_dir/__P.part.bam" &
sleep 200
samtools index "$res_dir/__F.part.bam" &
samtools index "$res_dir/__plasma.part.bam" &
wait
