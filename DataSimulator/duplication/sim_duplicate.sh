#!/bin/bash

# 1-chromosome; 2-CNV start position; 3-CNV length; 4-origin(M/P); 5-haplotype(A/B); 6-path to data;
# 7-path where to store the simulated plasma BAM; 8-path to dir with scripts

data_path=$6
results_path=$7
exec_path=$8
phase_sites=$data_path/trio.phase.vcf
bamfile=$data_path/__$4.part.bam
local_bam=/tmp/fcnv/__$4.part.bam
plasmaFile=$data_path/__plasma.part.bam
local_plasma=/tmp/fcnv/__plasma.part.bam
readLength=100
plasmaFetusRate=0.13
tmp_path=/tmp
if [[ -n "$TMPDIR" ]]; then
    tmp_path=$TMPDIR
fi
if [[ -r "$local_plasma" ]]; then
    plasmaFile=$local_plasma
    echo $plasmaFile
fi
if [[ -r "$local_bam" ]]; then
    bamfile=$local_bam
fi

chromosome=$1
begin=$2
end=$(($2 + $3))
haplotype=$5
source='I'$4'1'

region=$chromosome':'$begin'-'$end

pid=$$
logfile=$results_path/log_simDup.$pid.log
exec > $logfile 2>&1

#temp files names
tmp_plasma=$tmp_path/plasma$pid.bam
tmp_pileup_file=$tmp_path/plasma_pileup$pid.txt
raw_reads_file=$tmp_path/raw_reads$pid.sam
inside_sam_file=$tmp_path/inside$pid.sam
inside_bam_file=$tmp_path/inside$pid.bam
filtered_res_file=$tmp_path/filteredResults$pid.sam
plasma_file_prefix=$tmp_path/$source-$haplotype-$region-duplicate


echo "SimDuplication: $region $haplotype $source $3"
date

#check centromiers for begin 
wbegin=$(($begin - 20000000))
wend=$(($end + 20000000))
if ((wbegin <= 2300000))
then   
    ((wend=$wend+2300000-$wbegin))
    wbegin=2300000
fi

if ((120600000 <= wbegin && wbegin <= 147000000))
then   
    ((wend=$wend+147000000-$wbegin))
    wbegin=147000000
fi

#check centromiers for wend
if ((wend >= 243700000))
then   
    (( wbegin= $wbegin - ($wend-243700000) ))
    wend=243700000
fi

if ((120600000 <= wend && wend <= 147000000))
then   
    (( wbegin= $wbegin - ($wend-120600000) ))
    wend=120600000
fi

window=chr1:$wbegin-$wend
echo "window: $window"

samtools view -b $plasmaFile $window -o $tmp_plasma
samtools index $tmp_plasma

echo "Pileuping plasma reads in the region..."
tmp_region=$chromosome':'$(($begin - 1000))'-'$(($end + 1000))
samtools mpileup $tmp_plasma -q10 -Q10 -r $tmp_region | awk '{print $2, $4}' > $tmp_pileup_file

echo "Copying all the reads in the region..."
samtools view $bamfile $region -o $raw_reads_file
# Copying the header
samtools view -H $bamfile $region -o $inside_sam_file

echo "Filtering for the correct haplotype..."
pypy $exec_path/duplication.py $raw_reads_file $phase_sites $haplotype > $filtered_res_file

echo "Down sampling the results..."
numReads=`wc -l $filtered_res_file | awk '{print $1}'`
pypy $exec_path/duplication_down_sampler.py $filtered_res_file $numReads $readLength $tmp_pileup_file $plasmaFetusRate $region >> $inside_sam_file
samtools view -S -b $inside_sam_file -o $inside_bam_file

echo "Merging"
# Outside has the headers
samtools merge -f $plasma_file_prefix.sort.bam $tmp_plasma $inside_bam_file
#echo "Sorting"
#samtools sort $plasma_file_prefix.bam $plasma_file_prefix.sort
echo "Indexing"
samtools index $plasma_file_prefix.sort.bam

#move to results_path
mv $plasma_file_prefix.sort.bam* $results_path/

#rm $plasma_file_prefix.bam
rm $tmp_pileup_file $raw_reads_file $filtered_res_file $inside_sam_file $inside_bam_file $tmp_plasma*

echo "sim_duplicate script - DONE."

