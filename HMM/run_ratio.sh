#!/bin/bash

input_dir=$1
output_dir=$2

export SGE_O_PATH=$PATH

for doc_file in $input_dir/*.alleles_doc.txt
do
    tgt_file=`echo $doc_file | sed -e "s/alleles_doc/target/g"`
    log_file=$output_dir/$(basename "$doc_file" ".alleles_doc.txt")".out"
    
    #echo $doc_file $tgt_file $log_file
    #time pypy fcnv.py $doc_file $tgt_file > $log_file 2>&1 &
    qsub -q all.q -R y -V -l h_vmem=5G -l h_rt=01:00:00 -wd $output_dir -o $log_file -e $log_file -b y ratioFCNV.py $doc_file $tgt_file
done
wait
