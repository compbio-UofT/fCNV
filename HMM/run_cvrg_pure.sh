#!/bin/bash

input_dir=$1
output_dir=$2
mix=$3
win_size=$4

export SGE_O_PATH=$PATH

#for doc_file in $input_dir/IM1-B-chr1-150000000-150500000-duplicate.alleles_doc.txt
for x in {1..22}
do   
    chr=chr$x
    doc_file=`find $input_dir/__plasma-chr$x-*alleles_doc.txt`
    tgt_file=`echo $doc_file | sed -e "s/alleles_doc/target/g"`
    plasma_pile_file=`echo $doc_file | sed -e "s/alleles_doc/1/g"`
    ref_pile_file=/dupa-filer/laci/I1/$chr/g.plasma.txt
    #ref_pile_file=/dupa-filer/laci/I1/$chr/M.pile.txt
    seq_fasta=/dupa-filer/laci/I1/$chr/$chr.fa
    log_file=$output_dir/$(basename "$doc_file" ".alleles_doc.txt")".cvrg.out"
    
    #echo $doc_file $tgt_file $log_file
    #time pypy fcnv.py $tgt_file $plasma_pile_file $ref_pile_file $seq_fasta > $log_file 2>&1 &
    qsub -q all.q -R y -V -l h_vmem=25G -l h_rt=02:30:00 -wd $output_dir -o $log_file -e $log_file -b y cvrgONLY.py $tgt_file $plasma_pile_file $ref_pile_file $seq_fasta $win_size --ff $mix
    
    sleep 2
done
wait
