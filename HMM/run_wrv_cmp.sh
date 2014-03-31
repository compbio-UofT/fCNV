#!/bin/bash

input_dir=$1
output_dir=$2
mix=$3

export SGE_O_PATH=$PATH

#for doc_file in $input_dir/IM1-B-chr1-150000000-150500000-duplicate.alleles_doc.txt
for x in 1
do
    doc_file=`find $input_dir/__plasma-chr$x-*alleles_doc.txt`
    tgt_file=`echo $doc_file | sed -e "s/alleles_doc/target/g"`
    #plasma_pile_file=`echo $doc_file | sed -e "s/alleles_doc/1/g"`
    plasma_pile_file=/dupa-filer/laci/I1/covChr1/i.plasma.txt
    ref_pile_file=/dupa-filer/laci/I1/covChr1/g.plasma.txt
    #ref_pile_file=/dupa-filer/laci/I1/chr$x/M.pile.txt
    seq_file=/dupa-filer/laci/I1/chr$x/chr$x.fa
    log_file=$output_dir/$(basename "$doc_file" ".alleles_doc.txt")"IpGp.out"
    
    echo $doc_file $tgt_file $log_file
    #time pypy fcnv.py $doc_file $tgt_file $plasma_pile_file $ref_pile_file $seq_file > $log_file 2>&1 &
    qsub -q all.q -R y -V -l h_vmem=25G -l h_rt=04:00:00 -wd $output_dir -o $log_file -e $log_file -b y fcnv.py $doc_file $tgt_file $plasma_pile_file $ref_pile_file $seq_file --ff $mix --useCvrg
    
    
    plasma_pile_file=/dupa-filer/laci/I1/covChr1/i.plasma.txt
    ref_pile_file=/dupa-filer/laci/I1/covChr1/i.maternal.txt
    log_file=$output_dir/$(basename "$doc_file" ".alleles_doc.txt")"IpIm.out"
    echo $doc_file $tgt_file $log_file
    qsub -q all.q -R y -V -l h_vmem=25G -l h_rt=04:00:00 -wd $output_dir -o $log_file -e $log_file -b y fcnv.py $doc_file $tgt_file $plasma_pile_file $ref_pile_file $seq_file --ff $mix --useCvrg
    
    plasma_pile_file=/dupa-filer/laci/I1/covChr1/g.plasma.txt
    ref_pile_file=/dupa-filer/laci/I1/covChr1/i.plasma.txt
    log_file=$output_dir/$(basename "$doc_file" ".alleles_doc.txt")"GpIp.out"
    echo $doc_file $tgt_file $log_file
    qsub -q all.q -R y -V -l h_vmem=25G -l h_rt=04:00:00 -wd $output_dir -o $log_file -e $log_file -b y fcnv.py $doc_file $tgt_file $plasma_pile_file $ref_pile_file $seq_file --ff $mix --useCvrg
    
    plasma_pile_file=/dupa-filer/laci/I1/covChr1/g.plasma.txt
    ref_pile_file=/dupa-filer/laci/I1/covChr1/g.maternal.txt
    log_file=$output_dir/$(basename "$doc_file" ".alleles_doc.txt")"GpGm.out"
    echo $doc_file $tgt_file $log_file
    qsub -q all.q -R y -V -l h_vmem=25G -l h_rt=04:00:00 -wd $output_dir -o $log_file -e $log_file -b y fcnv.py $doc_file $tgt_file $plasma_pile_file $ref_pile_file $seq_file --ff $mix --useCvrg
    
done
wait
