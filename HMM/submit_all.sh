#!/bin/bash
for i in 1 2 3 4 5 6 7 8 9 10
do
  qsub -q all.q -l h_vmem=25G -S /bin/bash ~/FetalCNV/fcnv/run.sh ~/FetalCNV/fcnv/out_$i.txt $i
done
