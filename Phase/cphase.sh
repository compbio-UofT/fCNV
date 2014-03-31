#!/bin/bash

if [ $# -ne 3 ]; then
    echo "Wrong number of arguments: got $#, expected 3"
    exit
fi

for (( x=$1; x<=$2; x++ ))
do
    echo chr$x
    qsub -q all.q -R y -V -l h_vmem=28G -wd /dupa-filer/laci/I1/chr$x -l h_rt=30:00:00 -o log_phaseTrio_chr$x.txt -e log_phaseTrio_chr$x.txt -S /bin/bash /dupa-filer/laci/bin/phase_trio.sh M P F /filer/hg19/hg19.fa chr$x;
    sleep $3;
    qsub -q all.q -R y -V -l h_vmem=28G -wd /dupa-filer/laci/I1/chr$x -l h_rt=30:00:00 -o log_phaseMP_chr$x.txt -e log_phaseMP_chr$x.txt -S /bin/bash /dupa-filer/laci/bin/phase_MP.sh M P /filer/hg19/hg19.fa chr$x;
    sleep $3;
done
