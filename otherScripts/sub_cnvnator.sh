#!/bin/bash
tmp_path=/tmp
if [[ -n "$TMPDIR" ]]; then
    tmp_path=$TMPDIR
fi

x=$1 #chr number

for G in M P F
do   
    echo "samtools view -b -L /dupa-filer/laci/I1/chr$x/regions.bed -o $tmp_path/$G.reg.bam /dupa-filer/laci/I1/chr$x/__$G.part.bam"
    samtools view -b -L /dupa-filer/laci/I1/chr$x/regions.bed -o $tmp_path/$G.reg.bam /dupa-filer/laci/I1/chr$x/__$G.part.bam
    samtools index $tmp_path/$G.reg.bam
    
    #logFile=/dupa-filer/laci/I1/chr$x/cnvnatorCalls/cnvnator_run.log

    cd /dupa-filer/laci/I1/chr$x/cnvnatorCalls/
    ./run_cnvnator.sh $tmp_path/$G.reg.bam $G.reg.100.root $x 100 > $G.cnvs.100.txt &
done
wait
echo "DONE."
