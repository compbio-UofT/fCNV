#!/bin/bash

export SGE_O_PATH=$PATH

targetRate=$1

job_count=0
for x in {1..21}
do   
    if [[ "$job_count" == 10 ]] ; then
        job_count=0
        sleep 500
    fi
    job_count=$(($job_count+1))
    
    logFile=/dupa-filer/laci/I1/chr$x/sim_downRate-chr$x-$targetRate.log

    echo "qsub -q all.q -V -R y -pe parallel 6 -l h_vmem=8G -l h_rt=15:00:00 -o $logFile -e $logFile -S /bin/bash /dupa-filer/laci/bin/sim_fetusRate.sh $targetRate chr$x"
    qsub -q all.q -V -R y -pe parallel 6 -l h_vmem=6G -l h_rt=15:00:00 -o $logFile -e $logFile -S /bin/bash /dupa-filer/laci/bin/sim_fetusRate.sh $targetRate chr$x
    
    sleep 30
    
done
echo "DONE."
