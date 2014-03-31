#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o pipefail

listFile=$1
chr=$2
cnvType=$3

data_path="/dupa-filer/laci/I1/$chr/"
plasma_path="/dupa-filer/laci/I1/$chr/sim_plasma07/"
results_path="/dupa-filer/laci/I1/$chr/ff10/"
exec_path="/dupa-filer/laci/bin/"
if [[ ! -d "$data_path" || ! -d "$plasma_path" ]]; then
    echo "ERROR: data_path or plasma_path doesn't exist!"
    exit
fi
export SGE_O_PATH=$PATH

plasmaFile=__plasma-chr1-0.07rate.sort.bam
ratio=0.07

#================================== DELETIONS ==================================
if [[ "$cnvType" == DEL ]] ; then

#A-chr1-36415257-46415257-delete.target.txt
#fnms=(`cat $listFile | tr '\n' ' '`)
begs=(`awk -v FS="-" '{print $3}' $listFile | tr '\n' ' '`)
lens=(`awk -v FS="-" '{print $4-$3}' $listFile | tr '\n' ' '`)
haps=(`awk -v FS="-" '{print $1}' $listFile | tr '\n' ' '`)

echo "Starting CNV plasma files simulations: $cnvType"
for (( i=0; i < ${#begs[@]}; i++ )) 
do
    begin=${begs[$i]}
    length=${lens[$i]}
    haplo=${haps[$i]}
    #src=${srcs[$i]}
    
    echo $begin' '$length' '$haplo
    echo "qsub -q all.q -V -R y -l h_vmem=6G -l h_rt=15:00:00 -S /bin/bash $exec_path/sim_deletion.sh $chr $begin $length $haplo $data_path $plasma_path $exec_path $plasmaFile $ratio"
    qsub -q all.q -V -R y -l h_vmem=6G -l h_rt=15:00:00 -S /bin/bash $exec_path/sim_deletion.sh $chr $begin $length $haplo $data_path $plasma_path $exec_path $plasmaFile $ratio
    
    sleep 5
done

echo "DONE."
exit
fi

#================================= DUPLICATIONS ================================
if [[ "$cnvType" == DUP ]] ; then
#IP1-A-chr1-148349432-148449432-duplicate

#fnms=(`cat $listFile | tr '\n' ' '`)
begs=(`awk -v FS="-" '{print $4}' $listFile | tr '\n' ' '`)
lens=(`awk -v FS="-" '{print $5-$4}' $listFile | tr '\n' ' '`)
haps=(`awk -v FS="-" '{print $2}' $listFile | tr '\n' ' '`)
srcs=(`awk -v FS="-" '{print substr($1,2,1)}' $listFile | tr '\n' ' '`)

echo "Starting CNV plasma files simulations: $cnvType"
for (( i=0; i < ${#begs[@]}; i++ )) 
do
    begin=${begs[$i]}
    length=${lens[$i]}
    haplo=${haps[$i]}
    src=${srcs[$i]}
    
    echo $begin' '$length' '$haplo
    echo "qsub -q all.q -V -R y -l h_vmem=6G -l h_rt=15:00:00 -S /bin/bash $exec_path/sim_duplicate.sh $chr $begin $length $src $haplo $data_path $plasma_path $exec_path $plasmaFile $ratio"
    qsub -q all.q -V -R y -l h_vmem=6G -l h_rt=15:00:00 -S /bin/bash $exec_path/sim_duplicate.sh $chr $begin $length $src $haplo $data_path $plasma_path $exec_path $plasmaFile $ratio
    
    sleep 5
done

echo "DONE."
exit
fi

#$exec_path/sim_duplicate.sh 10000000 100000 P B $data_path $plasma_path $exec_path &
#$exec_path/sim_deletion.sh 10000000 100000 A $data_path $plasma_path $exec_path &
#wait


#beg=(98974829 229790173 202780071 206112235 187315211 179758438 150984880 201220013 148349432 212168453 161331296 234162187)
#len=(5000000 500000 100000 100000 100000 500000 1000000 500000 100000 1000000 100000 500000)
#srccc=(P M P M M P P P P M P M)
#hap=(B B B A B A B B A B B A)

#for i in {0..11}; do
#    begin=${beg[$i]}
#    length=${len[$i]}
#    src=${srccc[$i]}
#    haplo=${hap[$i]}
#    echo "qsub -q all.q -V -R y -l h_vmem=6G -l h_rt=15:00:00 -S /bin/bash $exec_path/sim_duplicate.sh $chr $begin $length $src $haplo $data_path $plasma_path $exec_path"
#done



