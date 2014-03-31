#!/bin/bash

data_path='/dupa-filer/laci/I1/chr20/'
plasma_path='/dupa-filer/laci/I1/chr20/test/'
results_path='/dupa-filer/laci/I1/chr20/test/'
exec_path='/dupa-filer/laci/bin/'

time $exec_path/sim_duplicate.sh 35000000 100000 M A $data_path $plasma_path $exec_path

#time ./sim_duplicate.sh 10000000 100000 P A $data_path $plasma_path $exec_path &
#time ./sim_duplicate.sh 10000000 1000000 P B $data_path $plasma_path $exec_path &
#time ./sim_duplicate.sh 10000000 10000000 P A $data_path $plasma_path $exec_path &

#time ./sim_duplicate.sh 10000000 100000 M A $data_path $plasma_path $exec_path &
#time ./sim_duplicate.sh 10000000 1000000 M B $data_path $plasma_path $exec_path &
#time ./sim_duplicate.sh 10000000 10000000 M A $data_path $plasma_path $exec_path &

