#!/bin/bash

for i in {1..2}
do
   mkdir set$i
done

for i in {1..2}
do
   cp simIndex=$((i-1)).in set$i
   cp interactions set$i
   cp initial.py set$i
   cp input_data_generate_f1.py set$i
   cd set$i/
   python3 initial.py > initial_conformation.txt
   python3 input_data_generate_f1.py > param
   sbatch -c 1 --qos=specific --job-name=L2P1100C8$i  -t 5-00:00:00 --mem 2G --wrap='lmp_serial -in simIndex='$((i-1))'.in &> chain.txt'
   cd ../
done
