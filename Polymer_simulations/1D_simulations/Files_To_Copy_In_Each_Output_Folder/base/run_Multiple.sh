#!/bin/bash

for i in {1..50}
do
   mkdir set$i
done

for i in {1..50}
do
   cp simIndex=$((i-1)).in set$i
   cp interactions set$i
   cp initial.py set$i
   cp input_data_generate_f1.py set$i
   cd set$i/
   python3 initial.py > initial_conformation.txt
   python3 input_data_generate_f1.py > param
   #nohup lmp_serial -in lammps_input_file.in &> sim$i.txt &
   #nohup python Pf_lammps.py &> hu_$i.txt &
   sbatch -c 1 --qos=specific --job-name=PX_$i  -t 10-00:00:00 --mem 2G --wrap='lmp_serial -in simIndex='$((i-1))'.in &> chain.txt'
   cd ../
done
