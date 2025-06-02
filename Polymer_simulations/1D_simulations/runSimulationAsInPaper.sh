#!/bin/bash




#this script runs in 1 min if showPlot="False" or 5 min if showPlot="False"

number_simul=2




#current folder:
pathdata="./example"

sizeKb_list="345 566 918 576 "
processivityKb_list="330"
cohesinNumberPerMb_list="12"
avg_CTCF_attach_time_list="150"
extrusion_speedKbS="0.25"
proba_CTCF_occupancy_list="0.5"
simulation_path=$pathdata'/simulation'
showPlot="False"  # "True" or "False"

now="$(date +"%T")"
echo "Start time : $now"


for i in `seq 0 $number_simul`
do
  echo "run sim "$i

  for proc in $processivityKb_list; do

      for coh in $cohesinNumberPerMb_list; do
      	for ctcf in $avg_CTCF_attach_time_list; do
              for sizeKb in $sizeKb_list; do
                  for extr in $extrusion_speedKbS; do
                      for proba_CTCF_occupancy in  $proba_CTCF_occupancy_list; do

                                python  runSimulationBeforeLammps.py $proc $coh $ctcf $sizeKb $extr  $proba_CTCF_occupancy $i $simulation_path $showPlot

                      done
                  done
              done
          done
      done
  done
done

now="$(date +"%T")"
echo "End time : $now"
