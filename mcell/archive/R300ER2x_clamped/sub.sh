#!/bin/bash
ISI=($(seq 120 20 200))
freq=($(seq 2 2 24))

for i in "${freq[@]}"
do
    qsub -N 2xRSI20d${i} -v I=$i jryr.sh
done
