#!/bin/bash
ISI=($(seq 1 1 25))

for i in "${ISI[@]}"
do
    qsub -N 2xNSI20d${i}V50 -v I=$i jns.sh
done
