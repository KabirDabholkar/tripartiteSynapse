#!/bin/bash
ISI=($(seq 120 20 200))
freq=($(seq 50 10 70))

for i in "${freq[@]}"
do
    qsub -N SI40V${i} -v I=$i jserca.sh
done
