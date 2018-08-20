#!/bin/bash
isi=($(seq 20 10 100))
freq=($(seq 70 10 80))


for i in "${freq[@]}"
do
    qsub -N IRSI40V${i} -v I=$i jirs.sh
done
