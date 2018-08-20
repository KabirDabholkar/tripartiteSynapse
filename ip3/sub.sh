#!/bin/bash
isi=($(seq 40 10 50))
freq=($(seq 100 10 160))

for i in "${freq[@]}"
do
    qsub -N IP3SI40V${i} -v I=$i jip3.sh
done
