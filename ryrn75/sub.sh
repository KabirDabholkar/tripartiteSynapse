#!/bin/bash
ISI=($(seq 120 20 200))
freq=($(seq 110 10 160))

for i in "${freq[@]}"
do
	 qsub -N R75SI20V${i} -v I=$i jryr.sh
done
