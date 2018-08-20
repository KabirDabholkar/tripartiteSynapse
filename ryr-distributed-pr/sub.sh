#!/bin/bash

i=250
while read -a line
do
    n="${line[1]}"
    vd="${line[3]}"
    n=$((n*3))
    #echo $vd
    
    #m=1000
    sed -i "s/-J.*/-J 1-$n/; s/V.*dist/V${vd}dist/" jryr.sh
	#m=$n
	#echo $m
	#sleep 1 
	
	qsub -N RSI${i}dist -v I=$i jryr.sh

done < pr-dist.dat

    

