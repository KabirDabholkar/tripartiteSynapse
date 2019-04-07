#!/bin/bash

ids=($(seq $1 $2))

for i in "${ids[@]}"
do
    qdel ${i}[]
done

