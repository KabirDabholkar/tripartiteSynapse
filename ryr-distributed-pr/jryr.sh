#!/bin/bash

#PBS -p 1
#PBS -j oe
##PBS -o /storage/nishant/garbage/
#PBS -J 1-3

/apps/bin/mcell /home/nishant/ryr-distributed-pr/RSI${I}V176dist.mdl -seed ${PBS_ARRAY_INDEX}
