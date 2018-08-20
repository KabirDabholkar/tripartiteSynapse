#!/bin/bash

#PBS -p -1
#PBS -j oe
#PBS -J 1-5000:1

/apps/bin/mcell /home/nishant/serca/SI40V${I}.mdl -seed ${PBS_ARRAY_INDEX}
