#!/bin/bash
#PBS -N yourjobname
#PBS -p -1
#PBS -j oe
#PBS -t 1-5

/apps/bin/mcell -input=file-${PBS_ARRAYID}
