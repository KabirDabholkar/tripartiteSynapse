sed --i 's/^#PBS -J.*/#PBS -J 1-1000/' pbs.py
qsub -N kleak0 -v I='/home/subhadra/kabir/tripartiteSynapse/AD/RSI20V40.mdl' pbs.py

