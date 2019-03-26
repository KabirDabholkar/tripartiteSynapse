sed --i 's/^#PBS -J.*/#PBS -J 1-50/' pbs.py
qsub -N pred_final750_eq -v I='/home/subhadra/kabir/tripartiteSynapse/final/pred_final750_eq.mdl' pbs.py

