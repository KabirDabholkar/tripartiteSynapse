sed --i 's/^#PBS -J.*/#PBS -J 1-2000/' pbs.py
qsub -N I20V100R150ER3x -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R150ER3x/RSI20V100.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-3000/' pbs.py
qsub -N I20V80R150ER3x -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R150ER3x/RSI20V80.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-5000/' pbs.py
qsub -N I20V60R150ER3x -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R150ER3x/RSI20V60.mdl' pbs.py

