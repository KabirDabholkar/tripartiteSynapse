sed --i 's/^#PBS -J.*/#PBS -J 1-2000/' pbs.py
qsub -N I20V90R300ER2x -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R300ER2x/RSI20V90.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-2000/' pbs.py
qsub -N I40V90R300ER2x -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R300ER2x/RSI40V90.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-2000/' pbs.py
qsub -N I60V90R300ER2x -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R300ER2x/RSI60V90.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-2000/' pbs.py
qsub -N I80V90R300ER2x -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R300ER2x/RSI80V90.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-2000/' pbs.py
qsub -N I100V90R300ER2x -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R300ER2x/RSI100V90.mdl' pbs.py

