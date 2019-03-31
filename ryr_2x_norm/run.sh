sed --i 's/^#PBS -J.*/#PBS -J 1-500/' pbs.py
qsub -N RSI20V90norm -v I='/home/subhadra/kabir/tripartiteSynapse/ryr_2x_norm/RSI20V90.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-500/' pbs.py
qsub -N RSI30V90norm -v I='/home/subhadra/kabir/tripartiteSynapse/ryr_2x_norm/RSI30V90.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-500/' pbs.py
qsub -N RSI40V90norm -v I='/home/subhadra/kabir/tripartiteSynapse/ryr_2x_norm/RSI40V90.mdl' pbs.py

