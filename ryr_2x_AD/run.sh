sed --i 's/^#PBS -J.*/#PBS -J 1-500/' pbs.py
qsub -N RSI20V90AD -v I='/home/subhadra/kabir/tripartiteSynapse/ryr_2x_AD/RSI20V90.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-500/' pbs.py
qsub -N RSI30V90AD -v I='/home/subhadra/kabir/tripartiteSynapse/ryr_2x_AD/RSI30V90.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-500/' pbs.py
qsub -N RSI40V90AD -v I='/home/subhadra/kabir/tripartiteSynapse/ryr_2x_AD/RSI40V90.mdl' pbs.py

