sed --i 's/^#PBS -J.*/#PBS -J 1-5000/' pbs.py
qsub -N I20V60R150control -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R150control/RSI20V60.mdl' pbs.py

