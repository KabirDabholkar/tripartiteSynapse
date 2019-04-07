sed --i 's/^#PBS -J.*/#PBS -J 1-2000/' pbs.py
qsub -N I20V90R150control -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R150control/RSI20V90.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-2000/' pbs.py
qsub -N I40V90R150control -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R150control/RSI40V90.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-2000/' pbs.py
qsub -N I60V90R150control -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R150control/RSI60V90.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-2000/' pbs.py
qsub -N I80V90R150control -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R150control/RSI80V90.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-2000/' pbs.py
qsub -N I100V90R150control -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R150control/RSI100V90.mdl' pbs.py

