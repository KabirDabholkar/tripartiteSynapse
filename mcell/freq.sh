sed --i 's/^#PBS -J.*/#PBS -J 1-1000/' pbs.py
qsub -N R150control/ -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R150control/RS20p20hz.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-1000/' pbs.py
qsub -N R150ER2x/ -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R150ER2x/RS20p20hz.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-1000/' pbs.py
qsub -N R300ER2x/ -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R300ER2x/RS20p20hz.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-1000/' pbs.py
qsub -N R150ER3x/ -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R150ER3x/RS20p20hz.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-1000/' pbs.py
qsub -N R300ER3x/ -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R300ER3x/RS20p20hz.mdl' pbs.py

