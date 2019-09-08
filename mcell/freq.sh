sed --i 's/^#PBS -J.*/#PBS -J 5000-10000/' pbs.py
qsub -N stores_blocked/RSI20V40.mdl -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/stores_blocked/RSI20V40.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 5000-10000/' pbs.py
qsub -N stores_blocked/RSI20V50.mdl -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/stores_blocked/RSI20V50.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 5000-10000/' pbs.py
qsub -N stores_blocked/RSI20V60.mdl -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/stores_blocked/RSI20V60.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 3001-10000/' pbs.py
qsub -N stores_blocked/RSI20V70.mdl -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/stores_blocked/RSI20V70.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 3001-10000/' pbs.py
qsub -N stores_blocked/RSI20V80.mdl -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/stores_blocked/RSI20V80.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 2001-10000/' pbs.py
qsub -N stores_blocked/RSI20V90.mdl -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/stores_blocked/RSI20V90.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 2001-10000/' pbs.py
qsub -N stores_blocked/RSI20V100.mdl -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/stores_blocked/RSI20V100.mdl' pbs.py

