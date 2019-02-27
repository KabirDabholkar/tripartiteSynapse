sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N sm50_ic200_fc500 -v I='/home/subhadra/kabir/tripartiteSynapse/brute750_9/RSnostim_sm50_ic200_fc500.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N sm100_ic200_fc500 -v I='/home/subhadra/kabir/tripartiteSynapse/brute750_9/RSnostim_sm100_ic200_fc500.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N sm125_ic200_fc500 -v I='/home/subhadra/kabir/tripartiteSynapse/brute750_9/RSnostim_sm125_ic200_fc500.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N sm150_ic200_fc500 -v I='/home/subhadra/kabir/tripartiteSynapse/brute750_9/RSnostim_sm150_ic200_fc500.mdl' pbs.py

