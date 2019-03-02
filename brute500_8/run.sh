sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N sm30_ic250 -v I='/home/subhadra/kabir/tripartiteSynapse/brute500_8/er_clamp_sm30_ic250_sf2_kleak0.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N sm60_ic250 -v I='/home/subhadra/kabir/tripartiteSynapse/brute500_8/er_clamp_sm60_ic250_sf2_kleak0.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N sm90_ic250 -v I='/home/subhadra/kabir/tripartiteSynapse/brute500_8/er_clamp_sm90_ic250_sf2_kleak0.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N sm120_ic250 -v I='/home/subhadra/kabir/tripartiteSynapse/brute500_8/er_clamp_sm120_ic250_sf2_kleak0.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N sm150_ic250 -v I='/home/subhadra/kabir/tripartiteSynapse/brute500_8/er_clamp_sm150_ic250_sf2_kleak0.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N sm180_ic250 -v I='/home/subhadra/kabir/tripartiteSynapse/brute500_8/er_clamp_sm180_ic250_sf2_kleak0.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N sm210_ic250 -v I='/home/subhadra/kabir/tripartiteSynapse/brute500_8/er_clamp_sm210_ic250_sf2_kleak0.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N sm240_ic250 -v I='/home/subhadra/kabir/tripartiteSynapse/brute500_8/er_clamp_sm240_ic250_sf2_kleak0.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N sm270_ic250 -v I='/home/subhadra/kabir/tripartiteSynapse/brute500_8/er_clamp_sm270_ic250_sf2_kleak0.mdl' pbs.py

