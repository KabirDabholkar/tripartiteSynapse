sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N kleak0 -v I='/home/subhadra/kabir/tripartiteSynapse/brute500_8/er_clamp_sm0_ic250_sf2_kleak0.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N kleak250 -v I='/home/subhadra/kabir/tripartiteSynapse/brute500_8/er_clamp_sm0_ic250_sf2_kleak250.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N kleak500 -v I='/home/subhadra/kabir/tripartiteSynapse/brute500_8/er_clamp_sm0_ic250_sf2_kleak500.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N kleak750 -v I='/home/subhadra/kabir/tripartiteSynapse/brute500_8/er_clamp_sm0_ic250_sf2_kleak750.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N kleak1000 -v I='/home/subhadra/kabir/tripartiteSynapse/brute500_8/er_clamp_sm0_ic250_sf2_kleak1000.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N kleak1250 -v I='/home/subhadra/kabir/tripartiteSynapse/brute500_8/er_clamp_sm0_ic250_sf2_kleak1250.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N kleak1500 -v I='/home/subhadra/kabir/tripartiteSynapse/brute500_8/er_clamp_sm0_ic250_sf2_kleak1500.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N kleak1750 -v I='/home/subhadra/kabir/tripartiteSynapse/brute500_8/er_clamp_sm0_ic250_sf2_kleak1750.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N kleak2000 -v I='/home/subhadra/kabir/tripartiteSynapse/brute500_8/er_clamp_sm0_ic250_sf2_kleak2000.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N kleak2250 -v I='/home/subhadra/kabir/tripartiteSynapse/brute500_8/er_clamp_sm0_ic250_sf2_kleak2250.mdl' pbs.py

