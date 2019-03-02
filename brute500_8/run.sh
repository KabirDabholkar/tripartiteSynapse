sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N er_clamp_sm0_ic260 -v I='/home/subhadra/kabir/tripartiteSynapse/brute500_8/er_clamp_sm0_ic260_sf2_kleak0.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N er_clamp_sm0_ic263 -v I='/home/subhadra/kabir/tripartiteSynapse/brute500_8/er_clamp_sm0_ic263_sf2_kleak0.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N er_clamp_sm0_ic266 -v I='/home/subhadra/kabir/tripartiteSynapse/brute500_8/er_clamp_sm0_ic266_sf2_kleak0.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N er_clamp_sm0_ic269 -v I='/home/subhadra/kabir/tripartiteSynapse/brute500_8/er_clamp_sm0_ic269_sf2_kleak0.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-20/' pbs.py
qsub -N er_clamp_sm0_ic272 -v I='/home/subhadra/kabir/tripartiteSynapse/brute500_8/er_clamp_sm0_ic272_sf2_kleak0.mdl' pbs.py

