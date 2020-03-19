sed --i 's/^#PBS -J.*/#PBS -J 1-500/' pbs.py
qsub -N R150control_20V180 -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R150control_clamped/RSI20V180.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-500/' pbs.py
qsub -N R150control_20V200 -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R150control_clamped/RSI20V200.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-500/' pbs.py
qsub -N R150ER2x_20V180 -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R150ER2x_clamped/RSI20V180.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-500/' pbs.py
qsub -N R150ER2x_20V200 -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R150ER2x_clamped/RSI20V200.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-500/' pbs.py
qsub -N R300ER2x_20V180 -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R300ER2x_clamped/RSI20V180.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-500/' pbs.py
qsub -N R300ER2x_20V200 -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R300ER2x_clamped/RSI20V200.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-500/' pbs.py
qsub -N R150ER3x_20V180 -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R150ER3x_clamped/RSI20V180.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-500/' pbs.py
qsub -N R150ER3x_20V200 -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R150ER3x_clamped/RSI20V200.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-500/' pbs.py
qsub -N R300ER3x_20V180 -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R300ER3x_clamped/RSI20V180.mdl' pbs.py

sed --i 's/^#PBS -J.*/#PBS -J 1-500/' pbs.py
qsub -N R300ER3x_20V200 -v I='/home/subhadra/kabir/tripartiteSynapse/mcell/R300ER3x_clamped/RSI20V200.mdl' pbs.py

