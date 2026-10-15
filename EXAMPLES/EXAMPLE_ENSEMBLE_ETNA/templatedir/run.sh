#!/bin/bash
echo $PWD

python create_input_ellipsoid.py > preprocessing.log

export OMP_NUM_THREADS=1
./IMEX_SfloW2D > run.log

# Optional per-run plots; rasters are stored in OUTPUT/<RUN_NAME>/.
#python plot_overlay.py > postprocessing1.log
