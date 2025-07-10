#!/bin/bash

# EDIT HERE

RECEP="5HT2B"
COND1="arr"
COND2="Gq"

# DONE EDITING

mkdir -p "plots"
mkdir -p "vispdb"
mkdir -p "results"

python ~/pensa_scripts/compare_features.py \
    --ref_file_a "traj/${RECEP}_${COND1}.gro" \
    --trj_file_a "traj/${RECEP}_${COND1}.xtc" \
    --ref_file_b "traj/${RECEP}_${COND2}.gro" \
    --trj_file_b "traj/${RECEP}_${COND2}.xtc" \
    --out_plots  "plots/${RECEP}" \
    --out_vispdb "vispdb/${RECEP}" \
    --out_results "results/${RECEP}" \
    --start_frame 0 \
    --print_num 12
