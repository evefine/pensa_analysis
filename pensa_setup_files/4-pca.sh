#!/bin/bash

NAME=$1

mkdir -p "$NAME/pca"
mkdir -p "$NAME/plots"
mkdir -p "$NAME/results"

# bb-torsions, bb-distances, or sc-torsions

python ~/pensa_scripts/pca.py \
	--ref_file_a "$NAME/traj/CB1_arr_v1.gro" \
	--trj_file_a "$NAME/traj/CB1_arr_v1.xtc" \
	--ref_file_b "$NAME/traj/CB1_Gi.gro" \
	--trj_file_b "$NAME/traj/CB1_Gi.xtc" \
	--out_plots "$NAME/plots/CB1_v1_dists" \
	--out_pc "$NAME/pca/CB1_v1_dists" \
    --out_results "$NAME/results/CB1_v1_dists" \
	--start_frame 0 \
	--feature_type 'bb-distances' \
	--num_eigenvalues 12 \
	--num_components 3 \
	--feat_threshold 0.4 



