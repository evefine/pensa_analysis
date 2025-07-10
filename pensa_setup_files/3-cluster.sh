#!/bin/bash

NAME=$1
mkdir -p "$NAME/plots"
mkdir -p "$NAME/clusters"
# change to --write when want to write it, but for wss analysis 
# dont need to 
python ~/pensa_scripts/cluster.py --write --wss \
	--ref_file_a "$NAME/traj/CB1_arr_v1.gro" \
	--trj_file_a "$NAME/traj/CB1_arr_v1.xtc" \
	--ref_file_b "$NAME/traj/CB1_Gi.gro" \
	--trj_file_b "$NAME/traj/CB1_Gi.xtc" \
	--label_a 'CB1-arr_v1' \
	--label_b 'CB1-Gi' \
	--out_plots "$NAME/plots/CB1_v1" \
    --out_results "$NAME/results/CB1_v1" \
	--out_frames_a "$NAME/clusters/CB1_arr_v1" \
	--out_frames_b "$NAME/clusters/CB1_Gi" \
	--start_frame 0 \
	--feature_type 'bb-distances' \
	--algorithm 'kmeans' \
	--max_num_clusters 12 \
	--write_num_clusters 3
