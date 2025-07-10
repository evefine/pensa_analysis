import pensa
from pensa.comparison import *
from pensa.features import *
from pensa.statesinfo import *
import numpy as np
import argparse
import os
#from define_system import *
from pensa.preprocessing import load_selection, \
    extract_coordinates, extract_coordinates_combined, \
    extract_aligned_coordinates


# === Create output directories ===
#name = "CB1_arr_Gi"
#os.makedirs(name, exist_ok=True)
#os.makedirs(f"{name}/output", exist_ok=True)
os.makedirs("traj", exist_ok=True)

# EDIT HERE!!!!!1

# === SIMULATION A ===
root_dir_a = "/oak/stanford/groups/rondror/projects/MD_simulations/amber/CB1R/apo_CB1_arr/v1"
ref_file_a = f"{root_dir_a}/prep/dabble/system_dabbled.psf"
pdb_file_a = f"{root_dir_a}/prep/dabble/system_dabbled.pdb"
trj_file_a = [f"{root_dir_a}/run_{i}/summary_traj_w_eq_stride5.nc" for i in range(1, 13)]

# === SIMULATION B ===
root_dir_b = "/oak/stanford/groups/rondror/projects/MD_simulations/amber/CB1R/apo_CB1_Gi/v1"
ref_file_b = f"{root_dir_b}/prep/dabble/system_dabbled.psf"
pdb_file_b = f"{root_dir_b}/prep/dabble/system_dabbled.pdb"
trj_file_b = [f"{root_dir_b}/run_{i}/summary_traj_w_eq_stride5.nc" for i in range(1, 13)]

# === Selection strings ===
resnums = "115:140 152:176 193:251 276:306 339:363 377:400"
sel_base_a = f"(not name H*) and protein and segid P1 and resnum {resnums}"
sel_base_b = f"(not name H*) and protein and segid P6 and resnum {resnums}"

out_name_a = f"traj/CB1_arr_v1"
out_name_b = f"traj/CB1_Gi"
out_name_combined = f"traj/CB1_arr_v1_Gi"

starting_frame = 200
num_reps_a = 12
num_reps_b = 12

# done editing 


# === EXTRACT COORDINATES ===

print("Extracting CONDITION A...")
extract_coordinates(ref_file_a,
    pdb_file_a,
    trj_file_a,
    out_name_a,
    sel_base_a,
    start_frame=starting_frame
)
print('complete')

print("Extracting CONDITION B...")

extract_coordinates(ref_file_b,
    pdb_file_b,
    trj_file_b,
    out_name_b,
    sel_base_b,
    start_frame=starting_frame
)

print('complete')

print("Extracting COMBINED...")
extract_coordinates_combined(
    [ref_file_a]*num_reps_a + [ref_file_b]*num_reps_b,
    trj_file_a + trj_file_b,
    [sel_base_a]*num_reps_a + [sel_base_b]*num_reps_b,
    out_name_combined + "_combined",
    start_frame=starting_frame
)
