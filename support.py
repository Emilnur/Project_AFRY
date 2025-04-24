# === Support Functions ===

import os
import glob
import sys
import csv

def create_output_structure():
    cwd = os.getcwd()
    doe_folders = sorted(glob.glob(os.path.join(cwd, 'doe_*')))
    next_doe_index = len(doe_folders) + 1
    doe_name = f"doe_{next_doe_index}"
    doe_path = os.path.join(cwd, doe_name)
    os.makedirs(doe_path)
    return doe_path, next_doe_index

def create_run_folder(base_path, run_index):
    run_name = f"run_{run_index:04d}"
    run_path = os.path.join(base_path, run_name)
    os.makedirs(run_path)
    return run_path

def log_run(doe_dir, doe_index, run_index, VR, NR, H, LN, RP, thickness):
    log_path = os.path.join(doe_dir, "doe_log.csv")
    header = ["DOE_Index", "Run_Index", "VR", "NR", "H", "LN", "RP", "Thickness_Vessel", "Thickness_Pad", "Thickness_Nozzle"]
    row = [doe_index, run_index, VR, NR, H, LN, RP, thickness[0], thickness[1], thickness[2]]

    write_header = not os.path.exists(log_path)

    with open(log_path, mode='a', newline='') as file:
        writer = csv.writer(file)
        if write_header:
            writer.writerow(header)
        writer.writerow(row)
