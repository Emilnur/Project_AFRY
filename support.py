# === Support Functions ===
import numpy as np
import os
import glob
import sys
import csv
import logging

def create_output_structure():
    cwd = os.getcwd()
    # Filter only directories that start with 'doe_'
    doe_folders = sorted(
        [d for d in os.listdir(cwd)
         if os.path.isdir(os.path.join(cwd, d)) and d.startswith('doe_')]
    )
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

def load_doe(filepath):
    designs = []
    with open(filepath, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            thickness = np.array([
                float(row["Thickness_Vessel"]),
                float(row["Thickness_Pad"]),
                float(row["Thickness_Nozzle"])
            ])
            VR = float(row["VR"])
            NR = float(row["NR"])
            H = float(row["H"])
            LN = float(row["LN"])
            RP = float(row["RP"])
            forces = np.array([
                float(row["Force_X"]),
                float(row["Force_Y"]),
                float(row["Force_Z"])
            ])
            moments = np.array([
                float(row["Moment_X"]),
                float(row["Moment_Y"]),
                float(row["Moment_Z"])
            ])
            P = float(row["Pressure"])

            designs.append((VR, NR, H, LN, RP, thickness, forces, moments, P))
    return designs

def del_abq_temp():
    # Define the patterns for the files to delete
    rec_pattern = 'abaqusi.rec*'  # Match all abaqusi.rec files
    rpy_pattern = 'abaqus.rpy.*'  # Match all abaqus.rpy.j files where j is an integer

    # Get a list of files matching the patterns
    rec_files = glob.glob(rec_pattern)
    rpy_files = glob.glob(rpy_pattern)

    # Delete the files
    for file in rec_files + rpy_files:
        try:
            os.remove(file)
            print(f"Deleted {file}")
        except Exception as e:
            print(f"Error deleting {file}: {e}")