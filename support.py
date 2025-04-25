# === Support Functions ===
import numpy as np
import os
import re
import glob
import sys
import csv
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import shutil

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
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Patterns to match
    patterns = [
        "abaqus.rpy",         # plain file
        "abaqus.rpy.*",       # numbered files
        "abaqus*.rec"         # rec files like abaqus1.rec, abaqus2.rec, etc.
    ]

    # Delete matching files
    for pattern in patterns:
        for file_path in glob.glob(os.path.join(script_dir, pattern)):
            try:
                os.remove(file_path)
                print(f"Deleted: {file_path}")
            except Exception as e:
                print(f"Failed to delete {file_path}: {e}")

    # Delete __pycache__ folder
    pycache_dir = os.path.join(script_dir, "__pycache__")
    if os.path.exists(pycache_dir) and os.path.isdir(pycache_dir):
        try:
            shutil.rmtree(pycache_dir)
            print(f"Deleted folder: {pycache_dir}")
        except Exception as e:
            print(f"Failed to delete folder {pycache_dir}: {e}")

def run_abaqus_job(folder_path, inp_file, job_name):
    try:
        print(f"Running {job_name} in {os.path.basename(folder_path)}...")
        run_job_script = os.path.join(os.path.dirname(__file__), "run_job.py")
        command = f'abaqus cae noGUI="{run_job_script}" -- "{inp_file}" "{job_name}"'
        subprocess.run(command, cwd=folder_path, shell=True, check=True)
        print(f"Done: {job_name}")
    except subprocess.CalledProcessError as e:
        print(f"Failed: {job_name} — {e}")

def parallel_exec(base_path, max_parallel_jobs):

    # Find all folders with .inp files and no .odb file
    jobs_to_run = []

    for folder_name in os.listdir(base_path):
        folder_path = os.path.join(base_path, folder_name)
        if not os.path.isdir(folder_path):
            continue

        inp_files = [f for f in os.listdir(folder_path) if f.endswith('.inp')]
        if not inp_files:
            continue

        inp_file = inp_files[0]
        job_name = os.path.splitext(inp_file)[0]
        odb_path = os.path.join(folder_path, job_name + ".odb")

        if not os.path.exists(odb_path):
            jobs_to_run.append((folder_path, inp_file, job_name))
        else:
            print(f"{job_name}.odb already exists — skipping.")

    with ThreadPoolExecutor(max_workers=max_parallel_jobs) as executor:
        futures = [executor.submit(run_abaqus_job, folder, inp, job)
                   for folder, inp, job in jobs_to_run]
        for future in as_completed(futures):
            _ = future.result()

    print("All jobs finished.")

    # Optional: clean up extra files
    extensions_to_keep = {'.inp', '.odb'}
    print("Cleaning up folders...")
    for folder_name in os.listdir(base_path):
        folder_path = os.path.join(base_path, folder_name)
        if not os.path.isdir(folder_path):
            continue
        for file_name in os.listdir(folder_path):
            file_path = os.path.join(folder_path, file_name)
            if os.path.isfile(file_path):
                ext = os.path.splitext(file_name)[1].lower()
                if ext not in extensions_to_keep:
                    try:
                        os.remove(file_path)
                        print(f"Deleted: {file_path}")
                    except Exception as e:
                        print(f"Could not delete {file_path}: {e}")

    print("Cleanup complete.")

def select_next_doe_csv():
    current_dir = os.getcwd()

    # Find all doe_*.csv files in the current directory
    csv_files = sorted(
        [f for f in os.listdir(current_dir) if re.match(r"doe_\d+\.csv$", f)],
        key=lambda x: int(re.search(r"\d+", x).group())
    )

    if not csv_files:
        print("No design of experiments file (doe_*.csv) found in the current directory.")
        return None

    # Find all folders named doe_*
    doe_folders = {
        f for f in os.listdir(current_dir)
        if os.path.isdir(f) and re.match(r"doe_\d+$", f)
    }

    # Look for the first csv file without a matching folder
    for csv_file in csv_files:
        folder_name = csv_file.rsplit('.', 1)[0]
        if folder_name not in doe_folders:
            print(f"Selected DOE file: {csv_file}")
            return csv_file

    print("No design of experiments left to run — all corresponding folders exist.")
    return None
