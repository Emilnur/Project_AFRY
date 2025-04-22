import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

# Setup
base_path = r"C:\Users\emilsjos\Documents\Code\Python\Project_Batch_mode\InputFiles"
max_parallel_jobs = 4  # Adjust this based on how many cores you want to use

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

# Define how to run a single job using Abaqus CAE
def run_abaqus_job(folder_path, inp_file, job_name):
    try:
        print(f"Running {job_name} in {os.path.basename(folder_path)}...")
        run_job_script = os.path.abspath("run_job.py")
        command = f'abaqus cae noGUI="{run_job_script}" -- "{inp_file}" "{job_name}"'
        subprocess.run(command, cwd=folder_path, shell=True, check=True)
        print(f"Done: {job_name}")
    except subprocess.CalledProcessError as e:
        print(f"Failed: {job_name} — {e}")

# Use threads to run jobs in parallel
with ThreadPoolExecutor(max_workers=max_parallel_jobs) as executor:
    futures = [executor.submit(run_abaqus_job, folder, inp, job)
               for folder, inp, job in jobs_to_run]
    for future in as_completed(futures):
        _ = future.result()

print("All jobs finished.")
