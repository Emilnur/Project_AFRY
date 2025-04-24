# main_pvessel.py

from abaqus import *
from abaqusConstants import *
import numpy as np
import itertools
from scipy.stats import qmc
import glob
import sys
import os
import csv

# Add current directory to sys.path so imports work inside Abaqus CAE GUI - not working
script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)

from pvessel_functions import *

# === Setup folder structure ===

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

# === Main Execution ===

doe_dir, doe_index = create_output_structure()

# Parameters
n_cols = 3          # Number of columns - meaning the thickness of each component
n_levels = 10       # Number of runs
lower_bound = 9     # Minimum thickness
upper_bound = 15    # Maximum thickness

# Create Sobol sampler
sampler = qmc.Sobol(d=n_cols, scramble=True)

# Generate samples in [0, 1]^3
samples = sampler.random(n=n_levels)

# Scale to desired bounds [9, 15]
thickness = qmc.scale(samples, l_bounds=[lower_bound]*n_cols, u_bounds=[upper_bound]*n_cols)

forces = np.array([62000.0, 32000.0, 68000.0])
moments = np.array([71000.0, 2000.0, 9000.0])
P = 0.3  # MPa

# Run a Full Factorial Combination of Geometric Parameters
vr = np.linspace(1600, 2000, 2)
nr = np.linspace(480, 508, 2)
h = np.linspace(2300, 2400, 2)
ln = np.linspace(250, 300, 2)
rp = np.linspace(630, 660, 2)

#Full factorial design: Cartesian product of all variables
designs = list(itertools.product(vr, nr, h, ln, rp,thickness))

for run_idx, (VR, NR, H, LN, RP, thickness) in enumerate(designs, start=1):
    run_dir = create_run_folder(doe_dir, run_idx)
    os.chdir(run_dir)
    jobname = f"pvessel_{run_idx:04d}"

    try:
        # Create new model
        Mdb()
        session.viewports['Viewport: 1'].setValues(displayedObject=None)
        mymodel = mdb.models[mdb.models.keys()[0]]

        # Validate and create geometry
        VR, NR, H, LN, RP, thickness, xsi_N = validate_geom(VR, NR, H, LN, RP, thickness)
        flatA = Parts(mymodel, VR, NR, H, LN, RP, thickness)
        Assembly(mymodel, thickness)
        Property(mymodel, thickness)
        Step(mymodel)
        Interactions(mymodel)
        Loads(mymodel, forces, moments, P, flatA, VR)
        Mesh(mymodel)
        InputFile(mymodel, job_name=jobname)
        
        # Log the run
        log_run(doe_dir, doe_index, run_idx, VR, NR, H, LN, RP, thickness)

    except:
       pass