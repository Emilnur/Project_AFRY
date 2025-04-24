# main_pvessel.py

from abaqus import *
from abaqusConstants import *
import numpy as np
import itertools
from scipy.stats import qmc
import traceback
from prepro_functions import *
from support import *

# Parameters
n_cols = 3          # Number of columns - meaning the thickness of each component
n_levels = 10       # Number of runs
lower_bound = 9     # Minimum thickness
upper_bound = 15    # Maximum thickness

# Create Sobol sampler
sampler = qmc.Sobol(d=n_cols, scramble=True)

# Generate samples in [0, 1]^3
samples = sampler.random(n=n_levels)

# Scale to desired bounds
thickness = qmc.scale(samples, l_bounds=[lower_bound]*n_cols, u_bounds=[upper_bound]*n_cols)

forces = np.array([62000.0, 32000.0, 68000.0])
moments = np.array([71000.0, 2000.0, 9000.0])
P = 0.3  # MPa

# Run a Full Factorial Combination of Geometric Parameters
vr = np.linspace(1600, 2000, 1)
nr = np.linspace(480, 508, 1)
h = np.linspace(2300, 2400, 1)
ln = np.linspace(250, 300, 1)
rp = np.linspace(630, 660, 1)

#Full factorial design: Cartesian product of all variables
designs = list(itertools.product(vr, nr, h, ln, rp, thickness))

# === Main Execution === #
doe_dir, doe_index = create_output_structure()

for run_idx, (VR, NR, H, LN, RP, thickness) in enumerate(designs, start=1):
    run_dir = create_run_folder(doe_dir, run_idx)
    os.chdir(run_dir)
    jobname = f"pvessel_{run_idx:04d}"

    try:
        PreProcess(VR, NR, H, LN, RP, thickness, forces, moments, P, jobname)

        # Log the run
        log_run(doe_dir, doe_index, run_idx, VR, NR, H, LN, RP, thickness)

    except Exception as e:
        print(f"Error in run {run_idx}: {e}")
        traceback.print_exc()

import run_jobs_parallel
run_jobs_parallel.main()
