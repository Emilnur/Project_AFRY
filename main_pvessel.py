# main_pvessel.py

import traceback
from prepro_functions import *
from support import *

designs = load_doe("doe_1.csv")

# === Main Execution === #
doe_dir, doe_index = create_output_structure()

for run_idx, (VR, NR, H, LN, RP, thickness, forces, moments, P) in enumerate(designs, start=1):
    run_dir = create_run_folder(doe_dir, run_idx)
    os.chdir(run_dir)
    jobname = f"pvessel_{run_idx:04d}"

    try:
        PreProcess(VR, NR, H, LN, RP, thickness, forces, moments, P, jobname)

    except Exception as e:
        print(f"Error in run {run_idx}: {e}")
        traceback.print_exc()


#del_abq_temp()