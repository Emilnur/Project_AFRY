# main_pvessel.py

import traceback
from support import *

csvfile = select_next_doe_csv()
designs = load_doe(csvfile)

## ====== Main Execution ====== ##

# === Pre-Processing === #
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

# === Running Job === #
parallel_exec(doe_dir, 8)

# === Post-Processing === #
odb_list = find_odb(doe_dir)

for odb_f, odb_dir in odb_list:
    try:
        PostProcess(odb_f, odb_dir, 1)
    except:
        continue

# === Cleanup === #
del_abq_temp()