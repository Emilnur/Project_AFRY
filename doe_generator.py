from scipy.stats import qmc
import numpy as np
import csv
import os
import re

# === Constant Parameters ===
forces = [62000.0, 32000.0, 68000.0]
moments = [71000.0, 2000.0, 9000.0]
P = 0.3  # MPa

# === Sobol Variable Bounds ===
bounds = [
    (9, 15),             # thickness 1
    (9, 15),             # thickness 2
    (9, 15),             # thickness 3
    (2000, 2000.00001),  # VR
    (508, 508.00001),    # NR
    (2400, 2400.00001),  # H
    (300, 300.00001),    # LN
    (660, 660.00001),    # RP
]

n_levels = 10
sampler = qmc.Sobol(d=len(bounds), scramble=True)
sample = sampler.random(n_levels)
scaled = qmc.scale(sample, [b[0] for b in bounds], [b[1] for b in bounds])

# === Determine next available DOE filename ===
def get_next_doe_filename(base_dir="."):
    existing_files = os.listdir(base_dir)
    doe_indices = []

    for f in existing_files:
        if f.startswith("doe_") and f.endswith(".csv"):
            match = re.match(r"doe_(\d+)\.csv", f)
            if match:
                idx = int(match.group(1))
                doe_indices.append(idx)

    next_idx = max(doe_indices, default=0) + 1
    return f"doe_{next_idx}.csv"

output_filename = get_next_doe_filename()

# === Write to CSV ===
with open(output_filename, "w", newline="") as f:
    writer = csv.writer(f)
    
    header = [
        "Thickness_Vessel", "Thickness_Pad", "Thickness_Nozzle",
        "VR", "NR", "H", "LN", "RP",
        "Force_X", "Force_Y", "Force_Z",
        "Moment_X", "Moment_Y", "Moment_Z",
        "Pressure"
    ]
    writer.writerow(header)

    for row in scaled:
        thickness_and_geom = list(row)
        full_row = (
            thickness_and_geom +
            forces +
            moments +
            [P]
        )
        writer.writerow(full_row)

print(f"✅ DOE file generated: {output_filename}")
