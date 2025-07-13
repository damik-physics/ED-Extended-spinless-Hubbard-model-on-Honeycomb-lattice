import json
import numpy as np
import glob
import os

for run_dir in glob.glob("output/run_*"):
    param_file = os.path.join(run_dir, "parameters.json")
    if not os.path.exists(param_file):
        continue
    with open(param_file) as f:
        params = json.load(f)
    if params["V"] == 0.0 and params["V2"] == 0.3:
        energy_file = os.path.join(run_dir, "energy.dat")
        if os.path.exists(energy_file):
            data = np.loadtxt(energy_file)
            print(f"Run {run_dir} energy shape:", data.shape)
