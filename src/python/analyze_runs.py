import json
import numpy as np
import glob
import os

# Iterate through all subdirectories matching the pattern 'output/run_*'
for run_dir in glob.glob("output/run_*"):
    param_file = os.path.join(run_dir, "parameters.json")

    # Skip if the parameters file is missing
    if not os.path.exists(param_file):
        continue

    # Open and parse the JSON file containing run parameters
    with open(param_file) as f:
        params = json.load(f)

    # Only print the parameters for specific runs, e.g., V = 0.0 and V2 = 0.3
    if params["V"] == 0.0 and params["V2"] == 0.3:
        print(f"Found match in {run_dir}:")
        print(json.dumps(params, indent=2))
