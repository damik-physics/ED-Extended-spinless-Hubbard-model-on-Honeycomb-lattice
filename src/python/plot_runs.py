import json
import glob
import os
import matplotlib.pyplot as plt

# Lists to store extracted data
v2_values = []
energy_values = []

# Iterate through all subdirectories matching the pattern 'output/run_*'
for run_dir in glob.glob("output/run_*"):
    param_file = os.path.join(run_dir, "parameters.json")
    result_file = os.path.join(run_dir, "results.json")  # Must exist

    if not os.path.exists(param_file) or not os.path.exists(result_file):
        continue

    with open(param_file) as f:
        params = json.load(f)
    with open(result_file) as f:
        results = json.load(f)

    # Optional filter: only include runs where V == 0.0
    if params.get("V") != 0.0:
        continue

    try:
        v2 = float(params["V2"])
        energy = float(results["ground_state_energy"])
    except(KeyError, ValueError):
        continue

    v2_values.append(v2)
    energy_values.append(energy)

# Sort data by V2
sorted_indices = sorted(range(len(v2_values)), key=lambda i: v2_values[i])
v2_sorted = [v2_values[i] for i in sorted_indices]
energy_sorted = [energy_values[i] for i in sorted_indices]

# Plot
plt.figure(figsize=(6, 4))
plt.plot(v2_sorted, energy_sorted, marker='o', linestyle='-')
plt.xlabel("V2")
plt.ylabel("Ground State Energy")
plt.title("Ground State Energy vs V2(V = 0.0)")
plt.grid(True)
plt.tight_layout()
plt.savefig("ground_state_energy_vs_v2.png", dpi=150)
plt.show()
