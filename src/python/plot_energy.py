#!/usr/bin/env python3
"""
plot_energy.py

Plots energy levels as a function of V2 for a given V1 and simulation parameters.

Usage:
    python plot_energy.py --output-dir OUTPUT_DIR \
                          --filter ucx=3 ucy=3 cluster=18C irrep=A1 filling=0.5 \
                          --v1 1.0 --v2min 0.0 --v2max 1.0 \
                          [--out FILE.png]

Arguments:
    --output-dir      Top-level folder containing run_* directories.
    --filter          Key=value pairs to identify the run (e.g., ucx=3 ucy=3 cluster=18C ...).
    --v1              The V1 value to plot.
    --v2min           Minimum V2 value to include.
    --v2max           Maximum V2 value to include.
    --out             Optional filename to save the figure instead of using run folder.

Example:
    python plot_energy.py --output-dir output \
                          --filter ucx=3 ucy=3 cluster=18C irrep=A1 filling=0.5 \
                          --v1 1.0 --v2min 0.0 --v2max 1.0 \
                          --out energy_plot.png
"""

import os
import json
import pandas as pd
import matplotlib.pyplot as plt
import argparse

def parse_filter(filter_list):
    filt = {}
    for item in filter_list:
        if '=' not in item:
            raise ValueError(f"Filter item {item} must be key=value")
        k, v = item.split('=', 1)
        try:
            v = float(v)
        except ValueError:
            pass
        filt[k] = v
    return filt

def find_run_folder(output_dir, filter_dict):
    for run in os.listdir(output_dir):
        run_path = os.path.join(output_dir, run)
        if not os.path.isdir(run_path) or not run.startswith('run_'):
            continue
        param_file = os.path.join(run_path, 'parameters', 'parameters.json')
        if not os.path.isfile(param_file):
            continue
        with open(param_file) as f:
            params = json.load(f)
        match = True
        for k, v in filter_dict.items():
            if k not in params or params[k] != v:
                match = False
                break
        if match:
            return run_path
    return None

def main():
    parser = argparse.ArgumentParser(description="Plot energy levels vs V2")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--filter", nargs='+', required=True, help="Key=value filters to select run")
    parser.add_argument("--v1", type=float, required=True)
    parser.add_argument("--v2min", type=float, default=None)
    parser.add_argument("--v2max", type=float, default=None)
    parser.add_argument("--out", type=str, default=None)
    args = parser.parse_args()

    filter_dict = parse_filter(args.filter)
    run_folder = find_run_folder(args.output_dir, filter_dict)
    if not run_folder:
        raise FileNotFoundError("No matching run folder found with the given parameters")

    csv_file = os.path.join(run_folder, 'spectra', 'energy.csv')
    if not os.path.isfile(csv_file):
        raise FileNotFoundError(f"energy.csv not found in {csv_file}")

    df = pd.read_csv(csv_file)
    
    # Determine V2 range
    v2 = df['V2']
    mask = (v2 >= (args.v2min if args.v2min is not None else v2.min())) & \
           (v2 <= (args.v2max if args.v2max is not None else v2.max()))
    df = df[mask]

    # Plot all energy columns (skip V2 column)
    energy_cols = [c for c in df.columns if c != 'V2']
    for col in energy_cols:
        plt.plot(df['V2'], df[col], label=col)
    
    plt.xlabel("V2")
    plt.ylabel("Energy")
    plt.title(f"Energy levels vs V2 for V1={args.v1}")
    plt.legend()

    # Default output: run_folder/plots/energy_plot.png
    if args.out:
        out_file = args.out
    else:
        plots_dir = os.path.join(run_folder, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        out_file = os.path.join(plots_dir, 'energy_plot.png')

    plt.savefig(out_file)
    print(f"Plot saved to {out_file}")
    plt.show()

if __name__ == "__main__":
    main()
