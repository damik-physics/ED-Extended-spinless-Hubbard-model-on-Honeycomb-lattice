# Exact Diagonalization of the Extended Hubbard Model

## Overview

This project computes the **exact many-body energy spectrum and eigenstates** of the **extended Hubbard model** on finite **2D honeycomb clusters** (regular and tilted). The numerical implementation includes:

- **Full** and **Lanczos diagonalization**
- Calculation of ground-state observables:
  - Current-current correlations
  - Density-density correlations
  - Charge and spin structure factors
- **Symmetry reduction** using the $C_{6v}$ space group:
  - Translations, rotations, mirror symmetries
  - Diagonalization within **irreducible representations**
- **OpenMP parallelization** for faster computation

---

## Project Structure

```text
.
├── Makefile               # Build script
├── input.nml              # Simulation input parameters
├── README.md              # This file
│
├── src/
│   ├── fortran/           # Fortran source files
│   │   ├── main.f90
│   │   ├── hamiltonian.f90
│   │   ├── basis.f90
│   │   └── ...
│   └── python/            # Post-processing scripts
│       └── analyze_runs.py
│
├── bin/                   # Compiled executable(s)
│   └── exe
│
├── build/                 # Compiler artifacts (ignored by Git)
│   ├── *.mod
│   ├── *.o
│   └── ...
│
├── output/                # Simulation output (not version-controlled)
│   └── run_YYYYMMDD_HHMMSS/
│       ├── input.nml
│       ├── lattice_data/
│       ├── parameters/
│       ├── logs/
│       ├── plots/
│       ├── spectra/
│       ├── states/
│       └── correlations/
│
├── doc/                   # Notes or documentation
│   └── README.md
│
└── archive/               # Deprecated/experimental code
    ├── old/
    └── dev/
