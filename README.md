# **Exact Diagonalization of the Extended Hubbard Model**  

## **Overview**  
This project computes the **exact many-body spectrum and eigenstates** of the **extended Hubbard model** in **2D** on both **regular** and **tilted honeycomb clusters**. The implementation includes:  

- **Lanczos** and **full diagonalization** methods  
- Calculation of various **ground-state observables**, such as:  
  - **Current-current correlation function**  
  - **Density-density correlation function (under construction)**  
  - **Charge and spin correlations (under construction)**  
  <!-- Extension to include entanglement entropy, level spacing and inverse partition ratio is under work. -->
- **Symmetrization under the** $C_{6v}$ **space group**(translations, rotations, and mirror symmetries)
- **Hamiltonian diagonalization within individual irreducible representations**  
- **Parallelization with OpenMP** for improved performance  

---

## 📁 Project Directory

```text
.
├── Makefile               # Main build file
├── input.nml              # List of input parameters
├── README.md              # Project documentation
│
├── src/                   # Source code
│   ├── fortran/           # Fortran modules
│   │   ├── basis.f90
│   │   ├── hamiltonian.f90
│   │   ├── ...
│   └── python/            # Plotting scripts
│       └── plot_energy.py
│
├── bin/                   # Compiled binaries
│   └── exe
│
├── build/                 # Build artifacts
│   ├── *.o                # Removed with 'make clean'
│   ├── *.mod
│   └── ...
│
├── output/                # Simulation results 
│   └── run_YYYYMMDD_HHMMSS/
│       ├── input.nml      # User input parameters
│       ├── parameters/parameters.json # Relevant run parameters 
│       ├── lattice_data/
│       ├── logs/
│       ├── plots/
│       ├── hamiltonians/
│       ├── correlations/
│       ├── spectra/
│       └── states/
│
└── README.md               # Project documentation
   

```
### **Modules and Functions**  

- **` main.f90`**: Program entry point
- **` basis.f90`**: Basis states and related operations
- **` diagonalization.f90`**: Diagonalization routines
- **` hamiltonian.f90`**: Hamiltonian definitions and matrix elements
- **` lattice.f90`**: Lattice geometry and connectivity
- **` symmetries.f90`**: Symmetry operations and irreps
- **` observables.f90`**: Correlation functions (currently only current-current correlations)
- **` types.f90`**: Derived types and data structures
- **` functions.f90`**: Helper functions specific to the physics
- **` params.f90`**: Simulation constants and variables not defined as derived types.
- **` core_utilities.f90`**: Essential utility routines used throughout simulation.
- **` io_utilities.f90`**: Subroutines used for I/O handling and data output. 
- **` test_utilities.f90`**: Contains routines for unit testing and validation of core modules.
- **` corr_writer.f90`**: Contains routines for storing correlation functions to CSV file.
- **` ham_hdf5_io.f90`**: Not in usage yet. Contains routines for storing Hamiltonian arrays in HDF5 format.


All files should be placed in the **same directory** before compilation.

---

## **Installation & Compilation**  

### **Prerequisites**  
Ensure you have:  
- A **Fortran compiler**(e.g., `gfortran`, Intel `ifort/ifx`, or `flang`)  
- **OpenMP** support(for parallelization)  
- **Make** installed  

### **Compiling the Code**  
To compile and link the Fortran files, simply run:  
```bash
make
```
This compiles and links modules, libraries and main file into an executable `./bin/exe`. Compiler options can be adjusted in the Makefile. Clean up build files with:
```bash
make clean
```

### **Running the Code**  

Prepare your simulation input in input.nml. Then execute:

./bin/exe

This creates a timestamped folder in output/ with results including:

- Lattice structure files
- Hamiltonian components
- Correlation data
- Spectra, eigenstates, and plots

## Plotting Energy Data

The `plot_energy.py` script allows you to **visualize energy levels as a function of V2** for a chosen V1 and simulation parameters.

**Features:**
- Automatically selects the correct `run_*` folder based on the parameters provided (`parameters.json` in each run).
- Reads `energy.csv` from `run_YYYYMMDD_HHMMSS/spectra/`.
- Saves the plot by default in `run_YYYYMMDD_HHMMSS/plots/`.
- Optional `--out` argument allows saving the plot to a custom location or filename.

---

### Usage

```bash
python src/python/plot_energy.py \
    --output-dir output \
    --filter ucx=3 ucy=3 cluster=18C irrep=A1 filling=0.5 \
    --v1 1.0 --v2min 0.0 --v2max 1.0
```
    


- `--output-dir`: Top-level folder containing `run_*` directories.
- `--filter`: Key=value pairs to select the run (e.g., `ucx=3 ucy=3 cluster=18C irrep=A1 filling=0.5`).
- `--v1`: V1 value to plot.
- `--v2min` / `--v2max`: Optional V2 range to include in the plot.
- `--out`: Optional filename to save the figure elsewhere instead of the default run folder.

### **Roadmap** ###

Add support for:
- entanglement entropy
- Level spacing statistics
- Inverse participation ratio(IPR)
- Density-density correlations
- Spin-spin correlations
