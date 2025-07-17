# **Exact Diagonalization of the Extended Hubbard Model**  

## **Overview**  
This project computes the **exact many-body spectrum and eigenstates** of the **extended Hubbard model** in **2D** on both **regular** and **tilted honeycomb clusters**. The implementation includes:  

- **Lanczos** and **full diagonalization** methods  
- Calculation of various **ground-state observables**, such as:  
  - **Current-current correlation function**  
  - **Density-density correlation function**  
  - **Charge and spin correlations**  
  <!-- Extension to include entanglement entropy, level spacing and inverse partition ratio is under work. -->
- **Symmetrization under the** $C_{6v}$ **space group** (translations, rotations, and mirror symmetries)
- **Hamiltonian diagonalization within individual irreducible representations**  
- **Parallelization with OpenMP** for improved performance  

---

## **Code Structure**  

<!-- 📁 **Project Directory**   -->
<!-- 📂 Extended-Hubbard-model/ ├── 📜 Makefile # Compilation script ├── 📜 main.f90 # Main program ├── 📜 variables.f90 # User-defined parameters ├── 📜 routines.f90 # Core subroutines and functions ├── 📂 output/ # Directory for simulation results └── 📜 README.md # This documentation file -->
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
│   └── python/            # Analysis scripts
│       └── analyze_runs.py
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
│       ├── input.nml
│       ├── lattice_data/
│       ├── logs/
│       ├── plots/
│       ├── hamiltonians/
│       ├── correlations/
│       ├── spectra/
│       ├── states/

│
├── doc/                   # Project documentation
│   └── README.md
│
└── archive/               # Legacy and experimental code
    ├── old/
    └── dev/
```
### **Modules and Functions**  

- **` main.f90`**: Program entry point
- **` basis.f90`**: Basis states and related operations
- **` diagonalization.f90`**: Diagonalization routines
- **` hamiltonian.f90`**: Hamiltonian definitions and matrix elements
- **` lattice.f90`**: Lattice geometry and connectivity
- **` symmetries.f90`**: Symmetry operations and irreps
- **` types.f90`**: Derived types and data structures
- **` functions.f90`**: Helper functions specific to the physics
- **` parameters.f90`**: Input parameters and variable definitions
- **` utils.f90`**: Utility routines: I/O, file handling, tests


All files should be placed in the **same directory** before compilation.

---

## **Installation & Compilation**  

### **Prerequisites**  
Ensure you have:  
- A **Fortran compiler** (e.g., `gfortran`, Intel `ifort/ifx`, or `flang`)  
- **OpenMP** support (for parallelization)  
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

### **Postprocessing** ###

To analyze results or generate plots, use:

```bash
python3 src/python/analyze_runs.py output/run_YYYYMMDD_HHMMSS/
```

### **Roadmap** ###

Add support for:
- entanglement entropy
- Level spacing statistics
- Inverse participation ratio (IPR)
