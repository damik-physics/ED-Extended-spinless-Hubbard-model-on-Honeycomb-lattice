# Simplified Honeycomb Hubbard Model

This is a streamlined implementation of the honeycomb Hubbard model that solves the same core physics problem as the original complex codebase, but with a much simpler structure.

## Problem Description

The **Honeycomb Hubbard Model** studies quantum many-body systems of interacting electrons on a 2D honeycomb lattice (like graphene). The physics involves:

- **Hopping (t)**: Electrons can move between neighboring sites
- **Interactions (V1)**: Electrons repel each other when on neighboring sites  
- **Many-body quantum mechanics**: The full problem requires diagonalizing exponentially large matrices

## Physics Equation

The Hamiltonian is:
```
H = -t Σ (c†ᵢcⱼ + h.c.) + V1 Σ nᵢnⱼ
```

Where:
- First term: electron hopping between nearest neighbors
- Second term: density-density interactions
- The sums run over nearest-neighbor pairs on the honeycomb lattice

## Simplified vs Original Code

### Original Complex Implementation:
- **10+ modules** with thousands of lines
- **Symmetry handling**: Complex momentum sectors, irreducible representations
- **Scalability**: Sparse matrices, parallel computing, HDF5 I/O
- **Multiple algorithms**: ARPACK, MKL, FEAST diagonalization
- **Advanced features**: Disorder, current correlations, parameter sweeps

### This Simplified Version:
- **2 files**: Main program + utilities (< 300 lines total)
- **Core physics only**: Hopping + interactions on small lattices
- **Dense matrices**: Direct LAPACK diagonalization
- **Educational focus**: Clear, readable implementation of the physics

## Key Simplifications Made

1. **Small system sizes** (4×4 lattice) instead of arbitrary large clusters
2. **Dense matrices** instead of sparse matrix optimizations  
3. **Single-threaded** instead of OpenMP parallelization
4. **Direct diagonalization** instead of iterative eigensolvers
5. **No symmetries** - work in full Hilbert space
6. **Fixed parameters** instead of parameter sweeps
7. **Minimal I/O** instead of comprehensive data management

## Usage

```bash
# Compile
make

# Run  
make run

# Clean
make clean
```

## Understanding the Output

The program calculates:
- **Ground state energy**: Lowest eigenvalue of the many-body Hamiltonian
- **Energy per site**: Ground state energy divided by number of lattice sites
- **Energy gap**: Difference between ground and first excited state
- **Low-lying spectrum**: First few eigenvalues

## Physics Regimes

- **V1/t < 0.5**: Weakly interacting, metallic behavior
- **V1/t ~ 1**: Intermediate coupling, competing physics
- **V1/t > 2**: Strongly interacting, potential insulating behavior
- **Half filling**: Special case that may show Mott insulator physics

## Extensions

To make this more like the original codebase, you could add:

1. **Parameter sweeps**: Loop over different V1 values
2. **Larger lattices**: Increase Lx, Ly (but exponential scaling!)
3. **Sparse matrices**: Use CSR format for larger systems
4. **Symmetries**: Implement translation invariance to reduce Hilbert space
5. **Observables**: Calculate density correlations, current-current correlations
6. **Next-nearest neighbors**: Add V2 interactions
7. **Disorder**: Random on-site potentials

## Mathematical Background

The computational challenge is that the Hilbert space dimension grows as:
```
dim = C(N_sites, N_particles) 
```

For the honeycomb lattice with N_sites = 2×Lx×Ly sites and N_particles electrons.

For example:
- 4×4 lattice (32 sites), half-filled: ~600 million basis states!
- This is why the original code uses sophisticated sparse matrix techniques

This simplified version works for small systems where we can afford dense matrix storage and diagonalization.

## Dependencies

- **gfortran**: GNU Fortran compiler
- **LAPACK/BLAS**: Linear algebra libraries

Install on Ubuntu/Debian:
```bash
sudo apt install gfortran liblapack-dev libblas-dev
```

Install on macOS:
```bash
brew install gcc lapack
```