hamiltonian.mod : \
  src/hamiltonian.f90

build/hamiltonian.o : \
  src/hamiltonian.f90 build/parameters.mod build/parameters.mod build/input_variables.mod \
  build/variables.mod build/symmetries.mod build/parameters.mod \
  build/basis.mod build/file_utils.mod build/functions.mod \
  build/io_routines.mod build/printing_routines.mod

