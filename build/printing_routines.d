printing_routines.mod : \
  src/printing_routines.f90

build/printing_routines.o : \
  src/printing_routines.f90 build/parameters.mod build/variables.mod build/input_variables.mod

