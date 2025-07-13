utilities.mod : \
  src/utilities.f90

build/utilities.o : \
  src/utilities.f90 build/input_variables.mod build/variables.mod build/variables.mod \
  build/functions.mod

