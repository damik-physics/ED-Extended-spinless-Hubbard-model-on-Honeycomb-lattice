io_routines.mod : \
  src/io_routines.f90

build/io_routines.o : \
  src/io_routines.f90 build/input_variables.mod build/input_variables.mod build/file_utils.mod \
  build/functions.mod

