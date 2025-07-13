basis.mod : \
  src/basis.f90

build/basis.o : \
  src/basis.f90 build/symmetries.mod build/functions.mod build/parameters.mod

