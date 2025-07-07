SHELL := /bin/bash

FC = ifx
CURR_DIR = $(shell pwd)
FILES = ${CURR_DIR}/files_local
LIBP2 = /home/david/Coding/Fortran

# Compiler flags
FFLAGS = -g -qopenmp -O0 -traceback -debug all -heap-arrays -warn nointerfaces
OPTIONS = -fpp

# MKL environment (make sure MKLROOT is set correctly in your environment)
MKLROOT ?= $(shell echo $$MKLROOT)

# MKL include flags
COPTIONS = -I"${MKLROOT}/include"

# MKL link flags for ifx with dynamic linking (from your mkl_link_tool output)
LIBMKL = -L${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

# Source files
SRC = ${FILES}/module_tmi.f90 ${FILES}/variables_tmi.f90 ${FILES}/tmi.f90

# Object files
OBJ = module_tmi.o variables_tmi.o tmi.o spblas.o solver.o

# MKL Fortran modules (compiled objects)
MKL_MOD1 = spblas.o
MKL_MOD2 = solver.o

.PHONY: all clean link compile

all: link

link: compile
	$(FC) $(FFLAGS) $(OPTIONS) -o ./exe $(OBJ) -L${LIBP2} -larpack ${LIBMKL}

compile: ${SRC} ${MKLROOT}/include/mkl_spblas.f90 ${MKLROOT}/include/mkl_solvers_ee.f90
	# Compile MKL modules first
	$(FC) $(FFLAGS) $(OPTIONS) ${COPTIONS} -c ${MKLROOT}/include/mkl_spblas.f90 -o spblas.o
	$(FC) $(FFLAGS) $(OPTIONS) ${COPTIONS} -c ${MKLROOT}/include/mkl_solvers_ee.f90 -o solver.o
	# Compile user sources
	$(FC) $(FFLAGS) $(OPTIONS) ${COPTIONS} -c ${FILES}/variables_tmi.f90 -o variables_tmi.o
	$(FC) $(FFLAGS) $(OPTIONS) ${COPTIONS} -c ${FILES}/module_tmi.f90 -o module_tmi.o
	$(FC) $(FFLAGS) $(OPTIONS) ${COPTIONS} -c ${FILES}/tmi.f90 -o tmi.o

clean:
	rm -f *.mod *.o exe
