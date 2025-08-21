# ======================================
# Makefile for structured Fortran project
# ======================================

# Compiler & flags
FC := ifx

# Environment variables
MKL_INC := $(MKLROOT)/include

# Build directories
SRC_DIR := src/fortran
BUILD_DIR := build
BIN_DIR := bin
ARPACK_DIR ?= /usr/local
HDF5_DIR ?= /opt/hdf5_ifx

# Compiler flags: include MKL sources, build dir, ARPACK, and HDF5
# FFLAGS ?= -O2 -g -check all -traceback -MMD -MF $(BUILD_DIR)/$*.d \
           -I$(MKL_INC) \
           -I$(BUILD_DIR) \
           -I$(ARPACK_DIR)/include \
           -I$(HDF5_DIR)/include
FFLAGS ?= -g -O0 -check all -traceback -fpe0 -debug full -MMD -MF $(BUILD_DIR)/$*.d \
           -I$(MKL_INC) \
           -I$(BUILD_DIR) \
           -I$(ARPACK_DIR)/include \
           -I$(HDF5_DIR)/include

# Libraries: MKL + ARPACK + HDF5
LIBS = -qmkl \
       -L$(ARPACK_DIR)/lib -larpack \
       -L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5 \
       -Wl,-rpath,$(HDF5_DIR)/lib

# MKL module source files to compile
MKL_MOD_SRC := $(MKL_INC)/mkl_spblas.f90 $(MKL_INC)/mkl_solvers_ee.f90
MKL_MOD_OBJS := $(patsubst $(MKL_INC)/%.f90,$(BUILD_DIR)/%.o,$(MKL_MOD_SRC))

# Your Fortran source files
SRC := src/fortran/params.f90 \
       src/fortran/functions.f90 \
       src/fortran/io_utils.f90 \
       src/fortran/corr_writer.f90 \
       src/fortran/basis.f90 \
       src/fortran/lattice.f90 \
       src/fortran/symmetries.f90 \
       src/fortran/ham_hdf5_io.f90 \
       src/fortran/hamiltonian.f90 \
       src/fortran/file_utils.f90 \
       src/fortran/diagonalization.f90 \
       src/fortran/observables.f90 \
       src/fortran/test_utils.f90 \
       src/fortran/main.f90

OBJ := $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(SRC))

# Manual module dependencies
$(BUILD_DIR)/io_utils.o: $(BUILD_DIR)/types.o $(BUILD_DIR)/functions.o $(BUILD_DIR)/file_utils.o
$(BUILD_DIR)/symmetries.o: $(BUILD_DIR)/types.o $(BUILD_DIR)/params.o $(BUILD_DIR)/functions.o
$(BUILD_DIR)/basis.o: $(BUILD_DIR)/types.o $(BUILD_DIR)/params.o $(BUILD_DIR)/functions.o $(BUILD_DIR)/symmetries.o
$(BUILD_DIR)/hamiltonian.o: $(BUILD_DIR)/params.o $(BUILD_DIR)/io_utils.o $(BUILD_DIR)/ham_hdf5_io.o $(BUILD_DIR)/basis.o 
$(BUILD_DIR)/lattice.o: $(BUILD_DIR)/types.o $(BUILD_DIR)/params.o $(BUILD_DIR)/functions.o $(BUILD_DIR)/io_utils.o
$(BUILD_DIR)/file_utils.o: $(BUILD_DIR)/params.o $(BUILD_DIR)/types.o $(BUILD_DIR)/symmetries.o   
$(BUILD_DIR)/test_utils.o: $(BUILD_DIR)/file_utils.o   
$(BUILD_DIR)/diagonalization.o: $(BUILD_DIR)/types.o $(BUILD_DIR)/functions.o $(BUILD_DIR)/params.o $(BUILD_DIR)/file_utils.o $(BUILD_DIR)/test_utils.o $(BUILD_DIR)/io_utils.o 
$(BUILD_DIR)/observables.o: $(BUILD_DIR)/types.o $(BUILD_DIR)/functions.o $(BUILD_DIR)/file_utils.o $(BUILD_DIR)/corr_writer.o $(BUILD_DIR)/symmetries.o

# Executable
EXE := $(BIN_DIR)/exe

# Default target
all: $(EXE)

# Link executable
$(EXE): $(MKL_MOD_OBJS) $(OBJ) | $(BIN_DIR)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

# Compile your Fortran sources
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -module $(BUILD_DIR) -o $@

# Compile MKL Fortran module sources
$(BUILD_DIR)/%.o: $(MKL_INC)/%.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -module $(BUILD_DIR) -o $@

# Create directories
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Run program
run: $(EXE)
	./$(EXE)

# Clean
clean:
	rm -rf $(BUILD_DIR)/* $(BIN_DIR)/*

distclean: clean
	rm -rf $(BUILD_DIR) $(BIN_DIR)

.PHONY: all clean distclean
