# # ======================================
# # Makefile for structured Fortran project
# # Builds from src/fortran/, outputs to build/ & bin/
# # ======================================

# # Compiler & flags
# FC := ifx

# # MKL include directory
# MKL_INC := $(MKLROOT)/include

# # Compiler flags: add MKL includes
# FFLAGS ?= -O2 -g -check all -traceback -MMD -MF build/$*.d -I$(MKL_INC)
# # DFLAGS = -O0 -g -check all -traceback

# # Libraries: MKL + ARPACK if needed
# LIBS = -qmkl -larpack


# # Directories
# # SRC_DIR := src
# SRC_DIR := src/fortran
# BUILD_DIR := build
# BIN_DIR := bin


# # Source & object files
# # SRC := $(wildcard $(SRC_DIR)/*.f90)
# SRC := src/fortran/params.f90 src/fortran/functions.f90 src/fortran/io_utils.f90 src/fortran/basis.f90 src/fortran/lattice.f90 src/fortran/symmetries.f90 src/fortran/hamiltonian.f90 src/fortran/file_utils.f90 src/fortran/diagonalization.f90 src/fortran/main.f90
# # OBJ := $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(SRC))
# OBJ := $(patsubst %.f90,$(BUILD_DIR)/%.o,$(notdir $(SRC)))


# # Manual module dependencies
# $(BUILD_DIR)/io_utils.o: $(BUILD_DIR)/types.o $(BUILD_DIR)/functions.o $(BUILD_DIR)/file_utils.o
# $(BUILD_DIR)/symmetries.o: $(BUILD_DIR)/params.o $(BUILD_DIR)/functions.o
# $(BUILD_DIR)/basis.o: $(BUILD_DIR)/types.o $(BUILD_DIR)/params.o $(BUILD_DIR)/functions.o $(BUILD_DIR)/symmetries.o
# $(BUILD_DIR)/hamiltonian.o: $(BUILD_DIR)/params.o $(BUILD_DIR)/io_utils.o $(BUILD_DIR)/basis.o 
# $(BUILD_DIR)/lattice.o: $(BUILD_DIR)/functions.o
# $(BUILD_DIR)/file_utils.o: $(BUILD_DIR)/params.o $(BUILD_DIR)/types.o $(BUILD_DIR)/symmetries.o   
# $(BUILD_DIR)/diagonalization.o: $(BUILD_DIR)/types.o


# # Automatically detect modules for ordering
# # MODSRC := $(shell grep -il '^\s*module' $(SRC_DIR)/*.f90 | grep -v 'program')
# # Manual module source list(hardcoded)
# MODSRC := src/fortran/params.f90 src/fortran/functions.f90 src/fortran/io_utils.f90 src/fortran/basis.f90 src/fortran/lattice.f90 src/fortran/symmetries.f90 src/fortran/hamiltonian.f90 src/fortran/file_utils.f90 src/fortran/diagonalization.f90

# # MODOBJ := $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(MODSRC))
# MODOBJ := $(patsubst %.f90,$(BUILD_DIR)/%.o,$(notdir $(MODSRC)))

# NOMODSRC := $(filter-out $(MODSRC),$(SRC))

# # NOMODOBJ := $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(NOMODSRC))
# NOMODOBJ := $(patsubst %.f90,$(BUILD_DIR)/%.o,$(notdir $(NOMODSRC)))

# # Executable name
# EXE := $(BIN_DIR)/exe

# # Default target
# all: $(EXE)

# # Link executable
# $(EXE): $(MODOBJ) $(NOMODOBJ) | $(BIN_DIR)
# 	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

# # Compile .f90 -> .o
# $(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BUILD_DIR)
# 	$(FC) $(FFLAGS) -c $< -module $(BUILD_DIR) -o $@


# # Debug build
# debug: FFLAGS=$(DFLAGS)
# debug: clean all

# # Run program
# run: $(EXE)
# 	./$(EXE)

# # Create directories if not exist
# $(BUILD_DIR):
# 	mkdir -p $(BUILD_DIR)

# $(BIN_DIR):
# 	mkdir -p $(BIN_DIR)

# # Cleaning
# clean:
# 	rm -rf $(BUILD_DIR)/* $(BIN_DIR)/*

# distclean: clean
# 	rm -rf $(BUILD_DIR) $(BIN_DIR)

# .PHONY: all run clean distclean debug

# # ------------------------------------
# # Automatically include dependency files
# # ------------------------------------
# -include $(OBJ:.o=.d)

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


# Compiler flags: include MKL sources and build directory for modules
FFLAGS ?= -O2 -g -check all -traceback -MMD -MF $(BUILD_DIR)/$*.d -I$(MKL_INC) -I$(BUILD_DIR) -I$(ARPACK_DIR)/include


# Libraries: MKL + ARPACK if needed
LIBS = -qmkl 
LIBS += -L$(ARPACK_DIR)/lib -larpack

# MKL module source files to compile
MKL_MOD_SRC := $(MKL_INC)/mkl_spblas.f90 $(MKL_INC)/mkl_solvers_ee.f90

# Compile MKL module objects
MKL_MOD_OBJS := $(patsubst $(MKL_INC)/%.f90,$(BUILD_DIR)/%.o,$(MKL_MOD_SRC))

# Your Fortran source files
SRC := src/fortran/params.f90 src/fortran/functions.f90 src/fortran/io_utils.f90 src/fortran/basis.f90 src/fortran/lattice.f90 src/fortran/symmetries.f90 src/fortran/hamiltonian.f90 src/fortran/file_utils.f90 src/fortran/diagonalization.f90 src/fortran/test_utils.f90 src/fortran/main.f90

# Object files for your source files (strip src/fortran prefix)
OBJ := $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(SRC))

# Manual module dependencies
$(BUILD_DIR)/io_utils.o: $(BUILD_DIR)/types.o $(BUILD_DIR)/functions.o $(BUILD_DIR)/file_utils.o
$(BUILD_DIR)/symmetries.o: $(BUILD_DIR)/types.o $(BUILD_DIR)/params.o $(BUILD_DIR)/functions.o
$(BUILD_DIR)/basis.o: $(BUILD_DIR)/types.o $(BUILD_DIR)/params.o $(BUILD_DIR)/functions.o $(BUILD_DIR)/symmetries.o
$(BUILD_DIR)/hamiltonian.o: $(BUILD_DIR)/params.o $(BUILD_DIR)/io_utils.o $(BUILD_DIR)/basis.o 
$(BUILD_DIR)/lattice.o: $(BUILD_DIR)/types.o $(BUILD_DIR)/params.o $(BUILD_DIR)/functions.o $(BUILD_DIR)/io_utils.o
$(BUILD_DIR)/file_utils.o: $(BUILD_DIR)/params.o $(BUILD_DIR)/types.o $(BUILD_DIR)/symmetries.o   
$(BUILD_DIR)/test_utils.o: $(BUILD_DIR)/file_utils.o   
$(BUILD_DIR)/diagonalization.o: $(BUILD_DIR)/types.o $(BUILD_DIR)/functions.o $(BUILD_DIR)/params.o $(BUILD_DIR)/file_utils.o $(BUILD_DIR)/test_utils.o $(BUILD_DIR)/io_utils.o 

# Executable name
EXE := $(BIN_DIR)/exe

# Default target
all: $(EXE)

# Link executable, including MKL module objects
$(EXE): $(MKL_MOD_OBJS) $(OBJ) | $(BIN_DIR)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

# Compile your Fortran sources
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -module $(BUILD_DIR) -o $@

# Compile MKL Fortran module sources
$(BUILD_DIR)/%.o: $(MKL_INC)/%.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -module $(BUILD_DIR) -o $@

# Create directories if they don't exist
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
