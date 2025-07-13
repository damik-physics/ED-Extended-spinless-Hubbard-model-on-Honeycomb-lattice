# ======================================
# Makefile for structured Fortran project
# Builds from src/, outputs to build/ & bin/
# ======================================

# Compiler & flags
FC := ifx
# FFLAGS ?= -O2 -g -MMD -MF build/$*.d
FFLAGS ?= -O2 -g -check all -traceback -MMD -MF build/$*.d
# DFLAGS = -O0 -g -check all -traceback

# Libraries
# For MKL:
LIBS = -qmkl
# For OpenBLAS / LAPACK / ARPACK:
# LIBS = -larpack -llapack -lblas

# Directories
SRC_DIR := src
BUILD_DIR := build
BIN_DIR := bin


# Source & object files
# SRC := $(wildcard $(SRC_DIR)/*.f90)
SRC := src/input_variables.f90 src/variables.f90 src/parameters.f90 src/functions.f90 src/io_routines.f90 src/printing_routines.f90 src/basis.f90 src/lattice.f90 src/symmetries.f90 src/utilities.f90 src/hamiltonian.f90 src/file_utils.f90 src/main.f90
OBJ := $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(SRC))

# Manual module dependencies
$(BUILD_DIR)/io_routines.o: $(BUILD_DIR)/functions.o $(BUILD_DIR)/file_utils.o
$(BUILD_DIR)/symmetries.o: $(BUILD_DIR)/parameters.o $(BUILD_DIR)/functions.o
$(BUILD_DIR)/basis.o: $(BUILD_DIR)/parameters.o $(BUILD_DIR)/functions.o $(BUILD_DIR)/symmetries.o 
$(BUILD_DIR)/hamiltonian.o: $(BUILD_DIR)/input_variables.o $(BUILD_DIR)/variables.o $(BUILD_DIR)/parameters.o $(BUILD_DIR)/printing_routines.o $(BUILD_DIR)/io_routines.o $(BUILD_DIR)/basis.o 
$(BUILD_DIR)/lattice.o: $(BUILD_DIR)/functions.o
$(BUILD_DIR)/printing_routines.o: $(BUILD_DIR)/input_variables.o $(BUILD_DIR)/variables.o $(BUILD_DIR)/parameters.o
$(BUILD_DIR)/io_routines.o: $(BUILD_DIR)/input_variables.o $(BUILD_DIR)/functions.o 
$(BUILD_DIR)/utilities.o: $(BUILD_DIR)/variables.o $(BUILD_DIR)/functions.o


# Automatically detect modules for ordering
MODSRC := $(shell grep -il '^\s*module' $(SRC_DIR)/*.f90 | grep -v 'program')
MODOBJ := $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(MODSRC))

NOMODSRC := $(filter-out $(MODSRC),$(SRC))
NOMODOBJ := $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(NOMODSRC))

# Executable name
EXE := $(BIN_DIR)/exe

# Default target
all: $(EXE)

# Link executable
$(EXE): $(MODOBJ) $(NOMODOBJ) | $(BIN_DIR)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

# Compile .f90 -> .o
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -module $(BUILD_DIR) -o $@


# Debug build
debug: FFLAGS=$(DFLAGS)
debug: clean all

# Run program
run: $(EXE)
	./$(EXE)

# Create directories if not exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Cleaning
clean:
	rm -rf $(BUILD_DIR)/* $(BIN_DIR)/*

distclean: clean
	rm -rf $(BUILD_DIR) $(BIN_DIR)

.PHONY: all run clean distclean debug

# ------------------------------------
# Automatically include dependency files
# ------------------------------------
-include $(OBJ:.o=.d)