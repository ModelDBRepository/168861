CC=g++ # GNU Compiler
NVCC=nvcc # CUDA Compiler
MKLFLAG=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core # Intel MKL libraries flags
NFLAGS= -Xcompiler -Wall -O3 -arch=sm_21 -x cu -lineinfo # General flags 
CFLAGS=-Wall -Wno-write-strings -O3 -g -DDEBUG=0  # -DDEBUG=X. This variable sets the printing level (0 or 1) of the simulation
LDFLAGS=-lGL -lGLU -lglut -lpthread # OpenGL flags
NLDFLAGS=-lcublas -lcusparse # CUDA flags

# Folder paths
INCLUDE_DIR=$(shell pwd)/inc
SOURCE_DIR=$(shell pwd)/src
BUILD_DIR=$(shell pwd)/build
BIN_DIR=$(shell pwd)/bin
DOC_DIR=$(shell pwd)/doc

# Files type
vpath %.h $(INCLUDE_DIR)
vpath %.cpp $(SOURCE_DIR)
vpath %.o $(BUILD_DIR)

# Common objects for both implementations
OBJECTS=graphics.o RgbImage.o input.o init.o coupling.o mechanical_model.o space.o discretization.o configuration.o
# Specific objects for the parallel implementations
OBJECTS_GPU=neurite-gpu.o solverEx-gpu.o solverIm-gpu.o cudaException.o
# Specific objects for the sequential implementations 
OBJECTS_CPU=neurite-cpu.o solverEx-cpu.o solverIm-cpu.o

# Name of the executable file
EXECUTABLE=Neurite

# Default objectives 
OBJECTIVES= cpu gpu

# General construction rules
# Default
def: $(OBJECTIVES)
# Parallel
gpu: $(OBJECTS) $(OBJECTS_GPU)
	cd $(BUILD_DIR); $(NVCC) $^ $(LDFLAGS) $(NLDFLAGS) -o $(BIN_DIR)/$(EXECUTABLE)_solver-gpu
# Sequential
cpu: $(OBJECTS) $(OBJECTS_CPU)
	cd $(BUILD_DIR); $(CC) $^ $(LDFLAGS) $(MKLFLAG) -o $(BIN_DIR)/$(EXECUTABLE)_solver-cpu
# Previous objects. Parallel
cudaException.o: cudaException.cpp cudaException.h 
	$(NVCC) $(NFLAGS) -I$(INCLUDE_DIR) -c $< -o $(BUILD_DIR)/$@
solverEx-gpu.o: solverEx-gpu.cpp solver.h solverEx-gpu.h
	$(NVCC) $(NFLAGS) -I$(INCLUDE_DIR) -c $< -o $(BUILD_DIR)/$@
solverIm-gpu.o: solverIm-gpu.cpp solver.h solverIm-gpu.h
	$(NVCC) $(NFLAGS) -I$(INCLUDE_DIR) -c $< -o $(BUILD_DIR)/$@
neurite-gpu.o: neurite.cpp neurite.h
	$(NVCC) $(NFLAGS) -DGPU -I$(INCLUDE_DIR) -c $< -o $(BUILD_DIR)/$@
# Previous objects. Sequential
solverEx-cpu.o: solverEx-cpu.cpp solver.h solverEx-cpu.h
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) -c $< -o $(BUILD_DIR)/$@
solverIm-cpu.o: solverIm-cpu.cpp solver.h solverIm-cpu.h
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) -c $< -o $(BUILD_DIR)/$@
neurite-cpu.o: neurite.cpp neurite.h
	$(CC) $(CFLAGS) -DCPU -I$(INCLUDE_DIR) -c $< -o $(BUILD_DIR)/$@
# Prevous objects. Common
%.o: %.cpp %.h	
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) -c $< -o $(BUILD_DIR)/$@

# Cleaning rules
.PHONY: clean clean-all 
clean-all: clean-doc clean
# Cleanning rules for the documentation
clean-doc: 
	-cd $(DOC_DIR); rm -rf html latex
# General cleaning rules
clean:
	-cd $(BUILD_DIR); rm -f $(OBJECTS) $(OBJECTS_OLD) $(OBJECTS_GPU2) $(OBJECTS_GPU) $(OBJECTS_CPU) $(OBJECTS_BASE)
	-cd $(BIN_DIR); rm -f  $(EXECUTABLE)_solver-cpu  $(EXECUTABLE)_solver-gpu

# Documentation rules
documentation:
	cd $(DOC_DIR); doxygen 
