# ==================================================
# Build mode (default = CPU)
# Use:
#   make cuda
#   make 
#   make avx512
# ==================================================
USE_CUDA ?= 0
USE_AVX512 ?= 0

# CUDA directory:
CUDA_ROOT_DIR=/usr/local/cuda

LCC=g++ #compiler for linking C and C++ code (CPU and CUDA)

# CC compiler options:
CC=gcc
CC_FLAGS= -O3 -fopenmp
CC_LIBS=

# NVCC compiler options:
NVCC=nvcc
NVCC_FLAGS= -O3 --ptxas-options=-v -Xcompiler "-DCUDA -fopenmp" # -fPIC
NVCC_LIBS=


# CUDA library directory:
CUDA_LIB_DIR= -L$(CUDA_ROOT_DIR)/lib64
# CUDA include directory:
CUDA_INC_DIR= -I$(CUDA_ROOT_DIR)/include 
# CUDA linking libraries:
CUDA_LINK_LIBS= -lcudart


## Project file structure ##
# Source file directory:
SRC_DIR = src/pultsito
# Object file directory:
OBJ_DIR = build
# Include header file diretory:
INC_DIR = include/pultsito
INC_DIR_LIBS = lib



## Make variables ##

# Target executable name:
BIN = bin
EXE = 

ifeq ($(USE_CUDA),1)
	EXE = cuda_simulator
else	
	EXE = cpu_simulator
endif

# Object files:
# Used by both CPU and GPU
COMMON_OBJS = \
$(OBJ_DIR)/main.o \
$(OBJ_DIR)/snn_library.o \
$(OBJ_DIR)/load_data.o \
$(OBJ_DIR)/helpers.o \
$(OBJ_DIR)/lif_neuron.o \
$(OBJ_DIR)/simulations.o \
$(OBJ_DIR)/stdp.o \

GPU_OBJS = \
$(OBJ_DIR)/GPU_lif_neuron.o \
$(OBJ_DIR)/cuda_helpers.o \
$(OBJ_DIR)/GPU_simulations.o


ifeq ($(USE_CUDA),1)
	OBJS = $(COMMON_OBJS) $(GPU_OBJS)
else	
	OBJS = $(COMMON_OBJS)
endif


ifeq ($(USE_CUDA),1)
	CC_FLAGS += -DCUDA
	CUDA_LINK = $(CUDA_INC_DIR) $(CUDA_LIB_DIR) $(CUDA_LINK_LIBS)
else
	CC_FLAGS += -ffast-math -fopenmp -funroll-loops -fprefetch-loop-arrays -flto -mtune=native

	ifeq ($(USE_AVX512),1)
		CC_FLAGS += -DAVX512
	endif

	CUDA_LINK =
endif



## Compile ##
# Link c and CUDA compiled object files to target executable:
$(BIN)/$(EXE) : $(OBJS)
	$(LCC) $(CC_FLAGS) $(OBJS) -o $@ $(CUDA_LINK) -no-pie -lm -fopenmp

# Compile main.c file to object files:
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.c
	$(CC) $(CC_FLAGS) -DREORDER -c $< -o $@ -I$(INC_DIR) -I$(INC_DIR_LIBS) 

# Compile C source files to object files:
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.c
	$(CC) $(CC_FLAGS) -c $< -o $@ -I$(INC_DIR) -I$(INC_DIR_LIBS) 

$(OBJ_DIR)/%.o : $(SRC_DIR)/neuron_models/%.c
	$(CC) $(CC_FLAGS) -c $< -o $@ -I$(INC_DIR) -I$(INC_DIR_LIBS) 

$(OBJ_DIR)/%.o : $(SRC_DIR)/simulations/%.c
	$(CC) $(CC_FLAGS) -c $< -o $@ -I$(INC_DIR) -I$(INC_DIR_LIBS) 

$(OBJ_DIR)/%.o : $(SRC_DIR)/training_rules/%.c
	$(CC) $(CC_FLAGS) -c $< -o $@ -I$(INC_DIR) -I$(INC_DIR_LIBS)

# Compile CUDA source files to object files:
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cu $(INC_DIR)/%.cuh
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@ $(NVCC_LIBS)

$(OBJ_DIR)/%.o : $(SRC_DIR)/neuron_models/%.cu 
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@ $(NVCC_LIBS) -I$(INC_DIR)

$(OBJ_DIR)/%.o : $(SRC_DIR)/cuda/%.cu 
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@ $(NVCC_LIBS) -I$(INC_DIR)



# build options: GPU (cuda), CPU (not vectorized), vectorized (AVX512)
gpu:
	$(MAKE) USE_CUDA=1
	$(MAKE) USE_AVX512=0

cpu:
	$(MAKE) USE_CUDA=0
	$(MAKE) USE_AVX512=0

avx512:
	$(MAKE) USE_CUDA=0
	$(MAKE) USE_AVX512=1

# Clean objects in object directory.
clean:
	$(RM) build/* *.o $(OBJS)