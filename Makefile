# ==================================================
# Build mode (default = CPU)
# Use:
#   make cuda
#   make
#   make avx512
# ==================================================

USE_CUDA ?= 0
USE_AVX512 ?= 0

CUDA_ROOT_DIR = /usr/local/cuda

# Compilers
CC   = gcc
CXX  = g++
NVCC = nvcc

# Flags
CC_FLAGS  = -g -O0 -fopenmp
CXX_FLAGS = -g -O0 -fopenmp

NVCC_FLAGS = -g -O0 --ptxas-options=-v -rdc=true -Xcompiler "-fopenmp -DCUDA"

# CUDA
CUDA_INC_DIR  = -I$(CUDA_ROOT_DIR)/include
CUDA_LIB_DIR  = -L$(CUDA_ROOT_DIR)/lib64
CUDA_LIBS     = -lcudart

# Project structure
SRC_DIR       = src
OBJ_DIR       = build
BIN_DIR       = bin
INC_DIR       = include
PRIV_INC_DIR  = src
INC_DIR_LIBS  = lib

# Executable
ifeq ($(USE_CUDA),1)
	EXE = cuda_simulator
else
	EXE = cpu_simulator
endif

# Common objects
COMMON_OBJS = \
$(OBJ_DIR)/main.o \
$(OBJ_DIR)/config_loader.o \
$(OBJ_DIR)/datasets.o \
$(OBJ_DIR)/snn.o \
$(OBJ_DIR)/neuron_models.o \
$(OBJ_DIR)/lif_neuron.o \
$(OBJ_DIR)/results.o \
$(OBJ_DIR)/simulations.o \
$(OBJ_DIR)/stdp.o \
$(OBJ_DIR)/utils.o \
$(OBJ_DIR)/snn_generator.o \
$(OBJ_DIR)/toml.o

# CUDA objects
GPU_OBJS = \
$(OBJ_DIR)/cuda_networks.o \
$(OBJ_DIR)/cuda_datasets.o \
$(OBJ_DIR)/cuda_results.o \
$(OBJ_DIR)/cuda_simulations.o \
$(OBJ_DIR)/cuda_simulations_conf.o \
$(OBJ_DIR)/cuda_lif.o \
$(OBJ_DIR)/cuda_utils.o

# Select objects
ifeq ($(USE_CUDA),1)
	OBJS = $(COMMON_OBJS) $(GPU_OBJS)
else
	OBJS = $(COMMON_OBJS)
endif

# Extra flags
ifeq ($(USE_CUDA),1)
	CC_FLAGS += -DCUDA
	CXX_FLAGS += -DCUDA
	LINKER = $(NVCC)
	LINK_FLAGS = $(NVCC_FLAGS) $(CUDA_INC_DIR) $(CUDA_LIB_DIR) $(CUDA_LIBS) -Xcompiler -fopenmp
else
	CC_FLAGS += -ffast-math -funroll-loops -fprefetch-loop-arrays -flto -mtune=native
	CXX_FLAGS += -ffast-math -funroll-loops -fprefetch-loop-arrays -flto -mtune=native
	LINKER = $(CXX)
	LINK_FLAGS = -fopenmp -lm

	ifeq ($(USE_AVX512),1)
		CC_FLAGS  += -DAVX512 -mavx512f -mavx512bw -mavx512vl -mfma -march=skylake-avx512
		CXX_FLAGS += -DAVX512 -mavx512f -mavx512bw -mavx512vl -mfma -march=skylake-avx512
	else
		CC_FLAGS  += -mavx2 -mfma
		CXX_FLAGS += -mavx2 -mfma
	endif
endif

# ==================================================
# Build rules
# ==================================================

all: $(BIN_DIR)/$(EXE)

# Link
$(BIN_DIR)/$(EXE): $(OBJS)
	mkdir -p $(BIN_DIR)
	$(LINKER) $(OBJS) -o $@ $(LINK_FLAGS)

# ==================================================
# C compilation
# ==================================================

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	mkdir -p $(OBJ_DIR)
	$(CC) $(CC_FLAGS) -c $< -o $@ -I$(INC_DIR) -I$(PRIV_INC_DIR) -I$(INC_DIR_LIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/config/%.c
	mkdir -p $(OBJ_DIR)
	$(CC) $(CC_FLAGS) -c $< -o $@ -I$(INC_DIR) -I$(PRIV_INC_DIR) -I$(INC_DIR_LIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/datasets/%.c
	mkdir -p $(OBJ_DIR)
	$(CC) $(CC_FLAGS) -c $< -o $@ -I$(INC_DIR) -I$(PRIV_INC_DIR) -I$(INC_DIR_LIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/networks/%.c
	mkdir -p $(OBJ_DIR)
	$(CC) $(CC_FLAGS) -c $< -o $@ -I$(INC_DIR) -I$(PRIV_INC_DIR) -I$(INC_DIR_LIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/neuron_models/%.c
	mkdir -p $(OBJ_DIR)
	$(CC) $(CC_FLAGS) -c $< -o $@ -I$(INC_DIR) -I$(PRIV_INC_DIR) -I$(INC_DIR_LIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/simulations/%.c
	mkdir -p $(OBJ_DIR)
	$(CC) $(CC_FLAGS) -c $< -o $@ -I$(INC_DIR) -I$(PRIV_INC_DIR) -I$(INC_DIR_LIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/training_rules/%.c
	mkdir -p $(OBJ_DIR)
	$(CC) $(CC_FLAGS) -c $< -o $@ -I$(INC_DIR) -I$(PRIV_INC_DIR) -I$(INC_DIR_LIBS)

# External lib
$(OBJ_DIR)/%.o: $(INC_DIR_LIBS)/toml_c/%.c
	mkdir -p $(OBJ_DIR)
	$(CC) $(CC_FLAGS) -c $< -o $@ -I$(INC_DIR_LIBS)/toml_c

# ==================================================
# CUDA compilation
# ==================================================

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cu
	mkdir -p $(OBJ_DIR)
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@ -I$(INC_DIR) -I$(PRIV_INC_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/cuda/%.cu
	mkdir -p $(OBJ_DIR)
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@ -I$(INC_DIR) -I$(PRIV_INC_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/cuda/neuron_models/%.cu
	mkdir -p $(OBJ_DIR)
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@ -I$(INC_DIR) -I$(PRIV_INC_DIR)

# ==================================================
# Build modes
# ==================================================

cuda:
	$(MAKE) USE_CUDA=1 USE_AVX512=0

cpu:
	$(MAKE) USE_CUDA=0 USE_AVX512=0

avx512:
	$(MAKE) USE_CUDA=0 USE_AVX512=1

# ==================================================
# Clean
# ==================================================

clean:
	rm -rf $(OBJ_DIR)/*.o $(BIN_DIR)/*