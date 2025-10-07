# CUDA directory:
CUDA_ROOT_DIR=/usr/local/cuda


# CC compiler options:
CC=gcc
CC_FLAGS=
CC_LIBS=

# NVCC compiler options:
NVCC=nvcc
NVCC_FLAGS=
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
EXE = cuda_simulation_no_learn 

# Object files:
OBJS = $(OBJ_DIR)/main.o $(OBJ_DIR)/snn_library.o $(OBJ_DIR)/load_data.o $(OBJ_DIR)/helpers.o $(OBJ_DIR)/lif_neuron.o $(OBJ_DIR)/simulations.o $(OBJ_DIR)/stdp.o $(OBJ_DIR)/GPU_lif_neuron.o

## Compile ##
# Link c and CUDA compiled object files to target executable:
$(BIN)/$(EXE) : $(OBJS)
	$(CC) $(CC_FLAGS) $(OBJS) -o $@ $(CUDA_INC_DIR) $(CUDA_LIB_DIR) $(CUDA_LINK_LIBS)

# Compile main.c file to object files:
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.c
	$(CC) $(CC_FLAGS) -c $< -o $@ -I$(INC_DIR) -I$(INC_DIR_LIBS)

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

# Clean objects in object directory.
clean:
	$(RM) build/* *.o $(OBJS)
