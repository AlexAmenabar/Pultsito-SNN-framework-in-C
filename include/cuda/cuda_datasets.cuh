#ifndef CUDA_DATASETS_CUH
#define CUDA_DATASETS_CUH

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GPU_dataset_t GPU_dataset_t;
typedef struct cuda_info_t cuda_info_t;

/// @brief Copy dataset structure from the CPU to the GPU(s)
/// @param cpu_dataset Dataset structure in the CPU
/// @param cuda_info Structure with information about how to perform the simulations in the GPU
/// @return Dataset structure in the GPU memory
GPU_dataset_t** cpy_dataset2GPU(GPU_dataset_t *cpu_dataset, cuda_info_t *cuda_info);


#ifdef __cplusplus
}
#endif

#endif