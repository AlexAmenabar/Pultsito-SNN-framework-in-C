#ifndef CUDA_NETWORKS_CUH
#define CUDA_NETWORKS_CUH

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GPU_SNN_t GPU_SNN_t;
typedef struct simulation_configuration_t simulation_configuration_t;
typedef struct cuda_info_t cuda_info_t;

/// @brief Copy SNN structure from the CPU to the GPU(s)
/// @param cpu_snn SNN structure in the CPU memory
/// @param cuda_info Structure with information about how to perform the simulations in the GPU
/// @param conf Structure with configuration information
/// @return SNN structure in the GPU memory
GPU_SNN_t** cpy_SNN2GPU(GPU_SNN_t *cpu_snn, cuda_info_t *cuda_info, simulation_configuration_t *conf);


#ifdef __cplusplus
}
#endif
#endif