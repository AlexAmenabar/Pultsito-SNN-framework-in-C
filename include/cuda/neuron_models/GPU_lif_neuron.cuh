#ifndef GPU_LIF_NEURON_CUH
#define GPU_LIF_NEURON_CUH

#include <stdbool.h>
#include <cuda.h>
#include <cuda_runtime.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GPU_SNN_t GPU_SNN_t;
typedef struct GPU_dataset_t GPU_dataset_t;
typedef struct GPU_results_t GPU_results_t;
typedef struct tmp_batch_cpu_t tmp_batch_cpu_t;
typedef struct simulation_configuration_t simulation_configuration_t;
typedef struct cuda_info_t cuda_info_t;


/// @brief Simulate all the batches in the GPU
/// @param gpu_snn SNN structure in the GPU
/// @param gpu_dataset Dataset structure in the GPU
/// @param gpu_results Structure that stores the results in the GPU
/// @param gpu_tmp_batch 
void simulate_batches_LIF_GPU(GPU_SNN_t **gpu_snn, GPU_dataset_t **gpu_dataset, GPU_results_t **gpu_results, tmp_batch_cpu_t **gpu_tmp_batch, simulation_configuration_t *conf, cuda_info_t *cuda_info, GPU_SNN_t *cpu_snn, GPU_dataset_t *cpu_dataset);
void simulate_LIF_batch_multi_GPU(GPU_SNN_t **gpu_snn, GPU_dataset_t **gpu_dataset, GPU_results_t **gpu_results, tmp_batch_cpu_t **gpu_tmp_batch, simulation_configuration_t *conf, cuda_info_t *cuda_info, size_t bidx, GPU_SNN_t *cpu_snn, GPU_dataset_t *cpu_dataset, GPU_results_t **cpu_results);
void simulate_LIF_batch_single_GPU(GPU_SNN_t *gpu_snn, GPU_dataset_t *gpu_dataset, GPU_results_t *gpu_results, tmp_batch_cpu_t *gpu_tmp_batch, simulation_configuration_t *conf, cuda_info_t *cuda_info, size_t bidx, GPU_SNN_t *cpu_snn, GPU_dataset_t *cpu_dataset, GPU_results_t *cpu_results);


#ifdef __cplusplus
}
#endif

#endif