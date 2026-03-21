#ifndef GPU_SIMULATIONS_H
#define GPU_SIMULATIONS_H


#ifdef __cplusplus
extern "C" {
#endif

// forwarded declarations
typedef struct GPU_SNN_t GPU_SNN_t;
typedef struct GPU_dataset_t GPU_dataset_t;
typedef struct simulation_configuration_t simulation_configuration_t;
typedef struct GPU_results_t GPU_results_t;
typedef struct tmp_batch_cpu_t tmp_batch_cpu_t;
typedef struct cuda_info_t cuda_info_t;


/*/// [USAGE NOT RECOMMENDED]
/// @brief Function for simulating all the batches in the GPU
/// @param gpu_snn SNN structure stored in the GPU memory
/// @param gpu_dataset Dataset structure stored in the GPU memory
/// @param gpu_results Results structure stored in the GPU memory
/// @param conf Structure with the configuration information
/// @param cuda_info Structure with information for performing the simulation in the GPU
/// @param cpu_snn SNN structure stored in the CPU memory
/// @param cpu_dataset Dataset structure stored in the CPU memory
/// @return Structure with the results
GPU_results_t** simulate_batches_GPU(GPU_SNN_t **gpu_snn, GPU_dataset_t **gpu_dataset, simulation_configuration_t *conf, cuda_info_t *cuda_info);*/

/// @brief Function for simulating all the batches in the GPU
/// @param gpu_snn SNN structure stored in the GPU memory
/// @param gpu_dataset Dataset structure stored in the GPU memory
/// @param gpu_results Results structure stored in the GPU memory
/// @param conf Structure with the configuration information
/// @param cuda_info Structure with information for performing the simulation in the GPU
/// @param cpu_snn SNN structure stored in the CPU memory
/// @param cpu_dataset Dataset structure stored in the CPU memory
/// @param bidx Index of the batch to simulate
/// @return Structure with the results
void simulate_batch_GPU(GPU_results_t *cpu_results, GPU_SNN_t **gpu_snn, GPU_dataset_t **gpu_dataset, simulation_configuration_t *conf, cuda_info_t *cuda_info, size_t bidx);


#ifdef __cplusplus
}
#endif

#endif