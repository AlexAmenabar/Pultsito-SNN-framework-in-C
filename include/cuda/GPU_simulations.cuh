#ifndef GPU_SIMULATIONS_H
#define GPU_SIMULATIONS_H


#ifdef __cplusplus
extern "C" {
#endif


typedef struct GPU_SNN_t GPU_SNN_t;
typedef struct GPU_dataset_t GPU_dataset_t;
typedef struct simulation_configuration_t simulation_configuration_t;
typedef struct GPU_results_t GPU_results_t;
typedef struct tmp_batch_cpu_t tmp_batch_cpu_t;
typedef struct cuda_info_t cuda_info_t;

// Simulations related functions
GPU_results_t* simulate_batches_GPU(cuda_info_t *cuda_info, GPU_SNN_t *cpu_snn, GPU_dataset_t *cpu_dataset, simulation_configuration_t *conf);


// CPU2GPU comunications
GPU_SNN_t** cpy_SNN2GPU(GPU_SNN_t *cpu_snn, cuda_info_t *cuda_info, simulation_configuration_t *conf);
GPU_dataset_t** cpy_dataset2GPU(GPU_dataset_t *cpu_dataset, cuda_info_t *cuda_info);


// GPU data initializations
GPU_results_t** initialize_results_str_in_GPU(size_t nDevices, size_t n, size_t s, cuda_info_t *cuda_info);
tmp_batch_cpu_t** allocate_batch_matrix_in_GPU(size_t nDevices, size_t iN, size_t T, cuda_info_t *cuda_info);
void cpy_batch_results_GPU2CPU(GPU_results_t **cpu_results, GPU_results_t **gpu_results, cuda_info_t *cuda_info, size_t N, size_t T);


#ifdef __cplusplus
}
#endif

#endif