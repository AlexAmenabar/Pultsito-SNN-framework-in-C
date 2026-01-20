#ifndef GPU_SIMULATIONS_H
#define GPU_SIMULATIONS_H


#ifdef __cplusplus
extern "C" {
#endif


// copy structure from CPU to GPU
GPU_SNN_t** cpy_SNN2GPU(GPU_SNN_t *cpu_snn, cuda_info_t *cuda_info, simulation_configuration_t *conf);
GPU_dataset_t** cpy_dataset2GPU(GPU_dataset_t *cpu_dataset, cuda_info_t *cuda_info);
GPU_results_t** initialize_results_str_in_GPU(size_t nDevices, size_t n, size_t s, cuda_info_t *cuda_info);
tmp_batch_cpu_t** allocate_batch_matrix_in_GPU(size_t nDevices, size_t iN, size_t T, cuda_info_t *cuda_info);
void cpy_batch_results_GPU2CPU(GPU_results_t **cpu_results, GPU_results_t **gpu_results, cuda_info_t *cuda_info, size_t N, size_t T);

GPU_results_t* simulate_batches_GPU(cuda_info_t *cuda_info, GPU_SNN_t *cpu_snn, GPU_dataset_t *cpu_dataset, simulation_configuration_t *conf);


// DEPREACTED
GPU_results_t* simulate_in_GPU_old(GPU_SNN_t **gpu_snn, GPU_dataset_t **gpu_dataset, simulation_configuration_t *conf, cuda_info_t *cuda_info, spiking_nn_t *cpu_snn, input_data_t *cpu_dataset);



#ifdef __cplusplus
}
#endif

#endif