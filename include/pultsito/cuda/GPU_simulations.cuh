#ifndef GPU_SIMULATIONS_H
#define GPU_SIMULATIONS_H


#ifdef __cplusplus
extern "C" {
#endif

// copy structure from CPU to GPU
GPU_SNN_t* cpy_SNN2GPU(GPU_SNN_t *cpu_snn, cuda_info_t *cuda_info);
GPU_dataset_t* cpy_dataset2GPU(GPU_dataset_t *cpu_dataset, cuda_info_t *cuda_info);
GPU_results_t* initialize_results_str_in_GPU(simulation_configuration_t *conf, int n, int s);

GPU_results_t* simulate_in_GPU(GPU_SNN_t *gpu_snn, GPU_dataset_t *gpu_dataset, simulation_configuration_t *conf, cuda_info_t *cuda_info, spiking_nn_t *cpu_snn, input_data_t *cpu_dataset);


#ifdef __cplusplus
}
#endif

#endif