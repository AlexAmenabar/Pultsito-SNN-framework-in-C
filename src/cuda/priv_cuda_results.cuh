#ifndef PRIV_CUDA_RESULTS_CUH
#define PRIV_CUDA_RESULTS_CUH

#ifdef __cplusplus
extern "C" {
#endif

__global__ void reinitialize_results_number_of_spikes_GPU(GPU_results_t *results, size_t N);


void cpy_batch_results_GPU2CPU(GPU_results_t *cpu_results, GPU_results_t *gpu_results, cuda_info_t *cuda_info, simulation_configuration_t *conf, size_t N, size_t dev);


/// @brief Copy results structure from the GPU(s) to the CPU
/// @param cpu_results Results structure stored in the CPU
/// @param gpu_results Results structure stored in the GPU
/// @param cuda_info Structure with information about how to perform the simulation in the GPU
/// @param N Number of neurons
/// @param T Time steps simulated
/// @return Structure with the results stored on each GPU
void cpy_batch_results_multiGPU2CPU(GPU_results_t *cpu_results, GPU_results_t **gpu_results, cuda_info_t *cuda_info, simulation_configuration_t *conf, size_t N);


#ifdef __cplusplus
}
#endif

#endif