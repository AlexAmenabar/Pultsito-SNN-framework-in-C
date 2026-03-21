#ifndef PRIV_CUDA_RESULTS_CUH
#define PRIV_CUDA_RESULTS_CUH

#ifdef __cplusplus
extern "C" {
#endif

__global__ void reinitialize_results_number_of_spikes_GPU(GPU_results_t *results, size_t N);


void cpy_batch_results_struct_GPU2CPU(GPU_results_t *cpu_results, GPU_results_t *gpu_results, cuda_info_t *cuda_info, size_t N, size_t T);


#ifdef __cplusplus
}
#endif

#endif