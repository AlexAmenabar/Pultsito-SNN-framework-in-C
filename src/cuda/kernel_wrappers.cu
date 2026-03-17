#include <cuda.h>
#include <cuda_runtime.h>

#include "kernel_wrappers.cuh"
#include "networks/snn.h"
#include "simulations/results.h"

#ifdef __cplusplus
extern "C" {
#endif

// Kernel wrappers for storing them in function pointers
void initialize_neurons_cuda(GPU_SNN_t *snn, dim3 grid, dim3 block){

}

void process_input_current_cuda(GPU_SNN_t *gpu_snn, size_t iN, size_t N, size_t t, size_t gt, size_t N_per_launch, size_t launch, dim3 grid, dim3 block);
void process_V_cuda(GPU_SNN_t *gpu_snn, size_t N, dim3 grid, dim3 block);
void process_firing_cuda(GPU_SNN_t *gpu_snn, GPU_results_t *gpu_results, size_t iN, size_t N, size_t t, dim3 grid, dim3 block);

#ifdef __cplusplus
}
#endif

#endif