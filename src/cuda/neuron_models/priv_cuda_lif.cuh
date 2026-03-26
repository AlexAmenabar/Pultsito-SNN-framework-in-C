#ifndef PRIV_CUDA_LIF_CUH
#define PRIV_CUDA_LIF_CUH

#ifdef __cplusplus
extern "C" {
#endif

// forwarded declarations
typedef struct GPU_SNN_t GPU_SNN_t;
typedef struct cuda_info_t cuda_info_t;
typedef struct GPU_results_t GPU_results_t;


/* Function pointers */
void set_LIF_cuda_kernels(cuda_info_t *cuda_info);

/* Memory management */
void allocate_LIF_memory_cuda(GPU_SNN_t *d_snn, size_t N, size_t batch);
void cpy_LIF_memory_cuda(GPU_SNN_t *d_snn, GPU_SNN_t *h_snn, size_t N);

/* Wrappers */
void wrap_initialize_LIF_neurons_batch_kernel(GPU_SNN_t *gpu_snn, cuda_info_t *cuda_info, size_t dev);
void wrap_process_input_current_batch_kernel(GPU_SNN_t *gpu_snn, size_t iN, size_t N, size_t t, size_t gt, cuda_info_t *cuda_info, size_t dev);
void wrap_process_V_batch_kernel(GPU_SNN_t *gpu_snn, size_t N, cuda_info_t *cuda_info, size_t dev);
void wrap_process_firing_batch_kernel(GPU_SNN_t *gpu_snn, GPU_results_t *gpu_results, size_t iN, size_t N, size_t t, size_t gt, cuda_info_t *cuda_info, size_t dev);

#ifdef __cplusplus
}
#endif

#endif