#ifndef PRIV_CUDA_SIMULATIONS_CUH
#define PRIV_CUDA_SIMULATIONS_CUH

// forwarded structure declarations
typedef struct GPU_SNN_t GPU_SNN_t;
typedef struct GPU_results_t GPU_results_t;
typedef struct tmp_batch_cpu_t tmp_batch_cpu_t;
typedef struct cuda_info_t cuda_info_t;
typedef struct simulation_configuration_t simulation_configuration_t;

// constant values for helping the simulations
#ifdef __CUDACC__
extern __constant__ size_t learn;
extern __constant__ size_t batch_size_per_block; // samples processed on each cuda block per neuron
extern __constant__ size_t blocks_per_batch; // number of cuda blocks for processing all neuron batches
extern __constant__ size_t dev_batch_size; // batch size proccessed by each device
extern __constant__ size_t dev_batch_offset; // offset to the sample proccessed by each device
extern __constant__ size_t batch_size;
extern __constant__ size_t thrN;
#endif


#ifdef __cplusplus
extern "C" {
#endif


void simulate_batch_single_GPU(GPU_SNN_t *gpu_snn, GPU_dataset_t *gpu_dataset, tmp_batch_cpu_t *gpu_tmp_batch, simulation_configuration_t *conf, cuda_info_t *cuda_info, size_t bidx);
void update_weights_single_GPU(GPU_SNN_t *gpu_snn, simulation_configuration_t *conf, cuda_info_t *cuda_info);
void simulate_batch_multi_GPU(GPU_SNN_t **gpu_snn, GPU_dataset_t **gpu_dataset, tmp_batch_cpu_t **gpu_tmp_batch, simulation_configuration_t *conf, cuda_info_t *cuda_info, size_t bidx);
void update_weights_multi_GPU(GPU_SNN_t **gpu_snn, simulation_configuration_t *conf, cuda_info_t *cuda_info);
void run_simulation_batch_GPU(GPU_SNN_t *gpu_snn, GPU_dataset_t *gpu_dataset, tmp_batch_cpu_t *gpu_tmp_batch, simulation_configuration_t *conf, cuda_info_t *cuda_info, size_t bidx, size_t dev);
tmp_batch_cpu_t** allocate_batch_matrix_in_GPU(size_t nDevices, size_t iN, size_t T, cuda_info_t *cuda_info);
void deallocate_batch_matrix_in_GPU(tmp_batch_cpu_t **batch_matrix, cuda_info_t *cuda_info);


/* Kernels */

// initialize synapses before starting the simulation



__global__ void initialize_synapses_batch(GPU_SNN_t *snn);
__global__ void initialize_batch_matrix_GPU(tmp_batch_cpu_t *tmp_batch, GPU_dataset_t *dataset, size_t bidx, size_t iN, size_t T);
__global__ void load_batch_time_step_in_SNN_GPU(GPU_SNN_t *gpu_snn, tmp_batch_cpu_t *tmp_batch, size_t iN, size_t N, size_t lt, size_t t, size_t bidx);
__global__ void compute_pre_fired_batch(GPU_SNN_t *gpu_snn, size_t N, size_t iN, size_t S, size_t t, size_t gt);
__global__ void compute_post_traces_batch(GPU_SNN_t *gpu_snn, size_t N);
__global__ void reinit_post_fired_batch(GPU_SNN_t *gpu_snn, size_t N);
__global__ void trace_based_STDP_batch(GPU_SNN_t *gpu_snn, size_t S);
__global__ void update_weights_batch(GPU_SNN_t *gpu_snn, size_t S);
__global__ void acc_weights_batch(GPU_SNN_t *gpu_snn, size_t S);
__global__ void update_weights_dev(GPU_SNN_t *gpu_snn, size_t S);

#ifdef __cplusplus
}
#endif

#endif