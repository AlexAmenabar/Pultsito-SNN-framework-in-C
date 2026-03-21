#ifndef CUDA_SIMULATIONS_CONF_H
#define CUDA_SIMULATIONS_CONF_H

/* Structures used for simulating networks */

#ifdef __cplusplus
extern "C" {
#endif

// forwarded declarations
typedef struct GPU_SNN_t GPU_SNN_t;
typedef struct GPU_dataset_t GPU_dataset_t;
typedef struct GPU_results_t GPU_results_t;
typedef struct simulation_configuration_t simulation_configuration_t;

/// @brief Structure that stores the information of the kernels to launch in CUDA
typedef struct cuda_info_t {

    // GPUs data: physical properties
    int multi_gpu_allowed; // wether multilple devices should be used
    size_t nDevices; // number of available devices
    double *gpu_total_mem; // free memory in the GPU
    double *gpu_free_mem; // memory available in the GPU
    double *shared_memory_mem; // Bytes in shared memory
    size_t maxThreadsPerMultiprocessor;
    size_t nMultiprocessor;
    size_t maxThreads;

    // cuda simulation details (network size, dataset size...)
    double dataset_size; // dataset size (bytes)
    double network_size; // network size (bytes)
    double network_cpy_size; // size of each network copy (bytes)
    double results_size; // results size (bytes)

    double *neuron_size; // size related to each neuron (bytes)
    double *batch_neuron_size; // size related to the entire neuron batch (bytes) 
    double *synapse_size; // size related to each synapse
    double *batch_synapse_size; // size related to the entire synapse batch

    size_t n_networks; // number of network copies in the GPU
    size_t n_samples; // number of samples in the dataset
    size_t batch_size; // batch_size of the simulation
    size_t time_steps; // time steps of the simulation
    
    size_t batch_size_per_block; // section of the batch processed on each cuda block
    size_t *blocks_per_batch; // number of blocks to simulate each batch
    size_t n_neuron_launchs; // number of kernels necessary to process all neurons
    size_t n_neurons_per_launch; // number of neurons processed on each launch
    size_t *thrN_per_launch; // thrN value for each launch

    // multigpu control variables
    size_t *dev_batch_size; // [nDevices]: number of samples processed by each device in the batch
    size_t *dev_batch_offset; // [nDevices]: index of the first sample proccessed by each device in the batch 
    size_t *n_networks_per_dev; // [nDevices]: number of networks stored on each device (limits the parallelization)

    size_t gpuId; // used to index the correct information on each GPU
    size_t thrN;

    // Function pointers
    void (*allocate_neuron_memory_cuda)(GPU_SNN_t *d_snn, size_t N, size_t batch);
    void (*cpy_neurons_CPU2GPU_cuda)(GPU_SNN_t *d_snn, GPU_SNN_t *h_snn, size_t N);

    void (*initialize_neurons_cuda)(GPU_SNN_t *snn, cuda_info_t *cuda_info, size_t dev);
    void (*process_input_current_cuda)(GPU_SNN_t *gpu_snn, size_t iN, size_t N, size_t t, size_t gt, cuda_info_t *cuda_info, size_t dev);
    void (*process_dynamics_cuda)(GPU_SNN_t *gpu_snn, size_t N, cuda_info_t *cuda_info, size_t dev);
    void (*process_firing_cuda)(GPU_SNN_t *gpu_snn, GPU_results_t *gpu_results, size_t iN, size_t N, size_t t, cuda_info_t *cuda_info, size_t dev);
    
    // helpers for weight updating
    size_t chunk_size; // chunk size used to manage how batches to process before moving to the CPU
    GPU_results_t **gpu_results; // [nDevices] used for storing results of the GPU 
    GPU_results_t **intermedaite_cpu_results; // [nDevices] used as intermediate step between GPU distributed results and CPU accumulated results
    GPU_SNN_t **tmp_snn; // tmp snn structure
    float **dw; // dw [nDevices, nSynapses]     


    // CPU data
    GPU_SNN_t *cpu_snn;
    GPU_dataset_t *cpu_dataset;


    // number of threads and blocks per kernel
    size_t *n_thr_per_blk_neurons_x, *n_thr_per_blk_neurons_y, *n_thr_per_blk_neurons_z;
    size_t *n_thr_per_blk_synapses_x, *n_thr_per_blk_synapses_y, *n_thr_per_blk_synapses_z;
    size_t *n_thr_per_blk_in_neurons_x, *n_thr_per_blk_in_neurons_y, *n_thr_per_blk_in_neurons_z;
    size_t *n_thr_per_blk_all_neurons_x, *n_thr_per_blk_all_neurons_y, *n_thr_per_blk_all_neurons_z;
    size_t *n_thr_per_blk_uw_x, *n_thr_per_blk_uw_y, *n_thr_per_blk_uw_z;
    size_t *n_thr_per_blk_is_x, *n_thr_per_blk_is_y, *n_thr_per_blk_is_z;

    size_t *n_blk_neurons_x, *n_blk_neurons_y, *n_blk_neurons_z;
    size_t *n_blk_synapses_x, *n_blk_synapses_y, *n_blk_synapses_z;
    size_t *n_blk_in_neurons_x, *n_blk_in_neurons_y, *n_blk_in_neurons_z;
    size_t *n_blk_all_neurons_x, *n_blk_all_neurons_y, *n_blk_all_neurons_z;
    size_t *n_blk_uw_x, *n_blk_uw_y, *n_blk_uw_z;
    size_t *n_blk_is_x, *n_blk_is_y, *n_blk_is_z; 

    size_t *n_threads_per_blk_synapses_x, *n_threads_per_blk_synapses_y, *n_threads_per_blk_synapses_z;

} cuda_info_t;



/* General function to network initialization */


/// @brief Function to configure the simulation in cuda: batch size processed by GPU, grid...
/// @param cuda_info structure to store the information for carrying out the simulation
/// @param snn structure storing the snn
/// @param dataset structure storing the dataset
/// @param conf structure storing the configuration of the simulation
cuda_info_t* configure_cuda_simulation(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf);


#ifdef __cplusplus
}
#endif

#endif