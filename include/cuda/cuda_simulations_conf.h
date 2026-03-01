#ifndef CUDA_SIMULATIONS_CONF_H
#define CUDA_SIMULATIONS_CONF_H

/* Structures used for simulating networks */

// forwarded declarations
typedef struct GPU_SNN_t GPU_SNN_t;
typedef struct GPU_dataset_t GPU_dataset_t;
typedef struct simulation_configuration_t simulation_configuration_t;

/// @brief Structure that stores the information of the kernels to launch in CUDA
typedef struct {

    // devices data
    int multi_gpu_allowed; // wether multilple devices should be used
    size_t nDevices; // number of available devices
    double *gpu_total_mem; // free memory in the GPU
    double *gpu_free_mem; // memory available in the GPU
    double *shared_memory_mem; // Bytes in shared memory
    size_t maxThreadsPerMultiprocessor;
    size_t nMultiprocessor;
    size_t maxThreads;

    // cuda simulation details (network size, dataset size...)
    double dataset_size; // dataset size
    double network_size; // network size
    double network_cpy_size; // size of each network copy
    double results_size; // results size

    double *neuron_size; // size related to each neuron
    double *batch_neuron_size; // size related to the entire neuron batch 
    double *synapse_size; // size related to each synapse
    double *batch_synapse_size; // size related to the entire synapse batch

    size_t n_networks; // DEPRECATED // n copies of the network in the GPU
    size_t n_samples; // n samples in the dataset
    size_t batch_size; // batch_size of the simulation
    size_t time_steps; // time steps of the simulation
    
    size_t batch_size_per_block;
    size_t *blocks_per_batch;
    size_t n_neuron_launchs;
    size_t n_neurons_per_launch;
    size_t *thrN_per_launch;

    
    // multigpu control variables
    size_t *dev_batch_size; // [nDevices] part of the batch simulated by each device
    size_t *dev_batch_offset; // [nDevices] the first sample in the batch to be simulated by device 
    size_t *n_networks_per_dev; // [nDevices]: number of copies on each device
    //size_t thrN;

    size_t gpuId; // used to index the correct information on each GPU
    
    
    // helpers for weight updating
    GPU_SNN_t **tmp_snn; // tmp snn structure
    float **dw; // dw [nDevices, nSynapses]     


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

///
#ifdef __cplusplus
extern "C" {
#endif

/// @brief Function to configure the simulation in cuda: batch size processed by GPU, grid...
/// @param cuda_info structure to store the information for carrying out the simulation
/// @param snn structure storing the snn
/// @param dataset structure storing the dataset
/// @param conf structure storing the configuration of the simulation
void configure_cuda_simulation(cuda_info_t *cuda_info, GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf);


#ifdef __cplusplus
}
#endif

#endif