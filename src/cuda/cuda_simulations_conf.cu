#include "config/config_loader.h"
#include "datasets/datasets.h"
#include "networks/snn.h"
#include "simulations/results.h"

#include "cuda/cuda_simulations.cuh"
#include "cuda/cuda_simulations_conf.cuh"
#include "cuda/cuda_utils.cuh"
#include "cuda/cuda_results.cuh"

#include "cuda/priv_cuda_simulations.cuh"
#include "cuda/neuron_models/priv_cuda_lif.cuh"
#include "cuda/priv_cuda_results.cuh"

// [CUDA]

/// @brief Function to get the number of devices available for the GPU simulation
/// @param cuda_info structure with information for configuring the GPU simualtion
/// @param conf configuration file with information about the network
/// @return SNN structure initialized
size_t get_number_of_devices(cuda_info_t *cuda_info, simulation_configuration_t *conf){

    size_t available_devices, nDevices;

    available_devices = cuda_info->nDevices; // available devices
    nDevices = conf->multi_gpu; // requested number of devices

    // check whether requested devices are available, if not, correct
    if(available_devices < nDevices)
        nDevices = available_devices;

    // store the configured number of devices in the cuda information structure
    cuda_info->nDevices = nDevices;

    // return the number of devices
    return nDevices;
}

/// @brief Function to split the batch simulation between several GPUs
/// @param cuda_info structure with information for configuring the GPU simualtion
/// @param conf configuration file with information about the network
void split_batch_in_GPUs(cuda_info_t *cuda_info, simulation_configuration_t *conf){

    size_t batch_per_dev, rem_samples;
    size_t i, off = 0, nDevices, batch_size;

    // store the number of devices
    nDevices = cuda_info->nDevices;
    printf("number of devices = %zu\n", nDevices);
    fflush(stdout);
    batch_size = conf->batch_size;

    // compute batch_size per device
    batch_per_dev = batch_size / nDevices; 
    rem_samples = batch_size % nDevices; // compute how much samples do not fit on the devices

    // allocate memory for each device data
    cuda_info->dev_batch_size = (size_t*)calloc(nDevices, sizeof(size_t)); // number of samples per batch on each device
    cuda_info->dev_batch_offset = (size_t*)calloc(nDevices, sizeof(size_t)); // offset of the first sample in the device
    
    // set the batch_per_dev to each device
    for(i = 0; i<nDevices; i++)    
        cuda_info->dev_batch_size[i] = batch_per_dev;

    // update batch_per_dev using the remaining samples
    for(i = 0; i<nDevices; i++)
        if(i < rem_samples)
            cuda_info->dev_batch_size[i] ++;


    // compute offsets: the position of the first sample in the batch processed by the device
    for(i = 0; i<nDevices; i++){

        cuda_info->dev_batch_offset[i] = off; 
        off += cuda_info->dev_batch_size[i]; 
    }
}

/// @brief Function to check if the GPUs have enough memory for running the simulation
/// @param cuda_info structure with information for configuring the GPU simualtion
/// @param snn SNN structure
/// @param dataset dataset to be simulated
/// @param conf configuration file with information about the network
void check_GPU_memory(cuda_info_t *cuda_info, GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf){

    double available_memory, mem;

    // get necessary memory of several data structures
    cuda_info->network_size = get_snn_size(snn);
    cuda_info->network_cpy_size = get_snn_cpy_size(snn);
    cuda_info->dataset_size = get_dataset_size(dataset); 
    cuda_info->results_size = get_results_size(snn->n_neurons, dataset->n_samples, conf->time_steps);

    // compute available memory in the GPU without copies
    mem = 
        cuda_info->gpu_free_mem[0] - cuda_info->dataset_size - cuda_info->results_size - cuda_info->network_size - 500 * 1024 * 1024;

    // TODO: temp, if not enough memory do something instead of leaving
    // check wether there is enough memory in the GPU
    if((mem - (cuda_info->network_cpy_size * cuda_info->dev_batch_size[0])) <= 0){

        // compute maximum number of networks copies that can be allocated in the GPU
        printf(" > Not enough memory in the GPU! Reduce the batch size. Exiting\n");
        exit(1);
    }

    printf(" > Available memory for SNN copies in GPU = %.2fMB\n", mem / 1024.0 / 1024.0);
    fflush(stdout);
}

/// @brief Function to allocate several auxiliary structure for the GPUs synchronization
/// @param cuda_info structure with information for configuring the GPU simualtion
/// @param snn SNN structure
void allocate_helper_structures_for_GPU_simulation(cuda_info_t *cuda_info, GPU_SNN_t *snn){

    size_t dev;

    // allocate memory for helper structures 
    cuda_info->tmp_snn = (GPU_SNN_t**)malloc(cuda_info->nDevices * sizeof(GPU_SNN_t*)); // tmp snn structure
    cuda_info->dw = (float**)malloc(cuda_info->nDevices * sizeof(float*)); // dw [nDevices, nSynapses], used for weights sync
    
    // loop over GPUs to allocate memory for each device helper
    for(dev = 0; dev<cuda_info->nDevices; dev++){
        
        cuda_info->tmp_snn[dev] = (GPU_SNN_t*)malloc(sizeof(GPU_SNN_t));
        cuda_info->dw[dev] = (float*)malloc(snn->n_synapses * sizeof(float));
    }
}

/// @brief Function configure the cuda grids and blocks for the kernels
/// @param cuda_info structure with information for configuring the GPU simualtion
/// @param snn SNN structure
/// @param conf configuration file with information about the network
void configure_cuda_grids(cuda_info_t *cuda_info, GPU_SNN_t *snn, simulation_configuration_t *conf){
    
    size_t dev, nDevices = cuda_info->nDevices;

    // allocate memory for number of threads on each device kernel
    cuda_info->n_thr_per_blk_neurons_x      = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_thr_per_blk_neurons_y      = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_thr_per_blk_neurons_z      = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_blk_neurons_x              = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_blk_neurons_y              = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_blk_neurons_z              = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_threads_per_blk_synapses_x = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_threads_per_blk_synapses_y = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_threads_per_blk_synapses_z = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_blk_synapses_x             = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_blk_synapses_y             = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_blk_synapses_z             = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_thr_per_blk_in_neurons_x   = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_thr_per_blk_in_neurons_y   = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_thr_per_blk_in_neurons_z   = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_blk_in_neurons_x           = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_blk_in_neurons_y           = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_blk_in_neurons_z           = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_thr_per_blk_all_neurons_x  = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_thr_per_blk_all_neurons_y  = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_thr_per_blk_all_neurons_z  = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_blk_all_neurons_x          = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_blk_all_neurons_y          = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_blk_all_neurons_z          = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_thr_per_blk_uw_x           = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_thr_per_blk_uw_y           = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_thr_per_blk_uw_z           = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_blk_uw_x                   = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_blk_uw_y                   = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_blk_uw_z                   = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_thr_per_blk_is_x           = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_thr_per_blk_is_y           = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_thr_per_blk_is_z           = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_blk_is_x                   = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_blk_is_y                   = (size_t*)calloc(nDevices, sizeof(size_t));
    cuda_info->n_blk_is_z                   = (size_t*)calloc(nDevices, sizeof(size_t));

    // allocate memory for stroing the number of blocks per batch
    cuda_info->blocks_per_batch = (size_t*)calloc(nDevices, sizeof(size_t)); 

    // [TODO]: revise the comments
    // loop over devices, and set the size of the cuda grids and blocks
    for(dev = 0; dev < cuda_info->nDevices; dev ++){
        
        printf(" Dev %zu batch size = %zu / %zu | offset = %zu\n", dev, cuda_info->dev_batch_size[dev], conf->batch_size, cuda_info->dev_batch_offset[dev]);
        fflush(stdout);

        // configuration for kernels that process neurons
        cuda_info->n_thr_per_blk_neurons_x[dev] = 512;
        cuda_info->n_thr_per_blk_neurons_y[dev] = 1;
        cuda_info->n_thr_per_blk_neurons_z[dev] = 1;
        cuda_info->n_blk_neurons_x[dev] = ((snn->n_neurons * cuda_info->n_networks_per_dev[dev]) / cuda_info->n_thr_per_blk_neurons_x[dev]) + 1;
        cuda_info->n_blk_neurons_y[dev] = 1;
        cuda_info->n_blk_neurons_z[dev] = 1;

        // configuration for kernels that process synapses
        cuda_info->n_threads_per_blk_synapses_x[dev] = 512;
        cuda_info->n_threads_per_blk_synapses_y[dev] = 1;
        cuda_info->n_threads_per_blk_synapses_z[dev] = 1;
        cuda_info->n_blk_synapses_x[dev] = (snn->n_synapses * cuda_info->n_networks_per_dev[dev]) / cuda_info->n_threads_per_blk_synapses_x[dev] + 1;
        cuda_info->n_blk_synapses_y[dev] = 1;
        cuda_info->n_blk_synapses_z[dev] = 1;

        // configuration for kernels that process input currents
        cuda_info->n_thr_per_blk_in_neurons_x[dev] = 512;
        cuda_info->n_thr_per_blk_in_neurons_y[dev] = 1;
        cuda_info->n_thr_per_blk_in_neurons_z[dev] = 1;
        cuda_info->n_blk_in_neurons_x[dev] = snn->n_input_neurons * cuda_info->n_networks_per_dev[dev] / cuda_info->n_thr_per_blk_in_neurons_x[dev] + 1;
        cuda_info->n_blk_in_neurons_y[dev] = 1;
        cuda_info->n_blk_in_neurons_z[dev] = 1;

        // configuration for kernels that process neurons
        cuda_info->n_thr_per_blk_all_neurons_x[dev] = 512;
        cuda_info->n_thr_per_blk_all_neurons_y[dev] = 1;
        cuda_info->n_thr_per_blk_all_neurons_z[dev] = 1;
        cuda_info->n_blk_all_neurons_x[dev] = (snn->n_input_neurons + snn->n_neurons) * cuda_info->n_networks_per_dev[dev] / cuda_info->n_thr_per_blk_all_neurons_x[dev] + 1;
        cuda_info->n_blk_all_neurons_y[dev] = 1;
        cuda_info->n_blk_all_neurons_z[dev] = 1;
        
        // configuration for kernels that process weights update
        cuda_info->n_thr_per_blk_uw_x[dev] = 512;
        cuda_info->n_thr_per_blk_uw_y[dev] = 1;
        cuda_info->n_thr_per_blk_uw_z[dev] = 1;
        cuda_info->n_blk_uw_x[dev] = snn->n_synapses / cuda_info->n_thr_per_blk_uw_x[dev] + 1;
        cuda_info->n_blk_uw_y[dev] = 1;
        cuda_info->n_blk_uw_z[dev] = 1;

        printf(" Grids set\n");
        fflush(stdout);


        // [TODO]: all the code below should be improved
        // configure the kernel for processing the input current

        // 
        cuda_info->n_neurons_per_launch = snn->n_neurons; // all neurons in one launch
        cuda_info->n_neuron_launchs = 1; // only one launch including all neurons
  


        // slicing the loop of synapses in trhN sections, compute how much cuda blocks are necessary
        cuda_info->batch_size_per_block = 1024 / conf->thrN; // on each cuda block 1024 / thrN samples are processed, not entire batches
        if(cuda_info->dev_batch_size[dev] < cuda_info->batch_size_per_block)
            cuda_info->batch_size_per_block = cuda_info->dev_batch_size[dev];


        // since sometimes it is impossible to compute the entire batch in only one block, compute how much samples in the
        // batch will be computed inside each cuda block
        cuda_info->blocks_per_batch[dev] = cuda_info->dev_batch_size[dev] / cuda_info->batch_size_per_block; // each block processes batch_size_per_block samples,
        
        
        // so batch_size / batch_size_per_block blocks are necessary to process the entire batch. Each neuron is processed by block_per_batch
        // blocks, since each batch processes one neuron

        // compute number of threads on each cuda block
        cuda_info->n_thr_per_blk_is_x[dev] = conf->thrN * cuda_info->batch_size_per_block; // the number of divisions to the array of synapses * the batch samples
        cuda_info->n_thr_per_blk_is_y[dev] = 1;
        cuda_info->n_thr_per_blk_is_z[dev] = 1;

        // compute the number of blocks
        //cuda_info->n_blk_is_x[dev] = snn->n_neurons * cuda_info->blocks_per_batch[dev];//(snn->n_neurons * cuda_info->n_thr_per_blk_is_x[dev]) / cuda_info->n_thr_per_blk_is_x[dev] + 1;
        
        cuda_info->n_blk_is_x[dev] = cuda_info->n_neurons_per_launch * cuda_info->blocks_per_batch[dev];
        cuda_info->n_blk_is_y[dev] = 1;
        cuda_info->n_blk_is_z[dev] = 1;
        
        printf(" Device %zu: \n - dev batch size = %zu \n - blocks per batch = %zu \n - batch size per block = %zu\n", dev, cuda_info->dev_batch_size[dev], cuda_info->blocks_per_batch[dev], cuda_info->batch_size_per_block);
    }
}


// [TODO]
void compute_input_synapses_subsets(cuda_info_t *cuda_info, GPU_SNN_t *snn, simulation_configuration_t *conf){

    size_t i, j;

    // get available bytes in shared memory per each block
    double shared_memory = cuda_info->shared_memory_mem[0];

    // store max synapses to be computed per thread
    size_t max_synapses_per_thr = 100; // temporal


    // define arrays to store offsets and numbers of offsets
    size_t *n_subsets_neuron; // number of subsets of input synapses for each neuron
    size_t *off_subsets_neuron; // first subset of the neuron (if neuron[i-1] has 3 subsets, the neuron[i] subset is 3)

    size_t *n_input_synapses_subset; // number of input synapses for each subset
    size_t *off_input_synapses_subset; // number of synapses processed before this subset

    size_t n_subsets = 0; // total number of subsets (the sum of the number of subsets of each neuron)
    size_t r_subset = 0; // helper variable

    // allocate memory
    off_subsets_neuron = (size_t*)calloc(snn->n_neurons, sizeof(size_t));
    n_subsets_neuron = (size_t*)calloc(snn->n_neurons, sizeof(size_t));


    // count total number of subsets first, and set number of subsets per neuron and the offset of the neuron
    for(i = 0; i<snn->n_neurons; i++){

        // subset offset is the number of previous offsets
        off_subsets_neuron[i] = n_subsets; 

        // compute number of subsets of the neuron and the number of synapses in the last subset
        n_subsets_neuron[i] = snn->n_neuron_input_synapses[i] / max_synapses_per_thr;

        // if there are remaining subsets, add one more subset
        if((snn->n_neuron_input_synapses[i] % max_synapses_per_thr) > 0)
            n_subsets_neuron[i] += 1;

        // compute total number of subsets
        n_subsets += n_subsets_neuron[i];
    }


    // allocate memory synapses information of each offset
    off_input_synapses_subset = (size_t*)calloc(n_subsets, sizeof(size_t)); // offset of the first synapse on each subset
    n_input_synapses_subset = (size_t*)calloc(n_subsets, sizeof(size_t)); // number of synapses on each subset
    
    // control variables
    size_t next_subset = 0;
    size_t n_synapses = 0; // total number of synapses processed to be used as offset

    // loop over input synapses of the neuron to form the subsets
    for(i = 0; i<snn->n_neurons; i++){

        // all subsets until the last has the maximum number of synapses
        for(j = 0; j<n_subsets_neuron[i]-1; j++){

            off_input_synapses_subset[next_subset] = n_synapses; // set the offset of the subset in the array of synapses
            n_input_synapses_subset[next_subset] = max_synapses_per_thr; // subset has the maximum number of synapses in it

            // update control variables
            n_synapses += max_synapses_per_thr; // processed synapses
            next_subset ++; // synapse counter
        }

        // set the last subset
        off_input_synapses_subset[next_subset] = n_synapses;

        // if r_subset > 0, the last subset of the neuron has less synapses than the max 
        r_subset = snn->n_neuron_input_synapses[i] % max_synapses_per_thr; 
        if(r_subset == 0){
            
            n_input_synapses_subset[next_subset] = max_synapses_per_thr;
            n_synapses += max_synapses_per_thr;
        }
        else{
            
            n_input_synapses_subset[next_subset] = r_subset;
            n_synapses += r_subset;
        }

        next_subset ++;
    }
}

void copy_constants_to_devices(cuda_info_t *cuda_info, simulation_configuration_t *conf){

    size_t dev;
    size_t nDevices = cuda_info->nDevices;
    cudaError_t err;

    // copy constants to all gpu devices
    for(dev = 0; dev<nDevices; dev++){
        
        cudaSetDevice(dev);
        
        err = cudaMemcpyToSymbol(learn, &conf->learn, sizeof(size_t));
        err = cudaMemcpyToSymbol(batch_size, &conf->batch_size, sizeof(size_t));
        err = cudaMemcpyToSymbol(batch_size_per_block, &cuda_info->batch_size_per_block, sizeof(size_t));
        err = cudaMemcpyToSymbol(blocks_per_batch, &cuda_info->blocks_per_batch[dev], sizeof(size_t));
        err = cudaMemcpyToSymbol(dev_batch_size, &cuda_info->dev_batch_size[dev], sizeof(size_t)); // how much samples simulates each device
        err = cudaMemcpyToSymbol(dev_batch_offset, &cuda_info->dev_batch_offset[dev], sizeof(size_t)); // offset to the first sample simulated by the device in the batch
        err = cudaMemcpyToSymbol(thrN, &conf->thrN, sizeof(size_t)); // offset to the first sample simulated by the device in the batch
    }
}


// [TODO]: think about this function, very caotic
cuda_info_t* configure_cuda_simulation(GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf){

    size_t i, dev, off, nDevices, batch_size = conf->batch_size;;
    
    // initialize cuda information structure
    cuda_info_t *cuda_info = getGPUProperties();

    // store information in cuda info struct
    cuda_info->n_samples = dataset->n_samples; // store number of samples in the dataset
    cuda_info->time_steps = conf->time_steps; // time steps of the simulation
    cuda_info->batch_size = batch_size; // store batch size in cuda info
    cuda_info->chunk_size = 10; // each 10 batches move


    // compute the number of devices available for the simulation
    nDevices = get_number_of_devices(cuda_info, conf);

    // compute the number of samples for each batch will simulate each device
    split_batch_in_GPUs(cuda_info, conf);

    printf(" Checking GPU memory\n");
    fflush(stdout);   
    
    /* Get memory occupation in GPU: improve the part of n_networks */
    check_GPU_memory(cuda_info, snn, dataset, conf);
    

    // [TODO]: if there is no enough memory for copying the network batch_size times, use each copy several times
    cuda_info->n_networks_per_dev = (size_t*)calloc(nDevices, sizeof(size_t)); // number of maximum networks per dev
    for(i=0; i<cuda_info->nDevices; i++){
        
        cuda_info->n_networks_per_dev[i] = cuda_info->dev_batch_size[i];
    }

    printf(" Allocating helper structures\n");
    fflush(stdout);   
    
    // allocate memory for auxiliary structures for the GPU simulation
    allocate_helper_structures_for_GPU_simulation(cuda_info, snn);


        printf(" Configuring cuda grids\n");
    fflush(stdout);   
    // configure cuda grids and blocks for kernels
    configure_cuda_grids(cuda_info, snn, conf);


    printf(" Copying constants\n");
    fflush(stdout);    
    // copy constants from the CPU to the GPUs
    copy_constants_to_devices(cuda_info, conf);

    printf(" Constants copied\n");
    fflush(stdout);  

    // set kernel pointer for neuron model
    switch(conf->neuron_type){
        
        case 1:
            set_LIF_cuda_kernels(cuda_info);
            break;
    }
    printf(" Pointers set\n");
    fflush(stdout);  
    printf(" Initializing results structures...\n");
    fflush(stdout);

    // struct to store the results in multiple GPUs and then map to multiple CPUs
    cuda_info->gpu_results = initialize_results_str_in_GPU(nDevices, snn->n_neurons, cuda_info);
    cuda_info->intermedaite_cpu_results = initialize_results_str_in_CPU(nDevices, snn->n_neurons, cuda_info);

    // TEMP
    cuda_info->thrN = conf->thrN;
    cuda_info->cpu_snn = snn;
    cuda_info->cpu_dataset = dataset;


    // return structure 
    return cuda_info;
}

void deallocate_cuda_info_str(cuda_info_t *cuda_info){

    printf(" TODO\n");

}