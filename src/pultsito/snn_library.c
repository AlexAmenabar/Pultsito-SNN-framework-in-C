#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "snn_library.h"
#include "training_rules/stdp.h"
#include "load_data.h"
#include "neuron_models/neuron_models.h"

/* [PRIVATE] */


/* [Helpers] */

// Function to get the maximum delay value in the network
size_t get_max_delay(GPU_SNN_t *snn){

    size_t i = 0;

    size_t max_delay = 0;
    for(i=0; i<snn->n_synapses; i++){

        if(snn->delay[i] > max_delay){
            max_delay = snn->delay[i];
        }
    }

    // set max delay
    snn->LT = max_delay;

    return max_delay;
}


// Function to set the function pointers for the neuron type in the network
void set_network_neuron_ptr(GPU_SNN_t *snn, simulation_configuration_t *conf){

    switch(conf->neuron_type){
        case 1: // LIF
            snn->neuron_allocator = &allocate_memory_for_LIF_neurons; // memory allocation for LIF neurons
            snn->init_neurons = &initialize_LIF_neurons; // LIF neurons initialization
            snn->reinit_neurons = &reinitialize_LIF_neurons_batch; // LIF neurons reinitializations
            snn->neuron_dynamics = &compute_LIF_V_batch; // LIF neuron dynamics computation
            break;
        default: // LIF
            snn->neuron_allocator = &allocate_memory_for_LIF_neurons; // memory allocation for LIF neurons
            snn->init_neurons = &initialize_LIF_neurons; // LIF neurons initialization
            snn->reinit_neurons = &reinitialize_LIF_neurons_batch; // LIF neurons reinitializations
            snn->neuron_dynamics = &compute_LIF_V_batch; // LIF neuron dynamics computation
            break;      
    }
}

// [Memory allocations]

void allocate_memory_for_SNN(GPU_SNN_t *snn, simulation_configuration_t *conf){

    // allocate memory for neurons and synapses
    snn->neuron_allocator(snn, conf);
    
    // allocate memory for synapses
    snn->w                             = (float*)malloc(snn->n_synapses * conf->batch_size * sizeof(float));
    snn->init_w                        = (float*)malloc(snn->n_synapses * conf->batch_size * sizeof(float));
    snn->dw                            = (float*)malloc(snn->n_synapses * conf->batch_size * sizeof(float));
    snn->delay                         = (int*)malloc(snn->n_synapses * conf->batch_size * sizeof(int));
    snn->lr                            = (int*)malloc(snn->n_synapses * conf->batch_size * sizeof(int));
    snn->pre_fired                     = (char*)malloc(snn->n_synapses * conf->batch_size * sizeof(char));
    snn->pre_trace                     = (float*)malloc(snn->n_synapses * conf->batch_size * sizeof(float));
    snn->pre_neuron_index              = (size_t*)malloc(snn->n_synapses * conf->batch_size * sizeof(size_t));
    snn->post_neuron_index             = (size_t*)malloc(snn->n_synapses * conf->batch_size * sizeof(size_t));
    

    // allocate memory for spike matrix spk matrix
    snn->spk_matrix                    = (char*)malloc((snn->n_neurons + snn->n_input_neurons) * snn->LT * conf->batch_size * sizeof(char)); // time steps is temporal
}

void reallocate_spk_matrix(GPU_SNN_t *snn, size_t N, size_t B, size_t T){

    free(snn->spk_matrix);
    snn->spk_matrix = (char*)calloc(N * B * T, sizeof(char));
}

// [Memory deallocations]

void deallocate_arrays_str(network_construction_lists_t *data){

    // deallocate internal arrays
    if(data->v_list)
        free(data->v_list);
    if(data->v_thres_list)
        free(data->v_thres_list);
    if(data->v_rest_list)
        free(data->v_rest_list);
    if(data->R_list)
        free(data->R_list);
    if(data->r_time_list)
        free(data->r_time_list);

    if(data->weight_list)
        free(data->weight_list);
    if(data->training_zones)
        free(data->training_zones);
    if(data->delay_list)
        free(data->delay_list);
    
    for(size_t i = 0; i<data->n_neurons; i++){
        if(data->synaptic_connections[i])
            free(data->synaptic_connections[i]);
    }
    if(data->synaptic_connections)
        free(data->synaptic_connections);

    // deallocate struct
    if(data)
        free(data);
}

void deallocate_snn_str(GPU_SNN_t *snn){

    // deallocate internal arrays
    if(snn->v)
        free(snn->v);
    if(snn->v_thresh)
        free(snn->v_thresh);
    if(snn->v_rest)
        free(snn->v_rest);
    if(snn->r_period)
        free(snn->r_period);
    if(snn->r_period_remain)
        free(snn->r_period_remain);
    if(snn->res)
        free(snn->res);
    if(snn->post_fired)
        free(snn->post_fired);
    if(snn->post_trace)
        free(snn->post_trace);
    if(snn->arrI)
        free(snn->arrI);
    if(snn->n_neuron_input_synapses)
        free(snn->n_neuron_input_synapses);
    if(snn->neuron_input_synapses_offset)
        free(snn->neuron_input_synapses_offset);
    
    if(snn->w)
        free(snn->w);
    if(snn->init_w)
        free(snn->init_w);
    if(snn->dw)
        free(snn->dw);
    if(snn->delay)
        free(snn->delay);
    if(snn->lr)
        free(snn->lr);
    if(snn->pre_neuron_index)
        free(snn->pre_neuron_index);
    if(snn->post_neuron_index)
        free(snn->post_neuron_index);
    if(snn->pre_fired)
        free(snn->pre_fired);
    if(snn->pre_trace)
        free(snn->pre_trace);

    if(snn->spk_matrix)
        free(snn->spk_matrix);

    // deallocate struct
    if(snn)
        free(snn);
}

void deallocate_dataset_str(GPU_dataset_t *dataset){
    
    // deallocate internal arrays
    if(dataset->n_spikes_per_feature)
        free(dataset->n_spikes_per_feature);
    if(dataset->sample_offset)
        free(dataset->sample_offset);
    if(dataset->feature_offset)
        free(dataset->feature_offset);
    if(dataset->spikes)
        free(dataset->spikes);
    if(dataset->freq)
        free(dataset->freq);
    if(dataset->first_spk)
        free(dataset->first_spk);

    // deallocate struct
    if(dataset)
        free(dataset);
}

void deallocate_results_str(GPU_results_t *results){

    if(results->n_spks)
        free(results->n_spks);
    if(results->gnt_spks)
        free(results->gnt_spks);

    if(results)
        free(results);
}


// [Network initialization]

/// @brief Initialize synaptic arrays
/// @param snn SNN structure to initialize synapses in
/// @param data structure with the values for initializing the synapses
/// @param conf structure with configuration information
void initialize_synapses_CPU(GPU_SNN_t *snn, network_construction_lists_t *data, simulation_configuration_t *conf){

    size_t i;

    for(i = 0; i<snn->n_synapses; i++){

        // read synapse parameters from lists
        snn->w[i] = data->weight_list[i];
        snn->init_w[i] = snn->w[i];
        snn->dw[i] = 0.0; 
        snn->delay[i] = data->delay_list[i];
        snn->lr[i] = data->training_zones[i]; // [TODO]: this will call the function in an array of function pointers

        // initialized later
        snn->pre_neuron_index[i] = 0;
        snn->post_neuron_index[i] = 0;
        
        // control variables
        snn->pre_fired[i] = 0;
        snn->pre_trace[i] = 0.0;
    }
}

/// @brief Compute the input synapses offsets for each neuron
/// @param snn SNN structure
/// @param data Structure containing helper values
void count_neurons_input_synapses_and_conpute_offsets(GPU_SNN_t *snn, network_construction_lists_t *data){
    
    size_t i, j, off = 0;

    // loop over neurons and count number of input synapses
    for(i = 0; i<snn->n_neurons; i++){

        // loop over input synapses of the neuron
        for(j = 0; j<(size_t)data->synaptic_connections[i][0]; j++){

            // set the number of input synapses
            // [i * 2 + 1] is the index of the neuron
            // [i * 2 + 2] the number of synapses from that neuron to the current 
            snn->n_neuron_input_synapses[i] += (size_t)data->synaptic_connections[i][j*2+2];
        }
    }

    // loop over neurons and compute offsets
    for(i = 0; i<snn->n_neurons; i++){

        // set offset
        snn->neuron_input_synapses_offset[i] = off;

        // update offset
        off += snn->n_neuron_input_synapses[i];
    }
}

/// @brief Set the input and output neurons of the synapses
/// @param snn SNN structure
/// @param data Structure containing helper values
void connect_synapses_to_network(GPU_SNN_t *snn, network_construction_lists_t *data){

    size_t i, j, l, next_syn = 0;

    // loop over the neurons
    for(i = 0; i<snn->n_neurons; i++){

        // loop over the input synapses of the neuron
        for(j = 0; j<(size_t)data->synaptic_connections[i][0]; j++){

            // set the input and output neurons of the synapse
            for(l = 0; l<(size_t)data->synaptic_connections[i][j * 2 + 2]; l++){
                
                snn->post_neuron_index[next_syn] = i + snn->n_input_neurons; // first [n_input_neurons] neurons are virtual input
                snn->pre_neuron_index[next_syn] = (size_t)data->synaptic_connections[i][j * 2 + 1];
                next_syn ++;
            }
        }
    }
}

/// @brief Connect the network following an input criteria, where, for each neuron, its input synapses are stored
/// @param snn SNN structure
/// @param data Structure containing helper values
/// @param conf structure containing the configuration of the simulation
void connect_network_input_criteria(GPU_SNN_t *snn, network_construction_lists_t *data, simulation_configuration_t *conf){

    // count the number of input synapses and the offsets for each neuron
    count_neurons_input_synapses_and_conpute_offsets(snn, data);

    // set pre and post neurons indexes for each synapse
    connect_synapses_to_network(snn, data);
}


// [CUDA]

/// @brief Function to get the number of devices available for the GPU simulation
/// @param cuda_info structure with information for configuring the GPU simualtion
/// @param conf configuration file with information about the network
/// @return SNN structure initialized
size_t get_number_of_devices(cuda_info_t *cuda_info, simulation_configuration_t *conf){

    size_t n_max_devices;

    // check wether multiGPU simulations are activated (0, no; 1: all devices; n > 1: number of devices)
    cuda_info->multi_gpu_allowed = conf->multi_gpu >= 1;
    
    // compute the number of GPUs available for the simulation
    if(conf->multi_gpu > 1){

        n_max_devices = conf->multi_gpu; // the number of GPUs indicated in the configuration file
    }
    else if(cuda_info->multi_gpu_allowed == 1){ 

        n_max_devices = cuda_info->nDevices; // use all GPUs available in the allocation
    }
    else{ 

        n_max_devices = 1; // no multiGPU
    }

    // we cannot use more than the available devices, correct number if necessary
    if(n_max_devices > cuda_info->nDevices){

        n_max_devices = cuda_info->nDevices;
        printf(" > There are no %zu devices available, %zu devices will be used.\n", n_max_devices, cuda_info->nDevices);
        fflush(stdout);

        if(n_max_devices == 0){
            
            printf(" >> No GPUs available!! Exitting!\n");
            fflush(stdout);
            exit(1);
        }
    }

    // store the configured number of devices in the cuda information structure
    cuda_info->nDevices = n_max_devices;

    // return the number of devices
    return n_max_devices;
}

/// @brief Function to split the batch simulation between several GPUs
/// @param cuda_info structure with information for configuring the GPU simualtion
/// @param conf configuration file with information about the network
void split_batch_in_GPUs(cuda_info_t *cuda_info, simulation_configuration_t *conf){

    size_t batch_per_dev, rem_samples;
    size_t i, off = 0, nDevices, batch_size;

    // store the number of devices
    nDevices = cuda_info->nDevices;
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


/* [PUBLIC FUNCTIONS] */

GPU_SNN_t* initialize_network_cpu(simulation_configuration_t *conf){

    // allocate memory for SNN structure
    GPU_SNN_t *snn = (GPU_SNN_t*)malloc(sizeof(GPU_SNN_t));

    // load network information into intermedaite arrays
    network_construction_lists_t *data = load_network_information_in_lists(conf);

    // initialize general info
    snn->n_neurons = data->n_neurons;
    snn->n_input_neurons = data->n_input_neurons;
    snn->n_output_neurons = data->n_output_neurons;
    snn->n_synapses = data->n_synapses;
    snn->LT = data->max_delay;

    // set function pointers for managing neurons
    set_network_neuron_ptr(snn, conf);

    // allocate memory for the SNN structure arrays
    allocate_memory_for_SNN(snn, conf);

    // initialize neurons (using the previously set function pointer)
    snn->init_neurons(snn, data, conf);
    
    // initialize synapses
    initialize_synapses_CPU(snn, data, conf);

    // connect the synapses and the neurons
    connect_network_input_criteria(snn, data, conf);

    // deallocate intermeadite arrays structure used for initializing the neurons and synapses
    deallocate_arrays_str(data);

    // return the initialized SNN structure
    return snn;
}


void cpy_snn(GPU_SNN_t *snn, simulation_configuration_t *conf){

    size_t i, j, t, p, B, LT, N, iN, S;

    // temporal variables to store values
    float *tmp_v, *tmp_arrI, *tmp_w, *tmp_dw, *tmp_pre_trace, *tmp_post_trace;
    int *tmp_rperiod_remain;
    char *tmp_pstfired, *tmp_prefired, *tmp_spk_matrix;

    // [TODO]: improve
    size_t padding = 32; // 32 bytes (256 bits) of padding for copied variables (avoid segfault in vectorized code)

    B = conf->batch_size;
    LT = snn->LT;
    N = snn->n_neurons;
    iN = snn->n_input_neurons;
    S = snn->n_synapses;

    // store all values for copying later
    tmp_v                         = snn->v;
    tmp_rperiod_remain            = snn->r_period_remain;
    tmp_pstfired                  = snn->post_fired;
    tmp_post_trace                = snn->post_trace;
    tmp_arrI                      = snn->arrI;

    tmp_w                         = snn->w;
    tmp_dw                        = snn->dw;
    tmp_pre_trace                 = snn->pre_trace;
    tmp_prefired                  = snn->pre_fired;

    tmp_spk_matrix                = snn->spk_matrix;

    // allocate memory for copied variables
    snn->v                            = (float*)malloc(N * B * sizeof(float) + padding);
    snn->arrI                         = (float*)malloc(N * B * sizeof(float) + padding);
    snn->r_period_remain              = (int*)malloc(N * B * sizeof(int) + padding);
    snn->post_fired                   = (char*)malloc(N * B * sizeof(char) + padding);
    snn->post_trace                   = (float*)malloc(N * B * sizeof(float) + padding);

    snn->w                            = (float*)malloc(S * B * sizeof(float) + padding);
    snn->dw                           = (float*)malloc(S * B * sizeof(float) + padding);
    snn->pre_trace                    = (float*)malloc(S * B * sizeof(float) + padding);
    snn->pre_fired                    = (char*)malloc(S * B * sizeof(char) + padding); 
    
    // reallocate spk matrix
    snn->spk_matrix                   = (char*)malloc((iN + N) * B * LT * sizeof(char) + padding);


    // initialize neurons using old values
    for(i = 0; i<N; i++){

        for(j = 0; j<B; j++){

            // copy values
            snn->v[i * B + j] = tmp_v[i];
            snn->r_period_remain[i * B + j] = tmp_rperiod_remain[i];
            snn->post_fired[i * B + j] = tmp_pstfired[i];
            snn->post_trace[i * B + j] = tmp_post_trace[i];
            snn->arrI[i * B + j] = tmp_arrI[i];
        }
    }

    // initialize synapses using old values
    for(i = 0; i<S; i++){

        for(j = 0; j<B; j++){

            // copy values
            snn->w[i * B + j] = tmp_w[i];
            snn->dw[i * B + j] = tmp_dw[i];
            snn->pre_trace[i * B + j] = tmp_pre_trace[i];
            snn->pre_fired[i * B + j] = tmp_prefired[i];
        }
    }

    // initialize spk matrix to 0
    for(t = 0; t<LT; t++){
        
        for(i = 0; i<iN + N; i++){
            
            for(j = 0; j<B; j++){

                snn->spk_matrix[(iN + N) * B * t + B * i + j] = 0;
            }
        }
    }


    // initialize padding positions to 0
    for(p = N * B ; p<padding / sizeof(float); p++){

        snn->v[p] = 0.0;                       
        snn->arrI[p] = 0.0;                         
        snn->post_trace[p] = 0.0;               
    }
    for(p = N * B ; p<padding / sizeof(int); p++){

        snn->r_period_remain[p] = 0;              
    }
    for(p = N * B ; p<padding / sizeof(char); p++){

        snn->post_fired[p] = 0;     
    }
    for(p = S * B ; p<padding / sizeof(float); p++){

        snn->w[p] = 0.0;                    
        snn->dw[p] = 0.0;  
        snn->pre_trace[p] = 0.0;                
    }
    for(p = S * B ; p<padding / sizeof(char); p++){

        snn->pre_fired[p] = 0;
    }

    for(p = (iN + N) * B * LT ; p<padding / sizeof(char); p++){

        snn->spk_matrix[p] = 0;
    }

    // deallocate original arrays
    if(tmp_v)
        free(tmp_v);
    if(tmp_rperiod_remain)
        free(tmp_rperiod_remain);
    if(tmp_pstfired)
        free(tmp_pstfired);
    if(tmp_post_trace)
        free(tmp_post_trace);
    if(tmp_arrI)
        free(tmp_arrI);
    if(tmp_w)
        free(tmp_w);
    if(tmp_dw)
        free(tmp_dw);
    if(tmp_pre_trace)
        free(tmp_pre_trace);
    if(tmp_prefired)
        free(tmp_prefired);
    if(tmp_spk_matrix)
        free(tmp_spk_matrix);
}

// [TODO]: rethink how to generalize for storing any result type
GPU_results_t* initialize_batch_results_cpu(simulation_configuration_t *conf, size_t N, size_t batch_size, size_t T){

    // allocate memory for the results structure
    GPU_results_t *results = (GPU_results_t*)malloc(sizeof(GPU_results_t));

    // allocate memory for storing the number of spikes // [TODO]: if...
    results->n_spks = (int*)calloc(N * batch_size, sizeof(int));


    // initialize execution times
    results->t = 0.0; // total execution time
    results->t_in = 0.0; // time processing input spikes
    results->t_v = 0.0; // time processing neurons dynamics
    results->t_out = 0.0; // time processing output spikes
    results->t_learn = 0.0; // time learning
    results->t_reinit = 0.0; // network reinitialization
    results->t_load = 0.0; // loading sample or batch in network

    // return results structure
    return results;
}

// [TODO]: the same mentioned for the previous one
void reinitialize_batch_results_cpu(GPU_results_t *results, simulation_configuration_t *conf, size_t N, size_t batch_size, size_t T){

    // reinitialize the number of spikes generated by each neuron
    for(size_t i = 0; i<N; i++){
        
        for(size_t j = 0; j<batch_size; j++){

            results->n_spks[i * batch_size + j] = 0;
        }
    }

    // reinitialize execution times
    results->t        = 0.0; // total execution time
    results->t_in     = 0.0; // time processing input spikes
    results->t_v      = 0.0; // time processing neurons dynamics
    results->t_out    = 0.0; // time processing output spikes
    results->t_learn  = 0.0; // time learning
    results->t_reinit = 0.0; // network reinitialization
    results->t_load   = 0.0; // loading sample or batch in network
}


void acc_batch_execution_times(GPU_results_t *results, GPU_results_t *batch_results){

    results->t        += batch_results->t;        // total execution time
    results->t_in     += batch_results->t_in;     // time processing input spikes
    results->t_v      += batch_results->t_v;      // time processing neurons dynamics
    results->t_out    += batch_results->t_out;    // time processing output spikes
    results->t_learn  += batch_results->t_learn;  // time learning
    results->t_reinit += batch_results->t_reinit; // network reinitialization
    results->t_load   += batch_results->t_load;   // loading sample or batch in network
}


// [TODO]: think about this function, very caotic
void configure_cuda_simulation(cuda_info_t *cuda_info, GPU_SNN_t *snn, GPU_dataset_t *dataset, simulation_configuration_t *conf){

    size_t i, dev, off, nDevices, batch_size = conf->batch_size;;
    
    // store information in cuda info struct
    cuda_info->n_samples = dataset->n_samples; // store number of samples in the dataset
    cuda_info->time_steps = conf->time_steps; // time steps of the simulation
    cuda_info->batch_size = batch_size; // store batch size in cuda info

    // compute the number of devices available for the simulation
    nDevices = get_number_of_devices(cuda_info, conf);

    // compute the number of samples for each batch will simulate each device
    split_batch_in_GPUs(cuda_info, conf);

    /* Get memory occupation in GPU: improve the part of n_networks */
    check_GPU_memory(cuda_info, snn, dataset, conf);
    

    // [TODO]: if there is no enough memory for copying the network batch_size times, use each copy several times
    /*cuda_info->n_networks_per_dev = (size_t*)calloc(nDevices, sizeof(size_t)); // number of maximum networks per dev
    for(i=0; i<cuda_info->nDevices; i++){
        
        cuda_info->n_networks_per_dev[i] = cuda_info->dev_batch_size[i];
    }*/

    // allocate memory for auxiliary structures for the GPU simulation
    allocate_helper_structures_for_GPU_simulation(cuda_info, snn);


    // configure cuda grids and blocks for kernels
    configure_cuda_grids(cuda_info, snn, conf);
}



// [TODO]: move?
/* Functions to get the size of several structures */

double get_snn_size(GPU_SNN_t *snn){
    
    size_t N = snn->n_neurons;
    size_t iN = snn->n_input_neurons;
    size_t S = snn->n_synapses;
    size_t LT = snn->LT;

    // n_input neurons or synapses?
    return (
        (sizeof(GPU_SNN_t) +
        sizeof(size_t) * 10 + // scalars
        
        // neurons 
        sizeof(float) * 5 * N + // neurons floats
        sizeof(int) * 3 * N + // neurons integers
        sizeof(char) * 1 * N + // neurons chars
        sizeof(size_t) * 2 * N + // neurons size_t
        
        // synapses
        sizeof(float) * 4 * S + // synapses floats
        sizeof(int) * 2 * S + // synapses integers
        sizeof(char) * 1 * S + // synapses chars
        sizeof(size_t) * 2 * S + // synapses size_t
        
        // spk matrix
        sizeof(char) * (iN + N) * LT) / 8.0
    );
}

double get_snn_cpy_size(GPU_SNN_t *snn){

    size_t N = snn->n_neurons;
    size_t iN = snn->n_input_neurons;
    size_t S = snn->n_synapses;
    size_t LT = snn->LT;

    return (

        // neurons
        (sizeof(float) * 3 * N + // copied float arrays
        sizeof(int) * 1 * N + // copied int arrays
        sizeof(char) * 1 * N + // copied char arrays

        // synapses
        sizeof(float) * 3 * S + // copied float arrays 
        sizeof(char) * 1 * S + // copied char arrays
        
        // spk matrix
        sizeof(char) * LT * (iN + N)) / 8.0
    );
}

double get_dataset_size(GPU_dataset_t *dataset){
    
    size_t nS = dataset->n_samples;
    size_t nF = dataset->n_features;
    size_t nSpks = dataset->n_spikes;

    return (

        // general and scalars
        (sizeof(GPU_dataset_t) +
        sizeof(int) + 
        sizeof(size_t) * 4 +

        // arrays
        sizeof(size_t) * (nS + nS * nF * 4 + nSpks)) / 8.0
    );
}

double get_results_size(size_t N, size_t nS, size_t T){

    return (

        (sizeof(GPU_results_t) + 
        sizeof(int) * N * nS + // number of spike generated by each neuron
        sizeof(double) * 7 // variables to store execution times
        ) / 8.0
    );
}

double get_tmp_batch_size(size_t iN, size_t nS, size_t T){

    return (

        (sizeof(tmp_batch_cpu_t) + 
        sizeof(size_t) * iN * nS * T 
        ) / 8.0
    );
}

///[TODO]: implement
double get_neuron_size(){
    
    double size = 0.0;
    // get memory related to one neuron
    

    return size;
}



// [TODO]: move?
/* Functions for printing information */

void print_network(GPU_SNN_t *snn){

    size_t i;

    printf(" > Printing network:\n\n");

    printf(" > > Network general information:\n");
    printf(" > >>  N: %zu\n", snn->n_neurons);
    printf(" > >> iN: %zu\n", snn->n_input_neurons);
    printf(" > >> oN: %zu\n", snn->n_output_neurons);
    printf(" > >>  S: %zu\n", snn->n_synapses);
    //printf(" > >> nN: %zu\n\n", snn->n_networks);

    printf(" > > Neurons data:\n");

    printf(" > >> V: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("%f, ", snn->v[i]);
    }
    printf("%f]\n", snn->v[i]);

    printf(" > >> Thresholds: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("%f, ", snn->v_thresh[i]);
    }
    printf("%f]\n", snn->v_thresh[i]);

    printf(" > >> Rests: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("%f, ", snn->v_rest[i]);
    }
    printf("%f]\n", snn->v_rest[i]);

    printf(" > >> Refract: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("%d (%d), ", snn->r_period[i], snn->r_period_remain[i]);
    }
    printf("%d (%d)]\n", snn->r_period[i], snn->r_period_remain[i]);

    printf(" > >> R: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("%d, ", snn->res[i]);
    }
    printf("%d]\n", snn->res[i]);

    printf(" > >> Post fired: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("%d, ", snn->post_fired[i]);
    }
    printf("%d]\n", snn->post_fired[i]);

    printf(" > >> I: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("%f, ", snn->arrI[i]);
    }
    printf("%f]\n", snn->arrI[i]);

    printf(" > >> N input synapses: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("%zu, ", snn->n_neuron_input_synapses[i]);
    }
    printf("%zu]\n", snn->n_neuron_input_synapses[i]);

    printf(" > >> Input synapses off: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("%zu, ", snn->neuron_input_synapses_offset[i]);
    }
    printf("%zu]\n", snn->neuron_input_synapses_offset[i]);



    printf(" > > Synapses data:\n");

    printf(" > >> W: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("%f, ", snn->w[i]);
    }
    printf("%f]\n", snn->w[i]);

    printf(" > >> Init W: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("%f, ", snn->init_w[i]);
    }
    printf("%f]\n", snn->init_w[i]);

    printf(" > >> dW: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("%f, ", snn->dw[i]);
    }
    printf("%f]\n", snn->dw[i]);

    printf(" > >> Delay: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("%d, ", snn->delay[i]);
    }
    printf("%d]\n", snn->delay[i]);

    printf(" > >> Lr: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("%d, ", snn->lr[i]);
    }
    printf("%d]\n", snn->lr[i]);

    printf(" > >> Pre neuron: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("%zu, ", snn->pre_neuron_index[i]);
    }
    printf("%zu]\n", snn->pre_neuron_index[i]);

    printf(" > >> Post neuron: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("%zu, ", snn->post_neuron_index[i]);
    }
    printf("%zu]\n", snn->post_neuron_index[i]);

    printf(" > >> Pre fired: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("%d, ", snn->pre_fired[i]);
    }
    printf("%d]\n", snn->pre_fired[i]);

    printf(" > >> Pre trace: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("%f, ", snn->pre_trace[i]);
    }
    printf("%f]\n", snn->pre_trace[i]);

    printf(" > >> Post trace: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("%f, ", snn->post_trace[i]);
    }
    printf("%f]\n", snn->post_trace[i]);

    fflush(stdout);
}

void print_networks(GPU_SNN_t *snn, simulation_configuration_t *conf){

    size_t i, b, B;
    B = conf->batch_size;

    printf(" > Printing network:\n\n");

    printf(" > > Network general information:\n");
    printf(" > >>  N: %zu\n", snn->n_neurons);
    printf(" > >> iN: %zu\n", snn->n_input_neurons);
    printf(" > >> oN: %zu\n", snn->n_output_neurons);
    printf(" > >>  S: %zu\n", snn->n_synapses);
    //printf(" > >> nN: %zu\n\n", snn->n_networks);

    printf(" > > Neurons data:\n");

    printf(" > >> V: [");
    for(i = 0; i<snn->n_neurons-1; i++){
        
        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%f, ", snn->v[i * B + b]);

        printf("%f]", snn->v[i * B + B - 1]);
    }
    printf("%f]\n", snn->v[i * B + B-1]);

    printf(" > >> Thresholds: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%f, ", snn->v_thresh[i * B + b]);

        printf("%f]", snn->v_thresh[i * B + B - 1]);    
    }
    printf("%f]\n", snn->v_thresh[i * B + B - 1]);

    printf(" > >> Rests: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%f, ", snn->v_rest[i * B + b]);

        printf("%f]", snn->v_rest[i * B + B - 1]);      
    }
    printf("%f]\n", snn->v_rest[i * B + B - 1]);

    printf(" > >> Refract: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%d (%d), ", snn->r_period[i * B + b], snn->r_period_remain[i * B + b]);

        printf("%d (%d)]", snn->r_period[i * B + B-1], snn->r_period_remain[i * B + B-1]);      
    }
    printf("%d (%d)]\n", snn->r_period[i * B + B-1], snn->r_period_remain[i * B + B-1]);

    printf(" > >> R: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%d, ", snn->res[i * B + b]);

        printf("%d]", snn->res[i * B + B-1]);      
    }
    printf("%d]\n", snn->res[i * B + B-1]);

    printf(" > >> Post fired: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%d, ", snn->post_fired[i * B + b]);

        printf("%d]", snn->post_fired[i * B + B-1]);      
    }
    printf("%d]\n", snn->post_fired[i * B + B-1]);

    printf(" > >> I: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%f, ", snn->arrI[i * B + b]);

        printf("%f]", snn->arrI[i * B + B-1]);      
    }
    printf("%f]\n", snn->arrI[i * B + B-1]);

    printf(" > >> N input synapses: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%zu, ", snn->n_neuron_input_synapses[i * B + b]);

        printf("%zu]", snn->n_neuron_input_synapses[i * B + B-1]);      
    }
    printf("%zu]\n", snn->n_neuron_input_synapses[i * B + B-1]);

    printf(" > >> Input synapses off: [");
    for(i = 0; i<snn->n_neurons-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%zu, ", snn->neuron_input_synapses_offset[i * B + b]);

        printf("%zu]", snn->neuron_input_synapses_offset[i * B + B-1]);      
    }
    printf("%zu]\n", snn->neuron_input_synapses_offset[i * B + B-1]);



    printf(" > > Synapses data:\n");

    printf(" > >> W: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%f, ", snn->w[i * B + b]);

        printf("%f]", snn->w[i * B + B-1]);      
    }
    printf("%f]\n", snn->w[i * B + B-1]);

    printf(" > >> Init W: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%f, ", snn->init_w[i * B + b]);

        printf("%f]", snn->init_w[i * B + B-1]);      
    }
    printf("%f]\n", snn->init_w[i * B + B-1]);

    printf(" > >> dW: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%f, ", snn->dw[i * B + b]);

        printf("%f]", snn->dw[i * B + B-1]);      
    }
    printf("%f]\n", snn->dw[i * B + B-1]);

    printf(" > >> Delay: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%d, ", snn->delay[i * B + b]);

        printf("%d]", snn->delay[i * B + B-1]);      
    }
    printf("%d]\n", snn->delay[i * B + B-1]);

    printf(" > >> Lr: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%d, ", snn->lr[i * B + b]);

        printf("%d]", snn->lr[i * B + B-1]);      
    }
    printf("%d]\n", snn->lr[i * B + B-1]);

    printf(" > >> Pre neuron: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%zu, ", snn->pre_neuron_index[i * B + b]);

        printf("%zu]", snn->pre_neuron_index[i * B + B-1]);      
    }
    printf("[");
    for(b = 0; b<B-1; b++)
        printf("%zu, ", snn->pre_neuron_index[i * B + b]);
    printf("%zu]\n", snn->pre_neuron_index[i * B + B-1]);

    printf(" > >> Post neuron: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%zu, ", snn->post_neuron_index[i * B + b]);

        printf("%zu]", snn->post_neuron_index[i * B + B-1]);      
    }
    printf("%zu]\n", snn->post_neuron_index[i * B + B-1]);

    printf(" > >> Pre fired: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%d, ", snn->pre_fired[i * B + b]);

        printf("%d]", snn->pre_fired[i * B + B-1]);      
    }
    printf("%d]\n", snn->pre_fired[i * B + B-1]);

    printf(" > >> Pre trace: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%f, ", snn->pre_trace[i * B + b]);

        printf("%f]", snn->pre_trace[i * B + B-1]);      
    }
    printf("%f]\n", snn->pre_trace[i * B + B-1]);

    printf(" > >> Post trace: [");
    for(i = 0; i<snn->n_synapses-1; i++){

        printf("[");
        for(b = 0; b<B-1; b++)
            printf("%f, ", snn->post_trace[i * B + b]);

        printf("%f]", snn->post_trace[i * B + B-1]);      
    }
    printf("%f]\n", snn->post_trace[i * B + B-1]);

    fflush(stdout);
}

void print_dataset(GPU_dataset_t *dataset){

    size_t i, j, l, next = 0;

    printf(" > Printing dataset: \n");
    for(i = 0; i<dataset->n_samples; i++){

        printf(" > > Sample %zu (offset = %zu)\n", i, dataset->sample_offset[i]);

        // loop over features
        for(j = 0; j<dataset->n_features; j++){

            printf(" > >> Feature %zu (offset = %zu): [", j, dataset->feature_offset[i * dataset->n_features + j]);

            // print spikes in feature
            for(l = 0; l<dataset->n_spikes_per_feature[i * dataset->n_features + j]-1; l++){

                printf("%zu, ", dataset->spikes[next]);
                next++;
            }
            printf("%zu]\n", dataset->spikes[next]);
            next++;
        }
    }
}