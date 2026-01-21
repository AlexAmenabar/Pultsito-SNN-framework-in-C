#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>

#include "snn_library.h"
#include "cuda/GPU_simulations.cuh"
#include "neuron_models/GPU_lif_neuron.cuh"

// functions below should be moved to a more general file

// copy SNN to GPU (usable for both single-GPU and multi-GPU)
extern "C" GPU_SNN_t** cpy_SNN2GPU(GPU_SNN_t *cpu_snn, cuda_info_t *cuda_info, simulation_configuration_t *conf){
    
    size_t N, iN, B, LT, S;
    size_t i, dev, dev_batch_size;

    N = cpu_snn->n_neurons;
    iN = cpu_snn->n_input_neurons;
    LT = cpu_snn->LT;
    S = cpu_snn->n_synapses;
    
    // allocate memory for GPU SNN strcutre for each device
    GPU_SNN_t **d_GPU_SNN = (GPU_SNN_t **)malloc(cuda_info->nDevices * sizeof(GPU_SNN_t*));

    // tmp structure to 
    GPU_SNN_t *tmp_gpu_snn = (GPU_SNN_t*)malloc(sizeof(GPU_SNN_t)); // input neurons???
    cudaError_t err;

    // loop over devices and allocate memory for all of them
    for(dev = 0; dev < cuda_info->nDevices; dev++){

        // set dev as active one
        cudaSetDevice(dev);

        // get the batch size in the device
        dev_batch_size = cuda_info->dev_batch_size[dev];

        
        /* Allocate memory */

        err = cudaMalloc(&(tmp_gpu_snn->v), N * dev_batch_size * sizeof(float)); 
        if (err != cudaSuccess) 
            printf("v allocation failed: %s\n", cudaGetErrorString(err));
        
        err = cudaMalloc(&(tmp_gpu_snn->arrI), N * dev_batch_size * sizeof(float));
        if (err != cudaSuccess) 
            printf("arrI allocation failed: %s\n", cudaGetErrorString(err)); 
        
        err = cudaMalloc(&(tmp_gpu_snn->v_thresh), N * sizeof(float)); 
        if (err != cudaSuccess) 
            printf("v_thresh allocation failed: %s\n", cudaGetErrorString(err));
        
        err = cudaMalloc(&(tmp_gpu_snn->v_rest), N * sizeof(float)); 
        if (err != cudaSuccess) 
            printf("v_rest allocation failed: %s\n", cudaGetErrorString(err));
        
        err = cudaMalloc(&(tmp_gpu_snn->r_period), N * sizeof(int)); 
        if (err != cudaSuccess) 
            printf("r_period allocation failed: %s\n", cudaGetErrorString(err));
        
        err = cudaMalloc(&(tmp_gpu_snn->r_period_remain), N * dev_batch_size * sizeof(int));
        if (err != cudaSuccess) 
            printf("r_period_remain allocation failed: %s\n", cudaGetErrorString(err)); 
        
        err = cudaMalloc(&(tmp_gpu_snn->res), N * sizeof(int)); 
        if (err != cudaSuccess) 
            printf("res allocation failed: %s\n", cudaGetErrorString(err));
        
        err = cudaMalloc(&(tmp_gpu_snn->n_neuron_input_synapses), N * sizeof(size_t));
        if (err != cudaSuccess) 
            printf("n_neuron_input_synapses allocation failed: %s\n", cudaGetErrorString(err));
        
        err = cudaMalloc(&(tmp_gpu_snn->neuron_input_synapses_offset), N * sizeof(size_t));
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));

            


        err = cudaMalloc(&(tmp_gpu_snn->w), S * dev_batch_size * sizeof(float));
        if (err != cudaSuccess) 
            printf("w allocation failed: %s\n", cudaGetErrorString(err)); 
        
        err = cudaMalloc(&(tmp_gpu_snn->init_w), S * sizeof(float));
        if (err != cudaSuccess) 
            printf("init_w allocation failed: %s\n", cudaGetErrorString(err)); 
        
        err = cudaMalloc(&(tmp_gpu_snn->dw), S * dev_batch_size * sizeof(float));
        if (err != cudaSuccess) 
            printf("dw allocation failed: %s\n", cudaGetErrorString(err)); 

        err = cudaMalloc(&(tmp_gpu_snn->acc_dw), S * sizeof(float));
        if (err != cudaSuccess) 
            printf("acc_dw allocation failed: %s\n", cudaGetErrorString(err)); 

        err = cudaMalloc(&(tmp_gpu_snn->delay), S * sizeof(int));
        if (err != cudaSuccess) 
            printf("delay allocation failed: %s\n", cudaGetErrorString(err)); 
        
        err = cudaMalloc(&(tmp_gpu_snn->lr), S * sizeof(int)); 
        if (err != cudaSuccess) 
            printf("lr allocation failed: %s\n", cudaGetErrorString(err));
        
        err = cudaMalloc(&(tmp_gpu_snn->pre_neuron_index), S * sizeof(size_t)); 
        if (err != cudaSuccess) 
            printf("pre_neuron_index allocation failed: %s\n", cudaGetErrorString(err));
        
        err = cudaMalloc(&(tmp_gpu_snn->post_neuron_index), S * sizeof(size_t)); 
        if (err != cudaSuccess) 
            printf("post_neuron_index allocation failed: %s\n", cudaGetErrorString(err));
        
        // allocate memory for spk matrix
        err = cudaMalloc(&(tmp_gpu_snn->spk_matrix), (iN + N) * LT * dev_batch_size * sizeof(char)); // allocate memory for neurons
        if (err != cudaSuccess) 
            printf("spk matrix allocation failed: %s\n", cudaGetErrorString(err));


        if(conf->learn){
        
            err = cudaMalloc(&(tmp_gpu_snn->post_fired), N * dev_batch_size * sizeof(char));
            if (err != cudaSuccess) 
                printf("post_fired allocation failed: %s\n", cudaGetErrorString(err)); 
            
            err = cudaMalloc(&(tmp_gpu_snn->post_trace), N * dev_batch_size * sizeof(float));
            if (err != cudaSuccess) 
                printf("post_trace allocation failed: %s\n", cudaGetErrorString(err)); 

            err = cudaMalloc(&(tmp_gpu_snn->pre_fired), S * dev_batch_size * sizeof(char));
            if (err != cudaSuccess) 
                printf("pre_fired allocation failed: %s\n", cudaGetErrorString(err)); 
            
            err = cudaMalloc(&(tmp_gpu_snn->pre_trace), S * dev_batch_size * sizeof(float));
            if (err != cudaSuccess) 
                printf("pre_trace allocation failed: %s\n", cudaGetErrorString(err)); 

        }


        /* Copy data from CPU to GPU */

        // cpy scalars
        tmp_gpu_snn->n_neurons = cpu_snn->n_neurons;
        tmp_gpu_snn->n_input_neurons = cpu_snn->n_input_neurons;
        tmp_gpu_snn->n_output_neurons = cpu_snn->n_output_neurons;
        tmp_gpu_snn->n_synapses = cpu_snn->n_synapses;
        tmp_gpu_snn->n_input_synapses = cpu_snn->n_input_synapses;
        tmp_gpu_snn->n_output_synapses = cpu_snn->n_output_synapses;
        tmp_gpu_snn->max_spikes = cpu_snn->max_spikes;
        tmp_gpu_snn->n_networks = cpu_snn->n_networks;
        tmp_gpu_snn->max_delay = cpu_snn->max_delay;
        tmp_gpu_snn->LT = cpu_snn->LT;
        
        // cpy not batched arrays, and then use them to initialize batched data
        cudaMemcpy(tmp_gpu_snn->v_thresh, cpu_snn->v_thresh, N * sizeof(float), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_snn->v_rest, cpu_snn->v_rest, N * sizeof(float), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_snn->r_period, cpu_snn->r_period, N * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_snn->res, cpu_snn->res, N * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_snn->n_neuron_input_synapses, cpu_snn->n_neuron_input_synapses, N * sizeof(size_t), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_snn->neuron_input_synapses_offset, cpu_snn->neuron_input_synapses_offset, N * sizeof(size_t), cudaMemcpyHostToDevice); // copy neurons information

        cudaMemcpy(tmp_gpu_snn->init_w, cpu_snn->init_w, S * sizeof(float), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_snn->delay, cpu_snn->delay, S * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_snn->lr, cpu_snn->lr, S * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_snn->pre_neuron_index, cpu_snn->pre_neuron_index, S * sizeof(size_t), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_snn->post_neuron_index, cpu_snn->post_neuron_index, S * sizeof(size_t), cudaMemcpyHostToDevice); // copy neurons information

        // final structure allocation and data copy
        err = cudaMalloc(&(d_GPU_SNN[dev]), sizeof(GPU_SNN_t));
        if (err != cudaSuccess) 
            printf("Network allocation failed: %s\n", cudaGetErrorString(err));
        cudaMemcpy(d_GPU_SNN[dev], tmp_gpu_snn, sizeof(GPU_SNN_t), cudaMemcpyHostToDevice); // copy neurons information
    }

    // free temporal GPU SNN structure
    free(tmp_gpu_snn);


    return d_GPU_SNN;
}

// copy dataset to GPU (the entire dataset)
extern "C" GPU_dataset_t** cpy_dataset2GPU(GPU_dataset_t *cpu_dataset, cuda_info_t *cuda_info){

    size_t nS, nF, nSpks;
    size_t dev;

    nS = cpu_dataset->n_samples;
    nF = cpu_dataset->n_features;
    nSpks = cpu_dataset->n_spikes;

    // allocate to temporal
    GPU_dataset_t *tmp_gpu_dataset = (GPU_dataset_t *)malloc(sizeof(GPU_dataset_t));
    GPU_dataset_t **d_gpu_dataset = (GPU_dataset_t**)malloc(cuda_info->nDevices * sizeof(GPU_dataset_t*));
    cudaError_t err;

    // loop over devices
    for(dev = 0; dev < cuda_info->nDevices; dev++){
        
        // set dev as active one
        cudaSetDevice(dev);

        
        // allocate memory in GPU
        err = cudaMalloc(&(tmp_gpu_dataset->spikes), nSpks * sizeof(size_t)); // allocate memory for neurons
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));
        
        err = cudaMalloc(&(tmp_gpu_dataset->n_spikes_per_feature), nS * nF * sizeof(size_t)); // allocate memory for neurons
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));

        err = cudaMalloc(&(tmp_gpu_dataset->sample_offset), nS * sizeof(size_t)); // allocate memory for neurons
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));

        err = cudaMalloc(&(tmp_gpu_dataset->feature_offset), nS * nF * sizeof(size_t)); // allocate memory for neurons
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));

        err = cudaMalloc(&(tmp_gpu_dataset->freq), nS * nF * sizeof(size_t)); // allocate memory for neurons
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));

        err = cudaMalloc(&(tmp_gpu_dataset->first_spk), nS * nF * sizeof(size_t)); // allocate memory for neurons
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));

        // copy data to GPU
        tmp_gpu_dataset->n_samples = cpu_dataset->n_samples;
        tmp_gpu_dataset->type = cpu_dataset->type;
        tmp_gpu_dataset->n_classes = cpu_dataset->n_classes;
        tmp_gpu_dataset->n_features = cpu_dataset->n_features;
        tmp_gpu_dataset->n_spikes = cpu_dataset->n_spikes;
        cudaMemcpy(tmp_gpu_dataset->spikes, cpu_dataset->spikes, nSpks * sizeof(size_t), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_dataset->n_spikes_per_feature, cpu_dataset->n_spikes_per_feature, nS * nF * sizeof(size_t), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_dataset->sample_offset, cpu_dataset->sample_offset, nS * sizeof(size_t), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_dataset->feature_offset, cpu_dataset->feature_offset, nS * nF * sizeof(size_t), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_dataset->freq, cpu_dataset->freq, nS * nF * sizeof(size_t), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_dataset->first_spk, cpu_dataset->first_spk, nS * nF * sizeof(size_t), cudaMemcpyHostToDevice); // copy neurons information

        // copy to GPU
        err = cudaMalloc(&(d_gpu_dataset[dev]), sizeof(GPU_dataset_t)); // allocate memory for neurons
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));

        cudaMemcpy(d_gpu_dataset[dev], tmp_gpu_dataset, sizeof(GPU_dataset_t), cudaMemcpyHostToDevice); // copy neurons information
    }

    // deallocate memory of temporal dataset
    free(tmp_gpu_dataset);


    return d_gpu_dataset;
}


// Initialize results struct in the GPU for batch_size samples
extern "C" GPU_results_t** initialize_results_str_in_GPU(size_t nDevices, size_t n, size_t s, cuda_info_t *cuda_info){

    size_t dev, dev_batch_size;
    GPU_results_t *tmp_r, **d_results;
    cudaError_t err;

    // allocate memory
    tmp_r = (GPU_results_t*)malloc(sizeof(GPU_results_t));
    d_results = (GPU_results_t**)malloc(nDevices * sizeof(GPU_results_t*));

    // loop over devices and allocate memory for results
    for(dev = 0; dev < nDevices; dev ++){
    
        // set dev as active one
        cudaSetDevice(dev);

        // batch size per dev
        dev_batch_size = cuda_info->dev_batch_size[dev];

        // allocate memory to store the number of spikes
        err = cudaMalloc(&(tmp_r->n_spks), n * dev_batch_size * sizeof(int));
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));


        // final structure
        err = cudaMalloc(&(d_results[dev]), sizeof(GPU_results_t)); // allocate memory for neurons
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));        
        
        // copy data to GPU
        cudaMemcpy(d_results[dev], tmp_r, sizeof(GPU_results_t), cudaMemcpyHostToDevice); // copy neurons information
    }
    
    // deallocate temporal struct memory
    free(tmp_r);

    return d_results;
}

extern "C" tmp_batch_cpu_t** allocate_batch_matrix_in_GPU(size_t nDevices, size_t iN, size_t T, cuda_info_t *cuda_info){

    size_t dev, dev_batch_size;
    tmp_batch_cpu_t *tmp_batch, **d_tmp_batch;
    cudaError_t err;

    // allocate memory
    tmp_batch = (tmp_batch_cpu_t*)malloc(sizeof(tmp_batch_cpu_t));
    d_tmp_batch = (tmp_batch_cpu_t**)malloc(nDevices * sizeof(tmp_batch_cpu_t*));

    // loop over devices
    for(dev = 0; dev < nDevices; dev ++){
    
        // set dev as active one
        cudaSetDevice(dev);

        // batch size per dev
        dev_batch_size = cuda_info->dev_batch_size[dev];

        // allocate memory to store the number of spikes
        err = cudaMalloc(&(tmp_batch->spikes), iN * dev_batch_size * T * sizeof(char));
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));

        // final structure
        err = cudaMalloc(&(d_tmp_batch[dev]), sizeof(tmp_batch_cpu_t)); // allocate memory for neurons
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));


        // copy data to GPU
        cudaMemcpy(d_tmp_batch[dev], tmp_batch, sizeof(tmp_batch_cpu_t), cudaMemcpyHostToDevice); // copy neurons information
    }
    
    // deallocate temporal struct memory
    free(tmp_batch);

    return d_tmp_batch;
}


// copy results from GPU to CPU
extern "C" void cpy_batch_results_GPU2CPU(GPU_results_t **cpu_results, GPU_results_t **gpu_results, cuda_info_t *cuda_info, size_t N, size_t T){

    size_t dev;
    size_t batch_size = cuda_info->batch_size;
    size_t dev_batch_size, dev_batch_offset;
    size_t nDevices = cuda_info->nDevices;

    // TODO: deallocate memory
    GPU_results_t *tmp_results = (GPU_results_t*)malloc(sizeof(GPU_results_t));
    
    // loop over devices
    for(dev = 0; dev < nDevices; dev ++){
    
        // set dev as active one
        cudaSetDevice(dev);

        // get batch size and offset of the device
        dev_batch_size = cuda_info->dev_batch_size[dev];
        dev_batch_offset = cuda_info->dev_batch_offset[dev];

        // copy from GPU to CPU
        cudaMemcpy(tmp_results, gpu_results[dev], sizeof(GPU_results_t), cudaMemcpyDeviceToHost);

        // copy results (REVISE)
        cudaMemcpy(cpu_results[dev]->n_spks, tmp_results->n_spks, N * dev_batch_size * sizeof(int), cudaMemcpyDeviceToHost);
    }

    free(tmp_results);
}


extern "C" GPU_results_t* simulate_batches_GPU(cuda_info_t *cuda_info, GPU_SNN_t *cpu_snn, GPU_dataset_t *cpu_dataset, simulation_configuration_t *conf){

    struct timespec start_gpu, end_gpu, start_helper, end_helper;
    double elpt_total = 0.0, elpt_simulation = 0.0, elpt_cpy = 0.0;

    clock_gettime(CLOCK_MONOTONIC, &start_gpu);

    // copy data (snn and dataset) to GPU
    clock_gettime(CLOCK_MONOTONIC, &start_helper);

    GPU_SNN_t **gpu_snn = cpy_SNN2GPU(cpu_snn, cuda_info, conf); // copied only once
    GPU_dataset_t **gpu_dataset = cpy_dataset2GPU(cpu_dataset, cuda_info);

    // initialize results in GPU
    GPU_results_t **gpu_results = initialize_results_str_in_GPU(cuda_info->nDevices, cpu_snn->n_neurons, cpu_dataset->n_samples, cuda_info);
    
    // initialize tmp_batch struct in gpu
    tmp_batch_cpu_t **gpu_tmp_batch = allocate_batch_matrix_in_GPU(cuda_info->nDevices, cpu_snn->n_input_neurons, conf->max_input_spikes, cuda_info); // refactorize to function

    clock_gettime(CLOCK_MONOTONIC, &end_helper);
    elpt_cpy += (end_helper.tv_sec - start_helper.tv_sec) + (end_helper.tv_nsec - start_helper.tv_nsec) / 1e9;

    // simulate
    clock_gettime(CLOCK_MONOTONIC, &start_helper);


    // TODO: adapt for several neuron types
    simulate_batches_LIF_GPU(gpu_snn, gpu_dataset, gpu_results, gpu_tmp_batch, conf, cuda_info, cpu_snn, cpu_dataset);

    /*switch (conf->neuron_type){
        // LIF neurons
        case 0:
            simulate_LIF_in_GPU(gpu_snn, gpu_dataset, gpu_results, conf, cuda_info, cpu_snn, cpu_dataset);
            break;
        default:
            simulate_LIF_in_GPU(gpu_snn, gpu_dataset, gpu_results, conf, cuda_info, cpu_snn, cpu_dataset);
            break;
    }*/
    
    clock_gettime(CLOCK_MONOTONIC, &end_helper);
    elpt_simulation += (end_helper.tv_sec - start_helper.tv_sec) + (end_helper.tv_nsec - start_helper.tv_nsec) / 1e9;
 

    // copy results to CPU again
    clock_gettime(CLOCK_MONOTONIC, &start_helper);

    GPU_results_t *cpu_results;

    clock_gettime(CLOCK_MONOTONIC, &end_helper);
    elpt_cpy += (end_helper.tv_sec - start_helper.tv_sec) + (end_helper.tv_nsec - start_helper.tv_nsec) / 1e9;

    clock_gettime(CLOCK_MONOTONIC, &end_gpu);
    elpt_total += (end_gpu.tv_sec - start_gpu.tv_sec) + (end_gpu.tv_nsec - start_gpu.tv_nsec) / 1e9;

    printf(" > Elapsed time copying data to GPU: %.3lf\n", elpt_cpy);
    printf(" > Elapsed time simulating: %.3lf\n", elpt_simulation);
    printf(" > Elapsed total time: %.3lf\n", elpt_total);

    return (cpu_results);
}