#include <cuda.h>
#include <cuda_runtime.h>

// 
#include "networks/snn.h"
#include "config/config_loader.h"


// cuda related includes
#include "cuda/cuda_networks.cuh"
#include "cuda/cuda_simulations_conf.cuh"


extern "C" GPU_SNN_t** cpy_SNN2GPU(GPU_SNN_t *cpu_snn, cuda_info_t *cuda_info, simulation_configuration_t *conf){
    
    size_t N, iN, LT, S;
    size_t dev, dev_batch_size;

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
        cuda_info->allocate_neuron_memory_cuda(tmp_gpu_snn, N, dev_batch_size);

        /*printf(" N = %zu, iN = %zu, dev_batch_size = %zu, S = %zu, LT = %zu\n", N, iN, dev_batch_size, S, LT);

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
            printf("res allocation failed: %s\n", cudaGetErrorString(err));*/


        // general neuron properties
        err = cudaMalloc(&(tmp_gpu_snn->n_neuron_input_synapses), N * sizeof(size_t));
        if (err != cudaSuccess) 
            printf("n_neuron_input_synapses allocation failed: %s\n", cudaGetErrorString(err));
        
        err = cudaMalloc(&(tmp_gpu_snn->neuron_input_synapses_offset), N * sizeof(size_t));
        if (err != cudaSuccess) 
            printf("neuron_input_synapses_offset allocation failed: %s\n", cudaGetErrorString(err));

        // synapses
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
        tmp_gpu_snn->n_networks = cpu_snn->n_networks;
        tmp_gpu_snn->LT = cpu_snn->LT;
        
        // cpy not batched arrays, and then use them to initialize batched data
        
        cuda_info->cpy_neurons_CPU2GPU_cuda(tmp_gpu_snn, cpu_snn, N);
        /*cudaMemcpy(tmp_gpu_snn->v_thresh, cpu_snn->v_thresh, N * sizeof(float), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_snn->v_rest, cpu_snn->v_rest, N * sizeof(float), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_snn->r_period, cpu_snn->r_period, N * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information
        cudaMemcpy(tmp_gpu_snn->res, cpu_snn->res, N * sizeof(int), cudaMemcpyHostToDevice); // copy neurons information*/

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